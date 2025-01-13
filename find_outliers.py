if True:
    from pdb import set_trace as stop
    from astropy.time import Time
    import pickle
    import os
    import time
    import numpy as np

##################################
########## Algorithm  ###########
##################################

# The data used in this program come from # get_swia_data_penh.pro
#-####################################################################-
# These contain all SWIA spectra below 200 km
# Reloading new archive spectra in March 2022
# Reload all IDL .sav files or quickload already-stored numpy arrays?
#-####################################################################-


######################    ALGORITHM    #####################################
###########################################################################
# This alg is updated to April 2022
# Algorithm to find orbits with evidence of direct solar wind penetration
# 1. Calculate orbit-avg spectrum below 200 km
# 2. If peak flux is unusually large then keep
# 3. If main peak is within SW energies then keep
# 4. Other thresholds used such as requiring a relatively large second He++ peak (look at cose below for details)
# 5. Look by eye and keep only those with H+,. He++, and He+ peaks
###########################################################################
###########################################################################
    
    
# Load Data
#swangspectot = np.load('detections_angspectot_REDO.npy',allow_pickle=True) #load detection angspectot
if True:
    if small_data:
        with open('small_file.pickle', "rb") as f: energy, juldate, alt, sza, elon, lat, Ls, Bcrust = pickle.load(f)
        print('********Only a subset of data is loaded - no dates or spectra!!!!!')
    
    else:
        with open('FILENAME_HERE', "rb") as f:bg, energy, angspecinc, angspecref, angspectot,fluxtot, juldate, datetimes, alt, sza, elon, lat, npen, nbpen, vpen, mmax, Ls, Bcrust = pickle.load(f)
        angspectot = np.array(angspectot)
        print('Full data loaded')


#quick load orbit indexx
with open('fake_orbit_numbers_REDO.pickle', "rb") as f: orbit = pickle.load(f)

    
# Solar wind flag array - 0 if no SW detection, 1 if solar wind detection
solarWindFlag = [] # indices of detections
solarWindFlagTime = []
solarWindFlagOrbit = [] 
proms = [] # store prominences of detections

#### Algorithm Loop to find cases ###############
if True:
    # Loop through each orbit
    for orbitPick in np.unique(orbit):
        # skip first orbit there is only one spectrum
        if orbitPick == 0: continue

        #if orbitPick>4000: break

        #Orbit 3975 is the Crismani paper example on Oct 10 2017
        
        #orbitPick = 513
        #orbitPick = 3975
        orbitPick = 3861
        

        # Pick an orbit based on date
        #thedate = '2016-11-30T7:00'
        #orbitPick = orbit[np.argmin(abs(juldate - Time(thedate).jd))]

        # Pick orbit
        pp = np.where(orbit==orbitPick)[0]

        # Only want dayside data
        #if np.max(sza[pp]) > 80:
        #    continue

        # Calculate orbit-averaged spectrum
        avgspec1 = np.mean(angspectot[pp].T,axis=1)


        # plot all spectra from orbit
        #plot_all_spectra(energy,angspectot[pp],dp[pp])

        # plot heat map spectra
        if False:
            x = mdates.date2num(dp[pp])
            y = energy
            z = np.log10(angspectot[pp])
            X,Y = np.meshgrid(x,y)
            fig, ax = plt.subplots()
            ax.pcolormesh(X.T,Y.T,z)
            ax.set_yscale('log')
            ax.xaxis_date()
            # We can use a DateFormatter to choose how this datetime string will look.
            # I have chosen HH:MM:SS though you could add DD/MM/YY if you had data
            # over different days.
            date_format = mdates.DateFormatter('%H:%M:%S')
            ax.xaxis.set_major_formatter(date_format)
            # This simply sets the x-axis data to diagonal so it fits better.
            # fig.autofmt_xdate()
        
        # Find peak in orbit-averaged spectrum
        # Only look above 100 ev for sc potential reasons
        kk = np.where(energy > 100)[0]


        # Update March 2022 - using archive data (when avilable) so skipping this step for now
        # Since this is coarse data two energy bins have the same flux so
        # Rebin spectrum by combing adjacent bins - only needed for coarse data which I am no longer using. 
        if False:
            newbins, newflux = [], []
            ind = 0
            for i in range(len(energy)):
                
                if energy[ind+1] == energy[-1]: break

                # Geometric mean for energy bins
                newbins.append(np.sqrt(energy[ind]*energy[ind+1]))

                # Mean for fluxes
                newflux.append(np.mean([avgspec1[ind],avgspec1[ind+1]]))

                ind = ind + 2
        newbins = energy
        newflux = avgspec1
        newbins, newflux = np.array(newbins), np.array(newflux)
        

        #################################################################################################################################
        #################################################################################################################################
        # Instead of looking for peaks using average spectrum from entire orbit below 200 km, want to window average over some time length
        # and flag data that has two potential peaks. This will be able to extend to highewr altitudes in the future.
        #################################################################################################################################
        #################################################################################################################################

        # Edges of time in orbit
        start, stop1 = juldate[pp[0]], juldate[pp[-1]]

        # Define window width
        width_minutes = 1.0 # width in minutes
        window_width = width_minutes/60./24. # width in days (to work with juldates)

        # Make windows
        windows = np.arange(start,stop1 + window_width ,window_width)
        

        # Loop through windows
        for wi, window1 in enumerate(windows):

            #print(wi)
            # break if last window
            if window1 == windows[-1]: continue

            # indices of time window
            windowInd = np.where( (orbit==orbitPick) & (juldate >= windows[wi]) & (juldate < windows[wi+1]) )[0]

            # Calculate average spectrum in time window - NORMALIZED TO MAX VALUE
            avgspecWindow = np.mean(angspectot[windowInd].T,axis=1) #change windowInd to pp to use full orbit-averaged specrtum like in the paper
            #avgspecWindow = avgspecWindow/max(avgspecWindow)
            
            
            
            # Find releative maxima in spectra. prominence is how high the peak is relative tot he "baseline" of the signal (lowest values of spectra)
            # Threshold is minimum required amount the peak is above adjacent data points
            peaks, properties = find_peaks(avgspecWindow, prominence=0)#,threshold=1e4) # originally had no threshold keyword
            prominences = properties['prominences']
            
            # Get rid of anything below 150 eV
            prominences  = prominences[peaks < 35]
            peaks = peaks[peaks < 35]

            # If no peaks below 150 eV then continue
            if len(peaks)  < 1:
                print('No peaks above above 150 ev')
                continue
        
            # sort peaks by prominence
            sortind = np.argsort(prominences) 
    
            # Is the most prominent peak within the typical solar wind energy range? -- *** THIS NOW
            # BOUNDs THE MAX AND MIN OF OBSERVED SW PARAMS AT MARS. So it only rejects peaks not in the solar wind range
            most_prom_energy = energy[peaks[np.argmax(prominences)]]

            # Pick which energy range to keep. By inspection of al Mars observations, low threshold of 300 eV captures everything, high threshold
            # of  500 ev captures everything
            # If most prominent energy is not in typical sw energy, then reject
            if (most_prom_energy < 300) or (most_prom_energy > 5000):
                print(str(orbitPick) + '  Most prominent peak is not within typical solar wind energy range')
                #solarWindFlag.append(0)
                #solarWindFlagTime.append(windows[wi]) #juldate of start of window
                continue
        
            # If most prominent peak is not above some differential flux threshold, reject
            # likely a noisy spectrum with no distinct H+ peak. And is just low counts in general
            # The noise level is on the order of 5e3 eflux, so if not a promient enough peak,
            # you get a bunch of false detections from just noise 
            if True:
                if (avgspecWindow[peaks[np.argmax(prominences)]] < 5e4):
                    #print(str(orbitPick) + '  Most prominent peak is not above flux threshold')
                    #solarWindFlag.append(0)
                    #solarWindFlagTime.append(windows[wi]) #juldate of start of window
                    continue

            # If the most prominent peak is not prominent enough then reject. If I don't do this
            # then some false detections arise that don't have any prominent H+ peak at all.
            if max(prominences) < 5e4: #1e4
                print(str(orbitPick) + '  No prominent peaks')
                continue

            # Is there a secondary peak at 2x the sw energy (the solar wind He++ alphas)?
            # All energy bins in SWIA are separated as E1/E2  = 1.16 where E1 and E2 are adjacent energy bins
            Esort = np.flip(newbins[peaks[sortind]]) # Energies of most prominent peaks listed in order, first is the most prominent
            Esort = Esort[Esort > 300] # remove lower energies that cant be alphas
            Esort = Esort[Esort < 7000] # Remove higher energies that aren't alphas for typical solar wind energies

            # Now check if in the remaining peaks, any of them are ~2x the solar wind energy
            twoTimesCheck = np.array([e1/Esort[0] for e1 in Esort]) # divide each prominant peak energy by the most prominent peak (which is the SW energy). 
            # fluxes of peaks separated by 2x in energy
            peakFluxes = avgspecWindow[[np.where(energy==e1)[0][0] for e1 in Esort]]
            switch = np.any((twoTimesCheck > 1.5) & (twoTimesCheck < 2.5 )) #check if any peaks close to 2x SW energy (most prominent peak) - if this is True then we might have a detection
            if switch: switchind = np.where((twoTimesCheck > 1.5) & (twoTimesCheck < 2.5 ))[0]
            
        
            # If no peak at 2x solar wind energy then reject
            if not switch:
                print(str(orbitPick) + '  No Peak at 2x solar wind energy')
                #solarWindFlag.append(0)
                #solarWindFlagTime.append(windows[wi]) #juldate of start of window
                continue
            
            # If the peak at 2x the SW energy is small (not priminent), then reject
            # Need to adjust this but for now using ratio of prominences peak flux
            # Don't like using ratio of prominences due to how the algorithm works
            # Better is to compare ratio of peak most prominent flux/secondary peak flux - I disagree now., need to chekc if secondary peak is prominent enough or get a ton of false detections
            indsort  = [np.where(energy[peaks] == e1)[0][0] for e1 in Esort] 
            # Note some times there are more than 2 peaks in the range because I use >1.5 <2.5. But since the peakFluxes array is in order of prominences, loop through all of them (usually not more than 2)
            promtest = [] #ratio of 2x energy prominence and SW peak prominence
            fluxtest = [] # even if promtest doesnt pass, if second peak > 1e5 flux then keep as detection
            for sss1 in switchind:
                promtest.append(prominences[indsort[0]]/prominences[indsort[sss1]]) # ration of prominences for (most prominent peak)/(prominence of secondayr peak)
                kkkk = np.where(energy == Esort[sss1])[0][0]
                fluxtest.append(avgspecWindow[kkkk])
                
            # Another check for prominence of second peak by comparing peak flux to surrounding pouints
            if max([avgspecWindow[kkkk]/avgspecWindow[kkkk+2],avgspecWindow[kkkk]/avgspecWindow[kkkk+2]]) < 2:
                print('Secondary Peak not large ENOUGH!!!')
                continue

            # Test prominence of second peak relative to main peak or if above 1e5
            # If passes this test then it is a detection - store it as one
            print(promtest,fluxtest)
            
            if (np.max(promtest) <= 5) or (np.max(fluxtest) > 1e5):  
                print(str(orbitPick) + ' DETECTION!')
                solarWindFlag.append(windowInd)
                solarWindFlagTime.append(windows[wi]) #juldate of start of window
                solarWindFlagOrbit.append(orbitPick)
                proms.append((prominences[indsort[0]],prominences[indsort[1]]))
                
                
                stop()
                # Plot for inspection
                plt.figure()
                plt.loglog(energy,avgspecWindow)
                plt.loglog(energy,avgspec1,ls='--',zorder=-32,alpha=0.5)
                plt.axvline(Esort[0],color='k') # mark most prominent peak
                plt.axvline(energy[kkkk],color='r') # mark secondary peak @ 2x energy
                plt.title(Time(windows[wi],format='jd').datetime)
                plt.text(300,1e5,'Orbit ' + str(orbitPick))
                plt.xlim(1e2,1e4)
                plt.ylim(1e3,1e7)
                #stop()
                
    # store detection data
    swflag = np.array(solarWindFlag)
    swtime = np.array(solarWindFlagTime) #juldate of start of window
    sworbit = np.array(solarWindFlagOrbit)
    orbs = np.unique(sworbit)
    swjd = []
    for orb1 in orbs:
        pp = np.where(orbit == orb1)[0]
        swjd.append(juldate[pp[0]])

    #with open('detections_list_REDO_Aug_2023.pickle', "wb") as f:  pickle.dump((swflag,swtime,sworbit,orbs,swjd), f)
    stop()
    
    # plot all algorithm detections - first make folder name
    foldername = os.getcwd() + '/Detections/' + time.strftime('%Hh%Mm%Ss_%m-%d-%Y')
    os.makedirs(foldername)
    if True:
        for orb1 in orbs:
            plot_time_series(orb1,foldername +'/')
