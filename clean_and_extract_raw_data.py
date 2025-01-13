if True:
    from scipy.io import idl
    from pdb import set_trace as stop
    from datetime import datetime
    from astropy.time import Time
    import math
    import sys
    import pickle
    import os
    import time
    import numpy as np
    import sys
    sys.path.append("path")
    import toolbox as tb
    


def load_crustal_contours():
    # Return Morschhauser crustal magnetic field model at 400 km in contour form w/ 1x1 degree resolution.
    # Output: Vector magnetic field [Br, Btheta, Bphi, Bradial, Bhorizontal]
    
    file =  '/fakeDir/Morsch_Contour_data_400_km_components.pickle'
    with open(file, "rb") as f:
       x, y, z, br, bh, info = pickle.load(f)
       
    return x, y, z, br, bh 

def get_ls_and_helio(julin,Lon1,Lat1):
    '''
    Purpose: 
        Given a location on Mars & and the date/time, find the Mars-Sun distance, solar longitude, local time, and location of the subsolar point. The code 
        borrows heavily from http://www.giss.nasa.gov/tools/mars24/help/algorithm.html (equations and constants) which is is based on Allison & McEwen (2000).
    
    Inputs: 
        julin: String Julian date
        Lon1:  Mars longitude *MUST BE WEAT LONGITUDE from 0-360 deg.
        Lat1: Mars latitude (deg.)

    Outputs: 
        Ls:         Solar longitude of Mars (deg.)
        dhelio:     Mars-Sun distance (AU)
        szacheck:   Solar zenith angle of exact location (deg.)
        dec:        Declination - latitude of subsolar point
        sublon:     Longitude of subsolar point 
        lst:        Local solar time
    '''
    Ls = []
    dhelio = []
    szacheck = []
    dec = []
    lst = []
    sublon=[]
    # Handles if non-array entered as input
    if isinstance(julin,float): 
        julin = [julin]
        Lon1 = [Lon1]
        Lat1 = [Lat1]

    # Loop through each date and calculate
    for ppp, j1 in enumerate(julin):
        deltaJul = j1-2451545.0

        # B-1 Mars Mean Anomaly
        MMA = 19.3871 + 0.52402073*deltaJul
        MMArad = math.radians(MMA)

        # B-2 Angle of fiction mean Sun
        FMS = 270.3871 + 0.524038496*deltaJul
        FMSrad = math.radians(FMS)

        # B-3 Determine Perturbers (PBS)
        Ai1 = [71,57,39,37,21,20,18]
        Ai =  [a*1e-4 for a in Ai1]
        Ti = [2.2353,2.7543,1.1177,15.7866,2.1354,2.4694,32.8493]
        Psi = [49.409,168.173,191.837,21.736,15.704,95.528,49.095]
        PBS1 =[]
        for i, a in enumerate(Ai):
            PBS1.append(a*np.cos(math.radians((0.985626*deltaJul/Ti[i]) + Psi[i])))
        PBS = sum(PBS1)

        # B-4 v-M
        vM = (10.691 + 3.0e-7*deltaJul)*np.sin(MMArad) + 0.623*np.sin(2.0*MMArad) + 0.050*np.sin(3.0*MMArad) + 0.005*np.sin(4.0*MMArad) + 0.0005*np.sin(5.0*MMArad) + PBS
        
        # B-5 Ls
        LsLs = (FMS+vM) % 360

        # C-1 Equation of time (EOT)
        EOT = 2.861*(np.sin(np.deg2rad(2.0*LsLs))) - 0.071*np.sin(np.deg2rad(4.0*LsLs)) + 0.002*np.sin(np.deg2rad(6.0*LsLs)) - vM  #EOT in degrees
        EOThrs = EOT*(1./15.0) #15 deg, per hour, multiply to get EOT in hours

        # C-2 Coordinated Mars Time ("Airy mean time")
        mmm = ((deltaJul - 4.5)/(1.0274912517) + 44796.0 - 0.0009626)*24
        MTC = mmm % 24

        # C-3 Local Mean Solar Time
        LMST = MTC - Lon1[ppp]*(1./15.)

        # C-4 Local True Solar Time (what I want)
        LTST = LMST + EOThrs
        if LTST < 0: LTST =  24. + LTST #correction for MAVEN convention????????  WHAT IS THIS!!!

        # C-5 Longitude at the subsolar point - this is West longitude as in the paper
        # "LW denotes the west longitude (measured westward from the prime meridian according to the planetary cartographic convention in the range 0±360)"
        sublon1 = ((MTC+EOThrs)*15. + 180.0) % 360

        # Hour Angle - equal to geographic longitude minus the longitude of the subsolar point
        HA = Lon1[ppp]-sublon1

        # D-1 Declination of the Sun
        decAlg = np.rad2deg(np.arcsin(0.42565*np.sin(np.deg2rad(LsLs)))) + 0.25*np.sin(np.deg2rad(LsLs))

        # Helio distance
        dhelio1 = (1.5236*(1.00436 - 0.09309*np.cos(MMArad) - 0.00436*np.cos(2.0*MMArad) - 0.00031*np.cos(3.0*MMArad)))

        # D-5 SZA - for checking
        SZAalg = np.rad2deg(np.arccos(np.sin(np.deg2rad(decAlg))*np.sin(np.deg2rad(Lat1[ppp])) + np.cos(np.deg2rad(decAlg))*np.cos(np.deg2rad(Lat1[ppp]))*np.cos(np.deg2rad(HA))))

        # Append arrays
        Ls.append(LsLs)
        dhelio.append(dhelio1)
        dec.append(decAlg)
        szacheck.append(SZAalg)
        sublon.append(sublon1)
        lst.append(LTST)
    
    # Convert to NumPys
    Ls = np.array(Ls)
    dhelio = np.array(dhelio)
    szacheck = np.array(szacheck)
    dec = np.array(dec)
    sublon = np.array(sublon)
    lst = np.array(lst)

    return Ls, dhelio, szacheck, dec, sublon, lst

def read_IDL(filename):
    ''' Read in IDL .sav file and return numpy dictionary of it '''
    idl_dat = idl.readsav(filename)
    numpy_dat = {}
    for i, idl_key in enumerate(idl_dat):
        numpy_dat[idl_key] = idl_dat[idl_key]
    
    print(numpy_dat.keys())
    return numpy_dat

def orbit_index(juldate):
    orbit = [0]
    index = 1
    for i, j1 in enumerate(juldate):
        if j1 == juldate[-1]: 
            break
        if (juldate[i+1] - juldate[i])*24 > 0.5:
            index = index + 1
        orbit.append(index)
    return np.array(orbit)



##################################################################
# Get list of filenames after fixing problem, then sort.
# There was a problem matching dates later in the program because 
# these filenames are not sorted by date, and they have 
# inconsistent "0" pads in the months.
##################################################################

# Data paths
filenames_dir = 'sample_data/'     # Location of raw data
new_filename_dir = 'cleaned_data/' # To save cleaned data files

# Filenames and dates lists
new_filenames = []         # Cleaned filename
original_filenames = []    # Original filename
jd_file =[]                # Juldate of filename

# Make dates consistent & fix "0" padding problem
for file1 in os.listdir(filenames_dir):
    if file1.endswith(".sav"):
        
        # if month has no leading zero, add one - this is necessary for sorting files
        mm = file1.split("-")[1]
        if len(mm) == 1:
            file2 = file1.replace("-" + mm + "-", "-" + "0" + mm + "-")
        else:
            file2 = file1

        original_filenames.append(file1)
        new_filenames.append(file2)
        jd_file.append(Time(file1.split('_')[-1].split('.')[0]).jd)

filenames, newfilenames, jd_file = np.array(original_filenames), np.array(new_filenames), np.array(jd_file)

# sort filenames
filenames = filenames[np.argsort(jd_file)]
filenames = [filenames_dir + f1 for f1 in filenames]



##############################################################
# Loop through each file and extract relevant data into
# numpy arrays
##############################################################

# Lists to append w/ extracted data
bg = []         # Background noise level (eV/eV/cm2/s/sr) 
angspecinc = [] # Incoming particles: angle-averaged differential energy flux (eV/eV/cm2/s/sr)
angspecref = [] # Reflected particles: angle-averaged differential energy flux (eV/eV/cm2/s/sr)
angspectot = [] # All particles: angle-averaged differential energy flux (eV/eV/cm2/s/sr)
fluxtot = []    # Total flux: angspectot intergrated over all angles      (eV/eV/cm2/s)
time = []       # UNIX time of each observation (byte strings)
alt = []        # Altitude above Mars (km)
sza = []        # Solar zenith angle (deg.)
elon = []       # East longitude on Mars (deg.)
lat = []        # Latitude on Mars (deg.)
npen = []       # Density flux of incoming particles (#/eV/cm2/s)
nbpen = []      # Density flux of reflected particles (#/eV/cm2/s)


# Loop through each file
for checkerind, file in enumerate(filenames):

    dat = tb.read_IDL(file) 
    print(file)

    time.extend(dat['time1'])
    sza.extend(dat['sza1'])
    alt.extend(dat['alt1'])
    elon.extend(dat['lon1'])
    lat.extend(dat['lat1'])
    npen.extend(dat['npen1'])
    nbpen.extend(dat['nbpen1'])
    fluxtot.extend(dat['fluxtot1'])
    bg.extend(dat['background'])

    # Energy spectra - there are 48 energies in each spectrum
    angspecinc.extend(dat['angspecinc1'].T)
    angspecref.extend(dat['angspecref1'].T)
    angspectot.extend(dat['angspectot1'].T)

# Decode dates: byte strings --> string, datetime, and juldate
datetimes = np.array([datetime.fromtimestamp(t1) for t1 in time])
juldate = Time(datetimes).jd

# Convert to numpys
bg, juldate, datetimes, alt, sza, elon, lat, npen, nbpen = np.array(bg), np.array(juldate), np.array(datetimes), np.array(alt), np.array(sza), np.array(elon), np.array(lat), np.array(npen), np.array(nbpen)

# Get crustal field strength at 400 km using 1x1 deg. resolution crustal field map
Bcrust = tb.load_crustal_contours(alt=400)

# Get solar longitude (Ls)
Ls = get_ls_and_helio(juldate,elon,lat)[0]

# Energy bins for each measurement - Always the same
energy1 = dat['energy1'] 
energy = energy1.T[0]

# Make "orbit #" index (Not! the official MAVEN orbit #)
orbit = orbit_index(juldate)

# Save dataset for quick loading for analysis
with open(new_filename_dir + 'sample_data_cleaned.pickle', "wb") as f:  
    pickle.dump(
        (bg, 
         energy, 
         angspecinc, 
         angspecref, 
         angspectot, 
         fluxtot, 
         juldate, 
         datetimes, 
         alt, 
         sza, 
         elon, 
         lat, 
         npen, 
         nbpen, 
         Ls, 
         Bcrust
         ), f)
stop()













