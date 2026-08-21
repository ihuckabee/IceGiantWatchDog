#Plotting given the 3D and 1D model spectra
import matplotlib.pyplot as plt
import pickle
import numpy as np
import os
os.environ["picaso_refdata"]="/Users/ihuckabee/Documents/picaso-master/reference"
os.environ['PYSYN_CDBS']="/Users/ihuckabee/Documents/picaso-master/grp/redcat/trds"
import sys
sys.path.append("..")
refdata = os.getenv("picaso_refdata")
from picaso import justdoit as pdi

#read in observed HST/SpeX spectrum, given in microns and observed reflectivity, the "Combined HST/STIS and IRTF/SpeX central-meridian-averaged I/F spectra" 
f = open("./RefFiles/obsrefspec.txt",'r')
lines = f.readlines()[110:] #data doesn't start til 0.3 microns
wv_obs = []; alb_obs = []
for i in range(len(lines)):
    wv_obs.append(float(lines[i].split()[0]))
    alb_obs.append(float(lines[i].split()[1]))

#read in 3D spectrum output, native resolution
sm = 4.514565e12 #m
rp = 2.47639e7 #m
fnames = ['./RefFiles/2x2_cldv3.pic']
wno_3d = []; alb_3d = []
for i in range(len(fnames)):
        infile = open(fnames[i], 'rb')
        out3d = pickle.load(infile)
        #pdb.set_trace()
        wno,fpfs = out3d[0],out3d[1]
        alb = fpfs*(sm/rp)**2
        wno_3d.append(wno)
        alb_3d.append(alb)

#read in 1D spectrum outputs, native resolution
wno_1d, fpfs_1d = np.load("./RefFiles/1D_cloudbasev3.npy")
wno_free, fpfs_free = np.load("./RefFiles/cloudfreespec.npy")
#wno_1dcld, fpfs_1dcld = np.load("./RefFiles/1Dcloudbase.npy")
#wno_1dtest, fpfs_1dtest = np.load("./RefFiles/1Dcloud_test.npy")
sm = 4.514565e12 #m
rp = 2.47639e7 #m
alb_1d = fpfs_1d*(sm/rp)**2
alb_free = fpfs_free*(sm/rp)**2


#plt.plot(1e4/wno_free,alb_free,label="1D Cloud Free")
#plt.plot(1e4/wno_1dtest,alb_1dtest,label="1D Cloud Base Test", alpha = 0.75)
#plt.plot(1e4/wno_1d,alb_1d,label="1D Cloud Base")
#plt.plot(wv_obs,alb_obs, 'k--', label="HST/SpeX Data", alpha = 0.75)
#plt.plot(1e4/wno_3d[0],alb_3d[0], label="10x10 Cloud Base", alpha = 0.5) #something is WRONG with this for sure
#plt.plot(1e4/wno_3d[1],alb_3d[1], label="10x10 Cloud Free", alpha = 0.5)
#plt.xlabel("Wavelength (microns)")
#plt.ylabel("Albedo")
#plt.legend()
#plt.xlim(0.3,1.6)
#plt.show()

wv_obs, alb_obs = pdi.mean_regrid(wv_obs, alb_obs, R=100)
wno_1d, alb_1d = pdi.mean_regrid(wno_1d, alb_1d, R=100)
wno_3d_v3, alb_3d_v3 = pdi.mean_regrid(wno_3d[0], alb_3d[0], R=100)
#wno_3d_v2, alb_3d_v2 = pdi.mean_regrid(wno_3d[1], alb_3d[1], R=100)


#residual_cf = alb_free - alb_obs
residual_1cb = alb_1d[::-1][:-1] - alb_obs
#residual_3cb_v2 = alb_3d_v2[::-1][:-1] - alb_obs
residual_3cb_v3 = alb_3d_v3[::-1][:-1] - alb_obs
# Create the figure and two subplots (vertically stacked)
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 6), sharex=True, 
                               gridspec_kw={'height_ratios': [3, 1]})

# Top plot: Albedo vs Wavelength
#ax1.plot(1e4/wno_1d, alb_free, label='1D Cloud Free')
ax1.plot(1e4/wno_1d[::-1][:-1], alb_1d[::-1][:-1], label='1D Cloud Base')
#ax1.plot(1e4/wno_3d_v2[::-1][:-1], alb_3d_v2[::-1][:-1], label='3D Cloud Base v2')
ax1.plot(1e4/wno_3d_v3[::-1][:-1], alb_3d_v3[::-1][:-1], label='3D Cloud Base v3')
ax1.plot(wv_obs, alb_obs, 'k--', label='HST/SpeX Data')
ax1.set_ylabel("Albedo")
ax1.legend()
ax1.grid(True)

# Bottom plot: Residuals
#ax2.plot(1e4/wno_free, residual_cf, label='Cloud Free - Data', color='blue')
ax2.plot(1e4/wno_1d[::-1][:-1], residual_1cb, label='1D- Data')
#ax2.plot(1e4/wno_3d_v2[::-1][:-1], residual_3cb_v2, label='3D - Data v2')
ax2.plot(1e4/wno_3d_v3[::-1][:-1], residual_3cb_v3, label='3D - Data v3')
ax2.axhline(0, color='k', linestyle='--', linewidth=1)
ax2.set_xlabel("Wavelength (microns)")
ax2.set_ylabel("Residual")
ax2.legend()
ax2.grid(True)

plt.tight_layout()
plt.show()