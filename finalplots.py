import numpy as np
import matplotlib.pyplot as plt
import pdb
from astropy.io import fits
wno1d, fpfs1d,fs  = np.load("./ToShare/RefFiles/spec_thrutest.npy")
h = 6.6261e-27 #erg s 
c = 3e10 #cm s-1 
brightclouds = np.load('highcadence_bright.npy')
#brightclouds = np.load('longdurhighcadence_darkonly.npy')
alb, wno = brightclouds
sm = 4.514565e12 #m semi major axis of planet
rp = 2.47639e7 #m radius of planet
d_to_nep_close = 29*1.496e+11 #m 
d_to_nep_far = 31*1.496e+11 #m 
d_avg = (d_to_nep_far+d_to_nep_close)/2 #distance from observer to planet
fpfs_bright = alb*(rp/sm)**2
rs = 6.96e8 #m radius of star
fs_obs = fs*(rs/(d_avg))**2
fp_bright = fpfs_bright*fs_obs

darkclouds = np.load('highcadence_darkonly.npy')
#darkclouds = np.load('longdurhighcadence_darkonly.npy')
alb, wno = darkclouds
fpfs_dark = alb*(rp/sm)**2
fp_dark = fpfs_dark*fs_obs

wavelength = 1e4/wno
fp = [fp_bright, fp_dark]
#galileo -> [.7465, .7655]
#TESS -> [0.6,1.]
#best_bands = [[0.6,1.], [0.850,0.920]]
best_bands = [[0.850,0.920]]
#best_bands = [[0.6,1.]]
#str_bands = ["TESS", "galileo strong methane (889nm)"]
str_bands = ["galileo strong methane (889nm)"]
color = [['lightcoral', 'cornflowerblue'], ['darkred', 'darkblue']]
#color = [['cornflowerblue', 'lightcoral'], ['darkblue', 'darkred']]
plt.figure(figsize=(10, 5))
for k in range(len(fp)):
    phot = []
    for band in best_bands:
        band_index = np.where(((wavelength[0] < band[1]) & (wavelength[0] > band[0])))[0]
        band_width = wavelength[0][band_index[0]] - wavelength[0][band_index[-1]]
        del_lam =abs(np.mean(np.diff(wavelength[0][band_index[0]:band_index[-1]+1])))
        phot_temp = [] 
        for i in range(len(wavelength)):
            band_sum = np.sum(fp[k][i][band_index[0] : band_index[-1]])
            phot_temp.append(band_sum/len(band_index))
        phot.append(phot_temp)
    phot = np.array(phot)
    order_sorted = np.arange(0,len(alb))
    del_per = []
    for i in range(len(phot)):
        mean = np.mean(phot[i])
        percent_change = (phot[i] - mean) / mean * 100
        del_per.append([order_sorted, percent_change])
    #pdb.set_trace()
    precision_cute = [0.08683402536, 0.636739] #TESS then strong methane
    precision_asteria = [100/632.7105, 100/86.298]
    plt.axhline(y=0.0, color='k', linestyle='--', alpha = 0.5)
    for i in range(len(str_bands)):
        if k == 0:
            plt.errorbar(del_per[0][0], del_per[i][1], yerr = precision_cute[i], fmt = 'o', label = str_bands[i], color = color[k][-1]) #*(16.11/32)
        else:
             plt.errorbar(del_per[0][0], del_per[i][1],yerr = precision_cute[i], fmt = 's', color = color[k][-1])

plt.ylabel("% change in albedo")
plt.xlabel("Time [hours]")
plt.title("Change in Albedo over One Rotation")
#plt.ylim(-7,4)
#plt.ylim(-2,2)

plt.savefig('TESS.pdf')
plt.show()
plt.close()
