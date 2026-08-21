import numpy as np
import matplotlib.pyplot as plt
import pdb
from astropy.io import fits

wno1d, fpfs1d, fs = np.load("./ToShare/RefFiles/spec_thrutest.npy")
h = 6.6261e-27  # erg s
c = 3e10  # cm s-1
longdur_clouds = np.load("longdurhighcadence_both_updatedwind.npy")
alb, wno = longdur_clouds
sm = 4.514565e12  # m semi major axis of planet
rp = 2.47639e7  # m radius of planet
d_to_nep_close = 29 * 1.496e11  # m
d_to_nep_far = 31 * 1.496e11  # m
d_avg = (d_to_nep_far + d_to_nep_close) / 2  # distance from observer to planet
fpfs_cloud = alb * (rp / sm) ** 2
rs = 6.96e8  # m radius of star
fs_obs = fs * (rs / (d_avg)) ** 2
fp_clouds = fpfs_cloud * fs_obs
wavelength = 1e4 / wno

#best_bands = [[0.6, 1.0], [0.850, 0.920]]
#best_bands = [[0.850,0.920]]
best_bands = [[0.6,1.]]
#best_bands = [[0.42,.9]]
str_bands = ["TESS", "galileo strong methane (889nm)"]
str_bands = ["galileo strong methane (889nm)"]
color = [["lightcoral", "cornflowerblue"], ["darkred", "darkblue"]]

phot = []
for band in best_bands:
    band_index = np.where(((wavelength[0] < band[1]) & (wavelength[0] > band[0])))[0]
    band_width = wavelength[0][band_index[0]] - wavelength[0][band_index[-1]]
    del_lam = abs(np.mean(np.diff(wavelength[0][band_index[0] : band_index[-1] + 1])))
    phot_temp = []
    for i in range(len(wavelength)):
        band_sum = np.sum(fp_clouds[i][band_index[0] : band_index[-1]])
        phot_temp.append(band_sum / len(band_index))
    phot.append(phot_temp)
phot = np.array(phot)
order_sorted = np.arange(0, len(alb))
del_per = []
for i in range(len(phot)):
    mean = np.mean(phot[i])
    percent_change = (phot[i] - mean) / mean * 100
    del_per.append([order_sorted, percent_change])
# pdb.set_trace()
precision_cute = [0.08683402536, 0.636739]  # TESS then strong methane
precision_asteria = [100 / 632.7105, 100 / 86.298]
#import pdb; pdb.set_trace()
rot_images_x = [del_per[0][0][3], del_per[0][0][19], del_per[0][0][72], del_per[0][0][105], del_per[0][0][123],del_per[0][0][171], del_per[0][0][210]]
rot_images_y = [del_per[i][1][3], del_per[i][1][19], del_per[i][1][72], del_per[i][1][105], del_per[i][1][123], del_per[i][1][171], del_per[i][1][210]]
plt.figure(figsize=(12,5))
plt.axhline(y=0.0, color="k", linestyle="--", alpha=0.5)

for i in range(len(str_bands)):
    del_per_x = del_per[0][0][::3]#[x for k, x in enumerate(del_per[0][0]) if (k + 1) % 3 != 0]

    del_per_y = del_per[i][1][::3]#[x for k, x in enumerate(del_per[i][1]) if (k + 1) % 3 != 0]

    plt.errorbar(
        del_per_x,
        del_per_y,
        #del_per[i][1],
        #yerr=precision_asteria[i],
        fmt = 'o-',
        label=str_bands[i],
        color='darkred',
    )  # *(16.11/32)
    plt.fill_between(
        del_per_x,
        np.array(del_per_y)-precision_cute[0], 
        y2 = np.array(del_per_y)+precision_cute[0],
        color = 'darkred',
        alpha = 0.5
    )

plt.plot(rot_images_x, rot_images_y, 'ro', markersize = 8, mec='k', zorder = 5)
plt.ylabel("% change in albedo")
plt.xlabel("Time [hours]")
#plt.title("Change in Albedo over One Rotation")
# plt.ylim(-7,4)
# plt.ylim(-2,2)

plt.savefig("StrongMethane_Long.pdf")
plt.show()
plt.close()
