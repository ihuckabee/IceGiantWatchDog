import numpy as np
import matplotlib.pyplot as plt
import pdb
from astropy.io import fits
wno1d, fpfs1d,fs  = np.load("./ToShare/RefFiles/spec_thrutest.npy")
#fp = fpfs_base*fs #erg/cm2/s/cm
#convert to photon flux - *delta_lam * lambda / hc
h = 6.6261e-27 #erg s 
c = 3e10 #cm s-1 
#sun radius ^2 / distance to neptune ^2 
solar_r = 6.957e+8 #m
kepler_benchmark = np.load('highcadence.npy')
alb, wno = kepler_benchmark
sm = 4.514565e12 #m semi major axis of planet
rp = 2.47639e7 #m radius of planet
d_to_nep_close = 29*1.496e+11 #m 
d_to_nep_far = 31*1.496e+11 #m 
d_avg = (d_to_nep_far+d_to_nep_close)/2 #distance from observer to planet
fpfs = alb*(rp/sm)**2
rs = 6.96e8 #m radius of star
fs_obs = fs*(rs/(d_avg))**2
fp = fpfs*fs_obs
pdb.set_trace()
#fp_1d = fpfs1d*fs_obs
#fp_1d = np.array([fp_1d, fp_1d, fp_1d, fp_1d, fp_1d, fp_1d, fp_1d, fp_1d, fp_1d])
#wno1d_reshape = np.array([wno1d, wno1d, wno1d, wno1d, wno1d, wno1d, wno1d, wno1d, wno1d])
wavelength = 1e4/wno
#pdb.set_trace()
#plt.plot(1e4/wno[0], fpfs[0], label = 'from kepler')
#plt.plot(1e4/wno1d, fpfs1d, label = 'base_1d')
#plt.show()
#plt.ylabel("Flux of Planet (Fp/Fs * Fs) [erg/cm^2/s/cm] ")
#plt.xlabel("Wavelength [micron]")
#plt.show()
pdb.set_trace()
boxcar_step = np.arange(0.4,1.7, 0.1)
str_boxcar = []
bands = []
for i in range(len(boxcar_step)):
    bands.append([(boxcar_step[i] - 0.2), (boxcar_step[i] + 0.2)])
    str_boxcar.append(str(np.round(boxcar_step[i],2)))

#SPITZER
#bands = [[3.551 - 0.750, 3.551 + 0.750], [4.493-1.01, 4.493+1.01]]
spitzer_bands = np.array([3.551, 4.493])
phot = []
#NARROW RED AND BIG ASS BLUE
#bands = [[.86,.91], [.3, .7], [0.43,0.89], [0.9227, 1.1877]]
#Galileo
#bands = [[0.722, 0.732], [.7465, .7655], [0.881,0.897]]
#cassini
#VIO BL1 GRN RED CB2 CB3
#wac_center = [0.4195, 0.4597, 0.5659, 0.6469, 0.7523, 0.9387]
#wac_delta = [0.0185, 0.0531, 0.1022, 0.1388, 0.0101, 0.0088]
#wac_bands = []
#best_bands = [[0.6,1.], [.7465, .7655],[0.881,0.897]]
#str_bands = ["TESS", "galileo continuum (756nm)", "galileo strong methane (889nm)"]
best_bands = [[0.6,1.], [0.850,0.920]]
str_bands = ["TESS", "galileo strong methane (889nm)"]
#for i in range(len(wac_center)):
    #wac_bands.append([wac_center[i]-wac_delta[i]/2, wac_center[i]+wac_delta[i]/2])
#galileo
#str_bands = ["weak methane band (727nm)", "continuum band (756nm)", "strong methane band (889nm)"]
#cassini
#str_bands = ['VIO', 'BL1', 'GRN', 'RED', 'CB2', 'CB3']


for band in best_bands:
    band_index = np.where(((wavelength[0] < band[1]) & (wavelength[0] > band[0])))[0]
    #band_sum = np.sum(alb[0][band_index[0] : band_index[-1]])
    band_width = wavelength[0][band_index[0]] - wavelength[0][band_index[-1]]
    del_lam =abs(np.mean(np.diff(wavelength[0][band_index[0]:band_index[-1]+1])))
    #wno_3d[0][band_index[0]] - wno_3d[0][band_index[1]]
    #pdb.set_trace()
    phot_temp = [] 
    for i in range(len(wavelength)):
        band_sum = np.sum(fp[i][band_index[0] : band_index[-1]])
        phot_temp.append(band_sum/len(band_index))
        #phot_temp.append(np.trapz(alb[i][band_index[0] : band_index[-1]], wavelength[0][band_index[0] : band_index[-1]]))
    phot.append(phot_temp)
phot = np.array(phot)
#sorted_indices = np.argsort(order)
order_sorted = np.arange(0,len(alb))  # now monotonic
#alb_sorted = alb.T[:, sorted_indices]

del_per = []
for i in range(len(phot)):
    initial = phot[i][0]
    percent_change = (phot[i] - initial) / initial * 100
    del_per.append([order_sorted, percent_change])


hdul = fits.open(name='/Users/ihuckabee/Downloads/neptune_solsys_surfbright_001.fits', memmap=False, cache=False, lazy_load_hdus=False)
mast_data = hdul[1].data
mast_wavelength = mast_data.WAVELENGTH #angstrom
mast_flux = mast_data.FLUX #erg s-1 cm-2 A-1 
hdul_sun = fits.open(name='/Users/ihuckabee/Downloads/solar_spec.fits', memmap=False, cache=False, lazy_load_hdus=False)
mast_data_sun = hdul_sun[1].data
mast_wavelength_sun = mast_data_sun.WAVELENGTH #angstrom
mast_flux_sun = mast_data_sun.FLUX #erg s-1 cm-2 A-1 
#pdb.set_trace()

#get fp from my models using the MAST solar spectrum 

#yinterp = np.interp(xvals, x, y) xvals = more values
fs_interp = np.interp(10000*wavelength[0][::-1], mast_wavelength_sun,mast_flux_sun) #A-1

fp_model_w_MAST = fpfs1d*fs_interp
wv_model_w_MAST = 10000*wavelength[0][::-1]


#narrow red band
narrow_red_center = 0.885 #microns
narrow_red_delta = 0.025 #micron
#big ass blue band
big_blue_center = 0.5 #micron
big_blue_delta = 0.2
#kepler bandpass
kepler_center = 0.66
kepler_delta = 0.23
f105w_center = 1.0552
f105_delta = 0.1325

#GALILEO
weak_center = 0.727
cont_center = 0.756
strong_center = 0.885
weak_delta = 0.05
cont_delta = 0.19/2
strong_delta = 0.07/2

#pdb.set_trace()
#center_lam = boxcar_step*1e-4
#center_lam = np.array([0.36, 0.438, 0.545, 0.641, 0.798])*1e-4 #U, B, V, R, I 
#center_lam = np.array([narrow_red_center, big_blue_center, kepler_center,f105w_center])
#center_lam = np.array([weak_center, cont_center,strong_center])
center_lam = np.array([0.6, strong_center])
#delta_lam_boxcar = 0.1
#delta_lam = np.array([narrow_red_delta, big_blue_delta, kepler_delta, f105_delta])
#delta_lam = np.array([weak_delta, cont_delta, strong_delta])
delta_lam = np.array([0.01,strong_delta])
#delta_lam = wac_delta
#delta_lam = np.array([0.06, 0.09, 0.085, 0.15, 0.15])*1e-4
effarea = np.pi*(6.05/2)**2 #in cm^2
throughput = 0.6*0.42
eff_thru =  24#effarea*throughput# = 24 cm^2 effarea*throughput for CUTE 
eff_thru = effarea*throughput #asteria 
int_time = 60 #s

noise_level = []; snr = []; signal_flux_per_bandpass = []; noise_per_bandpass = []
for i in range(len(center_lam)):
    signal_flux = []
    signal_noise = [] 
    #mask_model = (wavelength[0][::-1] > 1e4*(center_lam[i]-delta_lam[i]/2)) & (wavelength[0][::-1] < 1e4*(center_lam[i]+delta_lam[i]/2))
    mask_model = (wavelength[0][::-1] > (center_lam[i] - delta_lam[i])) & (wavelength[0][::-1] < (center_lam[i] + delta_lam[i]))
    #mask_model_spitz_test = (wavelength[0][::-1] > 1e4*(center_lam[i] - delta_lam[i])) & (wavelength[0][::-1] < 1e4*(center_lam[i] + delta_lam[i]))
    #mask_mast_spitz_test = (mast_wavelength > 1e8*(center_lam[i] - delta_lam[i])) & (mast_wavelength < 1e8*(center_lam[i] + delta_lam[i]))
    #mask_mast = (mast_wavelength > 1e4*(boxcar_step[i]-delta_lam_boxcar)) & (mast_wavelength < 1e4*(boxcar_step[i] + delta_lam_boxcar))
    #mask_mast_v_test = (mast_wavelength > 1e4*(0.55- 0.043)) & (mast_wavelength < 1e4*(0.55 + 0.043))
    #Phi = (1.0 / (h * c)) * center_lam[i]*(np.trapz(fp[0][::-1][mask],1e-4*wavelength[0][::-1][mask])) #center_lam is in cm, so wavelength should also be in cm. yes it is.
    #Phi = (1.0 / (h * c)) * 0.000055*(np.trapz(fp[0][mask],wno[0][mask])) 
    #Phi_mast = 4*(1.0 / (h * c)) * center_lam[i]*(np.trapz(mast_flux[mask_mast],mast_wavelength[mask_mast])) 
    #Phi_model = (1.0 / (h * c)) *  center_lam[i]*(np.trapz(fp_1d[0][::-1][mask_model],1e-4*wavelength[0][::-1][mask_model])) 
    #Phi_mast_v_test = 4*(1.0 / (h * c)) * center_lam[i]*(np.trapz(mast_flux[mask_mast_spitz_test],mast_wavelength[mask_mast_spitz_test])) 
    Phi_model_v_test = (1.0 / (h * c)) *  center_lam[i]*1e-4*(np.trapz(fp[0][::-1][mask_model],1e-4*wavelength[0][::-1][mask_model])) 
    signal = Phi_model_v_test*eff_thru*int_time#*np.pi*(rp/sm)**2
    #pdb.set_trace()
    snr_signal = signal/np.sqrt(signal)
    precision = 100/snr_signal
    #signal_flux.append(signal)
    #signal_noise.append(precision) 
    snr.append(signal/np.sqrt(signal))
    noise_level.append(precision)
        #pdb.set_trace()
    #signal_flux_per_bandpass.append(signal_flux)
    #noise_per_bandpass.append(signal_noise)

#print(str_boxcar)
#print("narrow red, big blue, kepler")
print(str_bands)
print(snr)
#print(noise_level)
pdb.set_trace()
snr_diagnostic_min = []; snr_diagnostic_max = []
#plt.fill_between(del_per[0][0], y1 = np.ones_like(del_per[0][0]), y2=-np.ones_like(del_per[0][0]), color = 'k', alpha = 0.3)
from matplotlib.pyplot import cm
#color = iter(cm.rainbow(np.linspace(0, 1, len(boxcar_step))))
#color = ['purple', 'blue', 'green', 'red', 'orange', 'brown']
color = ['blue', 'red', 'red']
precision = [0.04192,0.22551]
for i in range(len(str_bands)):
    #c = next(color)
    #plt.plot(del_per[0][0], del_per[i][1], label = str_boxcar[i], color = c)
    #plt.errorbar(del_per[0][0], signal_flux_per_bandpass[i], yerr = np.array(noise_per_bandpass[i]), label = str_boxcar[i])
    #snr_diagnostic_min.append(abs(del_per[i][1]/noise_level[i])[1])
    #snr_diagnostic_max.append(abs(del_per[i][1]/noise_level[i])[5])
    #plt.plot(del_per[0][0], del_per[i][1], label = str_bands[i], color = color[i])
    plt.errorbar(del_per[0][0], del_per[i][1], yerr = precision[i], fmt='.', label = str_bands[i], color = color[i])


plt.legend()
plt.ylabel("percent change in albedo")
#plt.yscale('log')
#plt.gca().invert_yaxis()
plt.xlabel("Time [hours]")
plt.title("Change in Albedo over One Rotation")
#plt.savefig("./Figures/albedochanges_spot_zoom_dark.pdf")
plt.show()
plt.close()
pdb.set_trace()

print(str_boxcar)
#print(snr_diagnostic)
snr_diagnostic_min = np.array(snr_diagnostic_min)/(1+np.array(snr_diagnostic_min))
snr_diagnostic_max = np.array(snr_diagnostic_max)/(1+np.array(snr_diagnostic_max))
from matplotlib.pyplot import cm
#color = iter(cm.rainbow(np.linspace(0, 1, len(boxcar_step))))
edgecolor = iter(cm.bwr(np.linspace(0, 1, len(boxcar_step))))
k_lam = 2190/1000
k_del=390/1000
#h in micron
h_lam = 1630/1000
h_del = 307/1000
#j in micron
j_lam = 1220/1000
j_del = 213/1000
#for i in range(len(boxcar_step)):  
#color = ['b', 'purple', 'orange', 'r']
#color = ['purple', 'blue', 'green', 'red', 'orange', 'brown']
#bands = [[.3, .7], [0.43,0.89], [0.9227, 1.1877],[.86,.91]]
bands = best_bands
color = ['cornflowerblue', 'lightcoral']
for i in range(len(str_bands)):  
    #c = next(color)
    plt.fill_between([bands[i][0],bands[i][1]], [0, 1], color = color[i],step="pre", label = str_bands[i])#, label = str_boxcar[i])
    #plt.fill_between([np.array(boxcar_step[i]-0.1), np.array(boxcar_step[i] + 0.1)], [snr_diagnostic_min[i], snr_diagnostic_min[i]], color = c, step="pre", alpha = 0.75)#, label = str_boxcar[i])
#plt.fill_between([(h_lam-h_del), (h_lam+h_del)], [0.25, 0.25], edgecolor = 'k',linestyle = '--', step="pre", label = 'H')
#plt.fill_between([(j_lam-j_del), (j_lam+j_del)], [1, 1], edgecolor = 'k', alpha = 0.2, linestyle = '-.', step="pre", label = 'J')


#plt.legend()
plt.xlabel("Wavelength (micron)")
plt.xlim(0.4,1.1)
plt.plot(1e4/wno[0], fpfs[0]/(max(fpfs[0])), 'k-')
#plt.legend()
#plt.plot(mast_wavelength*1e-4, mast_flux*1e12, 'k', label = 'MAST Neptune Flux')
#plt.legend()
#plt.ylabel("S/N Diagnostic")
#plt.title("System Throughputs per Filter")
#plt.xlim(0.3,1.8)
plt.show()

pdb.set_trace()

'''

snr_uvis = []; snr_ir = []
#******HST STUFF*******
#HST UVIS
F625W = np.array([0.4777, 0.7707])*1e-4  # (6242.6 ± 1464.6/2) Å → microns -> cm
F475W = np.array([0.4101, 0.5445])*1e-4   # (4773.1 ± 1343.5/2) Å → microns  
F410M = np.array([0.3938, 0.4281])*1e-4   # (4109 ± 172/2) Å → microns
F467M = np.array([0.4583, 0.4783])*1e-4   # (4682.6 ± 200/2) Å → microns
F547M = np.array([0.5123, 0.5773])*1e-4  # (5447.5 ± 650/2) Å → microns
F814W = np.array([0.8039 - 0.1565/2, 0.8039 + 0.1565/2])*1e-4 
F850LP = np.array([0.91761 - 0.11925/2, 0.91761 + 0.11925/2])*1e-4 

#HST IR
F139M = np.array([1.3511, 1.4154])*1e-4  # (1383.3 ± 64.3/2) nm → microns
F160W = np.array([1.4027, 1.6710])*1e-4   # (1536.9 ± 268.3/2) nm → microns
F105W = np.array([0.9227, 1.1877])*1e-4   # (1055.2 ± 265/2) nm → microns
F125W = np.array([1.2500, 1.5345])*1e-4   # (1392.3 ± 284.5/2) nm → microns

#KEPLER
Kepler = np.array([0.43, 0.89])*1e-4
uvis_bands = [F625W,F475W, F814W, F850LP, Kepler]
ir_bands = [F139M,F160W,F105W,F125W]
lam_arr = []; del_arr = []
lam_arr_uvis = []; del_arr_uvis = []
lam_arr_ir = []; del_arr_ir = []
#pdb.set_trace()
for i in range(len(uvis_bands)):
    lam_arr_uvis.append((uvis_bands[i][0]+uvis_bands[i][1])/2)
    del_arr_uvis.append(uvis_bands[i][1]-uvis_bands[i][0])
for i in range(len(ir_bands)):
    lam_arr_ir.append((ir_bands[i][0]+ir_bands[i][1])/2)
    del_arr_ir.append(ir_bands[i][1]-ir_bands[i][0])


#UVIS OR IR ********* MUST CHOOSE   
bands = uvis_bands
phot = []

for band in bands:
    band = band*1e4
    band_index = np.where(((wavelength[0] < band[1]) & (wavelength[0] > band[0])))[0]
    #band_sum = np.sum(alb[0][band_index[0] : band_index[-1]])
    band_width = wavelength[0][band_index[0]] - wavelength[0][band_index[-1]]
    del_lam =abs(np.mean(np.diff(wavelength[0][band_index[0]:band_index[-1]+1])))
    #wno_3d[0][band_index[0]] - wno_3d[0][band_index[1]]
    #pdb.set_trace()
    phot_temp = [] 
    for i in range(len(wavelength)):
        band_sum = np.sum(fp[i][band_index[0] : band_index[-1]])
        phot_temp.append(band_sum/len(band_index))
        #phot_temp.append(np.trapz(alb[i][band_index[0] : band_index[-1]], wavelength[0][band_index[0] : band_index[-1]]))
    phot.append(phot_temp)
phot = np.array(phot)
#sorted_indices = np.argsort(order)
order_sorted = np.arange(0,len(alb))  # now monotonic
#alb_sorted = alb.T[:, sorted_indices]

del_per = []
for i in range(len(phot)):
    initial = phot[i][0]
    percent_change = (phot[i] - initial) / initial * 100
    del_per.append([order_sorted, percent_change])


uvis_thru = [0.28, 0.27, 0.23, 0.11]
ir_thru = [0.54, 0.56, 0.52, 0.56]
uvis_loc_arr = [-7327, -4645, -9858, -11181]
ir_loc_arr = [5925, 4870, 8625, 5850]

for i in range(len(lam_arr_uvis)):
    mask = (1e-4*wavelength[0][::-1] > (lam_arr_uvis[i] - del_arr_uvis[i]/2)) & (1e-4*wavelength[0][::-1] < (lam_arr_uvis[i] + del_arr_uvis[i]/2))
    Phi = (1.0 / (h * c)) * lam_arr_uvis[i]*(np.trapz(fp[0][::-1][mask],1e-4*wavelength[0][::-1][mask])) #wavelength in cm 
    #^^^^ BRO WHY IS IT NEGATIVE SOMETIMES
    signal = Phi*effarea*throughput*int_time*np.pi*(rp/sm)**2
    snr_signal = signal/np.sqrt(signal)
    precision = 100/snr_signal
    #snr.append(signal/np.sqrt(signal))
    snr_uvis.append(precision)
    
for i in range(len(lam_arr_ir)):
    mask = (1e-4*wavelength[0][::-1] > (lam_arr_ir[i] - del_arr_ir[i]/2)) & (1e-4*wavelength[0][::-1] < (lam_arr_ir[i] + del_arr_ir[i]/2))
    Phi = (1.0 / (h * c)) * lam_arr_ir[i]*(np.trapz(fp[0][::-1][mask],1e-4*wavelength[0][::-1][mask])) #wavelength in cm 
    signal = Phi*effarea*throughput*int_time*np.pi*(rp/sm)**2
    snr_signal = signal/np.sqrt(signal)
    precision = 100/snr_signal
    #snr.append(signal/np.sqrt(signal))
    snr_ir.append(precision)
    
str_bands_uvis = ['F625W','F475W','F814W', 'F850LP', 'Kepler']
str_bands_ir = ['F139M','F160W','F105W','F125W']

print(str_bands_uvis)
print(snr_uvis)
print(str_bands_ir)
print(snr_ir)

noise_level = snr_uvis
#plt.fill_between(del_per[0][0], y1 = np.ones_like(del_per[0][0]), y2=-1*np.ones_like(del_per[0][0]), color = 'k', alpha = 0.3)
#plt.plot(del_per[0][0], del_per[0][1]/noise_level[0], label = str_bands_uvis[0])
#plt.plot(del_per[0][0], del_per[3][1]/noise_level[3], label = str_bands_uvis[3])
#plt.plot(del_per[0][0], del_per[1][1]/noise_level[1], label = str_bands_uvis[1])
#plt.plot(del_per[0][0], del_per[2][1]/noise_level[2], label = str_bands_uvis[2])

plt.errorbar(del_per[0][0], del_per[0][1], yerr= noise_level[0],  capsize=3, lw = 2, label = str_bands_uvis[0])
plt.errorbar(del_per[0][0], del_per[3][1], yerr = noise_level[3],  capsize=3, lw = 2, label = str_bands_uvis[3])
plt.errorbar(del_per[0][0], del_per[1][1], yerr= noise_level[1],  capsize=3, lw = 2, label = str_bands_uvis[1])
plt.errorbar(del_per[0][0], del_per[2][1], yerr = noise_level[2],  capsize=3, lw = 2, label = str_bands_uvis[2])
plt.errorbar(del_per[0][0], del_per[4][1], yerr = noise_level[4],  capsize=3, lw = 2, label = str_bands_uvis[4])
plt.legend()
#plt.ylabel("S/N")
plt.ylabel("Percent Change in Albedo")
plt.xlabel("Time [hours]")
plt.title("Change in Albedo over One Rotation")
plt.show()
#pdb.set_trace()
'''