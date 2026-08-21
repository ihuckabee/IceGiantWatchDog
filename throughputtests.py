import numpy as np
import matplotlib.pyplot as plt
import pdb
wno_base, fpfs_base, fs = np.load("./ToShare/RefFiles/spec_thrutest.npy")
sm = 4.514565e12 #m
rp = 2.47639e7 #m
fp = fpfs_base*fs #erg/cm2/s/cm
#convert to photon flux - *delta_lam * lambda / hc
h = 6.6261e-27 #cm2 g s-1 
c = 3e10 #cm s-1 
wavelength = 1e4/wno_base
#HST UVIS
F625W = np.array([0.4777, 0.7707])*1e-4  # (6242.6 ± 1464.6/2) Å → microns -> cm
F475W = np.array([0.4101, 0.5445])*1e-4   # (4773.1 ± 1343.5/2) Å → microns  
F410M = np.array([0.3938, 0.4281])*1e-4   # (4109 ± 172/2) Å → microns
F467M = np.array([0.4583, 0.4783])*1e-4   # (4682.6 ± 200/2) Å → microns
F547M = np.array([0.5123, 0.5773])*1e-4  # (5447.5 ± 650/2) Å → microns
F814W = np.array([0.8039 - 0.1565, 0.8039 + 0.1565])*1e-4 
F850LP = np.array([0.91761 - 0.11925, 0.91761 + 0.11925])*1e-4 

#HST IR
F139M = np.array([1.3511, 1.4154])*1e-4  # (1383.3 ± 64.3/2) nm → microns
F160W = np.array([1.4027, 1.6710])*1e-4   # (1536.9 ± 268.3/2) nm → microns
F105W = np.array([0.9227, 1.1877])*1e-4   # (1055.2 ± 265/2) nm → microns
F125W = np.array([1.2500, 1.5345])*1e-4   # (1392.3 ± 284.5/2) nm → microns

#uvis_bands = [F625W,F475W,F410M,F467M,F547M, F814W, F850LP]
uvis_bands = [F625W,F475W, F814W, F850LP]
ir_bands = [F139M,F160W,F105W,F125W]
lam_arr = []; del_arr = []
lam_arr_uvis = []; del_arr_uvis = []
lam_arr_ir = []; del_arr_ir = []
#lam_arr
#[6.242e-08, 4.773e-08, 4.1095e-08, 4.6829999999999996e-08, 5.4479999999999996e-08, 1.38325e-07, 1.5368500000000002e-07, 1.0552e-07, 1.3922499999999998e-07]
for i in range(len(uvis_bands)):
    lam_arr_uvis.append((uvis_bands[i][0]+uvis_bands[i][1])/2)
    del_arr_uvis.append(uvis_bands[i][1]-uvis_bands[i][0])
for i in range(len(ir_bands)):
    lam_arr_ir.append((ir_bands[i][0]+ir_bands[i][1])/2)
    del_arr_ir.append(ir_bands[i][1]-ir_bands[i][0])
#uvis_thru = [0.28, 0.27, 0.22,0.28,0.27, 0.23, 0.11]
#uvis_thru = [0.28, 0.27, 0.23, 0.11]
uvis_thru = [0.28, 0.27, 0.23, 0.11]
ir_thru = [0.54, 0.56, 0.52, 0.56]

#uvis_loc_arr = [-7327, -4645, -3148, -4450, -5966, -9858, -11181]
uvis_loc_arr = [-7327, -4645, -9858, -11181]
ir_loc_arr = [5925, 4870, 8625, 5850]
#loc_arr = [-3945, -6085, -7860, 1310, 4275, 7160]
#blue in cm 
b_lam = 445e-7
b_del = 94e-7
#v in cm 
v_lam = 551e-7
v_del = 88e-7
#red in cm 
r_lam = 658e-7
r_del = 138e-7
#k in cm 
k_lam = 2190e-7
k_del=390e-7
#h in cm
h_lam = 1630e-7
h_del = 307e-7
#j in cm 
j_lam = 1220e-7
j_del = 213e-7

effarea = np.pi*(6.05/2)**2 #in cm^2
throughput = 0.6*0.42 #lens throughput (0.6) * quantum efficiency (0.42)
asteria_thru = 0.336
eff_thru_CUTE = 24 #cm^2 effarea*throughput for CUTE 

int_time = 7200 #s
lam_arr = [b_lam, v_lam, r_lam, k_lam, h_lam, j_lam]
del_arr = [b_del, v_del, r_del, k_del, h_del, j_del]
loc_arr = [-3945, -6085, -7860, 1310, 4275, 7160]

snr_base = [] ; snr_uvis = []; snr_ir = []
for i in range(len(lam_arr_uvis)):
    mask = (1e-4*wavelength[::-1] > (lam_arr_uvis[i] - del_arr_uvis[i]/2)) & (1e-4*wavelength[::-1] < (lam_arr_uvis[i] + del_arr_uvis[i]/2))
    Phi = (1.0 / (h * c)) * lam_arr_uvis[i]*(np.trapz(fp[::-1][mask],1e-4*wavelength[::-1][mask])) #wavelength in cm 
    signal = Phi*effarea*throughput*int_time*np.pi*(rp/sm)**2
    snr_signal = signal/np.sqrt(signal)
    precision = 100/snr_signal
    #snr.append(signal/np.sqrt(signal))
    snr_uvis.append(precision)

for i in range(len(lam_arr_ir)):
    mask = (1e-4*wavelength[::-1] > (lam_arr_ir[i] - del_arr_ir[i]/2)) & (1e-4*wavelength[::-1] < (lam_arr_ir[i] + del_arr_ir[i]/2))
    Phi = (1.0 / (h * c)) * lam_arr_ir[i]*(np.trapz(fp[::-1][mask],1e-4*wavelength[::-1][mask])) #wavelength in cm 
    signal = Phi*effarea*throughput*int_time*np.pi*(rp/sm)**2
    snr_signal = signal/np.sqrt(signal)
    precision = 100/snr_signal
    #snr.append(signal/np.sqrt(signal))
    snr_ir.append(precision)


for i in range(len(lam_arr)):
    signal = fp[loc_arr[i]]*(lam_arr[i]*del_arr[i]/(h*c))*eff_thru_CUTE*int_time*np.pi*(rp/sm)**2
    snr_signal = signal/np.sqrt(signal)
    precision = 100/snr_signal
    #snr.append(signal/np.sqrt(signal))
    snr_base.append(precision)
#str_bands_uvis = ['F625W','F475W','F410M','F467M','F547M', 'F814W', 'F850LP']
str_bands_uvis = ['F625W','F475W','F814W', 'F850LP']
str_bands_ir = ['F139M','F160W','F105W','F125W']
#print("B V R K H J")
#print(snr_base)
print(str_bands_uvis)
print(snr_uvis)
print(str_bands_ir)
print(snr_ir)
pdb.set_trace()
filter_name = ['B', 'V', 'R', 'K', 'H', 'J']
colors = ['blue', 'green', 'red', 'black', 'purple', 'orange']

#for i in range(len(loc_arr)):
    #plt.plot([np.array(lam_arr[i]-del_arr[i]/2)*1e4, np.array(lam_arr[i]-del_arr[i]/2)*1e4, np.array(lam_arr[i]+del_arr[i]/2)*1e4, np.array(lam_arr[i]+del_arr[i]/2)*1e4], [0, throughput, throughput, 0], colors[i], label = filter_name[i], alpha = 0.5, ls = '-.')
    #plt.fill_between([np.array(lam_arr[i]-del_arr[i]/2)*1e4, np.array(lam_arr[i]+del_arr[i]/2)*1e4], [throughput, throughput],  step="pre", color = colors[i], alpha = 0.3, label = filter_name[i])

for i in range(len(lam_arr_uvis)):
    plt.plot([np.array(lam_arr_uvis[i]-del_arr_uvis[i]/2)*1e4, np.array(lam_arr_uvis[i]-del_arr_uvis[i]/2)*1e4, np.array(lam_arr_uvis[i]+del_arr_uvis[i]/2)*1e4, np.array(lam_arr_uvis[i]+del_arr_uvis[i]/2)*1e4], [0, asteria_thru, asteria_thru, 0],label = str_bands_uvis[i], ls = '--', lw = 2)
    #plt.fill_between([np.array(lam_arr_uvis[i]-del_arr_uvis[i]/2)*1e7, np.array(lam_arr_uvis[i]+del_arr_uvis[i]/2)*1e7], [uvis_thru[i], uvis_thru[i]], step="pre", alpha = 0.5, label = str_bands_uvis[i])

#for i in range(len(lam_arr_ir)):
    #plt.plot([np.array(lam_arr_ir[i]-del_arr_ir[i]/2)*1e4, np.array(lam_arr_ir[i]-del_arr_ir[i]/2)*1e4, np.array(lam_arr_ir[i]+del_arr_ir[i]/2)*1e4, np.array(lam_arr_ir[i]+del_arr_ir[i]/2)*1e4], [0, ir_thru[i], ir_thru[i], 0],label = str_bands_ir[i], ls = '--', lw = 2)
    #plt.fill_between([np.array(lam_arr_uvis[i]-del_arr_uvis[i]/2)*1e7, np.array(lam_arr_uvis[i]+del_arr_uvis[i]/2)*1e7], [uvis_thru[i], uvis_thru[i]], step="pre", alpha = 0.5, label = str_bands_uvis[i])

plt.legend()
plt.xlabel("Wavelength (micron)")
plt.plot(1e4/wno_base, fpfs_base*1e10, 'k')

plt.ylabel("System Throughput")
plt.title("System Throughputs per Filter")
plt.ylim(0,0.8)

plt.show()
plt.close()
pdb.set_trace()

#plt.plot(1e4/wno_base, fpfs_base, 'k')
