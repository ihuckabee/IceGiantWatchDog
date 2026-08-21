import numpy as np
#don't know if these are necessary for this code but picaso yells at me 
#when i try to import picaso and don't do this, regardless of if i use a ref input
#so i define the reference data location to not get yelled at
import os
os.environ["picaso_refdata"]="/Users/ihuckabee/Documents/picaso-master/reference"
os.environ['PYSYN_CDBS']="/Users/ihuckabee/Documents/picaso-master/grp/redcat/trds"
refdata = os.getenv("picaso_refdata")
import pdb
import sys
sys.path.append("..")
from virga import justdoit as vdi
import pandas as pd
import astropy.units as u
from picaso import justplotit as ppi
from picaso import justdoit as pdi
import matplotlib.pyplot as plt
import copy



#SETTING UP THE PLANET
opacity = pdi.opannection(wave_range=[0.3,2.5])
neptune = pdi.inputs()
neptune.phase_angle(0)
neptune.gravity(radius=0.3463884070945,radius_unit=pdi.u.Unit('R_jup'), mass=0.053740779768177, mass_unit=pdi.u.Unit('M_jup'))
#neptune.star(opacity, temp=5778, metal=0.00, logg=4.4374, semi_major = 50., semi_major_unit = u.Unit("au"), radius=1.0, radius_unit=pdi.u.R_sun,)
neptune.star(opacity, temp=5778, metal=0.00, logg=4.4374, semi_major = 30.178, semi_major_unit = u.Unit("au"), radius=1.0, radius_unit=pdi.u.R_sun,)
neptune.atmosphere(filename="/Users/ihuckabee/Documents/picaso-master/neptune_stuff/ToShare/RefFiles/nepforpicaso_salr_abridged.txt", delim_whitespace=True)

#neptune.chemeq_visscher(c_o = 1., log_mh = 0.) #solar but ummmmmmm... a quick google did not find neptune's actual c_o or met. probably an issue
#df = neptune.spectrum(opacity)

#cloud free spectrum
#wno, alb, fpfs = df['wavenumber'] , df['albedo'] , df['fpfs_reflected']
#wno, alb = pdi.mean_regrid(wno, alb , R=100)

#reading in cloud info calculated previously 
df_1 = pd.read_csv('/Users/ihuckabee/Documents/picaso-master/neptune_stuff/ToShare/RefFiles/cloud1_v3.csv')
df_2 = pd.read_csv('/Users/ihuckabee/Documents/picaso-master/neptune_stuff/ToShare/RefFiles/cloud2_v3.csv')
df_3 = pd.read_csv('/Users/ihuckabee/Documents/picaso-master/neptune_stuff/ToShare/RefFiles/cloud3_v3.csv')
df_4 = pd.read_csv('/Users/ihuckabee/Documents/picaso-master/neptune_stuff/ToShare/RefFiles/cloud4_v3.csv')
#df_32 = pd.read_csv('/Users/ihuckabee/Documents/picaso-master/ToShare/RefFiles/cloud32_v2.csv')
#pdb.set_trace()
#cloud base
df_base = copy.deepcopy(df_1[["pressure", "wavenumber"]])
df_base['g0'] = df_1['g0']+df_2['g0']+df_3['g0']#+df_4['g0']
df_base['w0'] = df_1['w0']+df_2['w0']+df_3['w0']#+df_4['w0']
df_base['opd'] = df_1['opd']+df_2['opd']+df_3['opd']#+df_4['opd']

#df_base['g0'] = df_2['g0']+df_3['g0']#+df_4['g0']
#df_base['w0'] = df_2['w0']+df_3['w0']#+df_4['w0']
#df_base['opd'] = df_2['opd']+df_3['opd']#+df_4['opd']

neptune.clouds(df=df_base.astype(float))
df_spec = neptune.spectrum(opacity, full_output = True)
wno_base, alb_base, fpfs_base = df_spec['wavenumber'] , df_spec['albedo'] , df_spec['fpfs_reflected']



#fs_base =df_spec['full_output']['star']['flux']
#wno_base, alb_base = pdi.mean_regrid(wno_base, alb_base, R=100)

#wno_base, fs_base = pdi.mean_regrid(wno_base_og, df_spec['full_output']['star']['flux'], R=100)
#can also get the spectrum of clouds individually if you want like this:
#neptune.clouds(df=df_base.astype(float))
#df_spec = neptune.spectrum(opacity)
#pdb.set_trace()
#ppi.show(ppi.plot_cld_input(24, 66,df=neptune.inputs['clouds']['profile']))
'''

#wno_4, alb_4, fpfs_4 = df_spec['wavenumber'] , df_spec['albedo'] , df_spec['fpfs_reflected']
#wno_4, alb_4 = pdi.mean_regrid(wno_4, fpfs_4, R=100)

#pdb.set_trace()
#save spectrum to read into plottingoutputs.py
#np.save("./RefFiles/1D_darkcloud.npy", [wno_base, alb_base])
#np.save("./RefFiles/1Dcloud4_v3.npy", [wno_base, alb_base])
#np.save("./RefFiles/MIRalbedospec.npy", [wno_base, fpfs_base, fs_base])
#np.save("./RefFiles/spec_thrutest.npy", [wno_base, fpfs_base, df_spec['full_output']['star']['flux']])

pdb.set_trace()

#np.save("./RefFiles/1Dtaucldtest.npy",np.array([wno_base,df_spec['full_output']['taucld'][52].T[0]]))
#pdb.set_trace() #enter c to continue or can look at variable values or whateva
sm = 4.514565e12 #m
rp = 2.47639e7 #m
fs = df_spec['full_output']['star']['flux']
fp = fpfs_base*fs #erg/cm2/s/cm
#convert to photon flux - *delta_lam * lambda / hc
h = 6.6261e-27 #cm2 g s-1 
c = 3e10 #cm s-1 

#HST UVIS
F625W = np.array([0.4777, 0.7707])*1e-7  # (6242.6 ± 1464.6/2) Å → microns -> cm
F475W = np.array([0.4101, 0.5445])*1e-7   # (4773.1 ± 1343.5/2) Å → microns  
F410M = np.array([0.3938, 0.4281])*1e-7   # (4109 ± 172/2) Å → microns
F467M = np.array([0.4583, 0.4783])*1e-7   # (4682.6 ± 200/2) Å → microns
F547M = np.array([0.5123, 0.5773])*1e-7  # (5447.5 ± 650/2) Å → microns

#HST IR
F139M = np.array([1.3511, 1.4154])*1e-7   # (1383.3 ± 64.3/2) nm → microns
F160W = np.array([1.4027, 1.6710])*1e-7   # (1536.9 ± 268.3/2) nm → microns
F105W = np.array([0.9227, 1.1877])*1e-7   # (1055.2 ± 265/2) nm → microns
F125W = np.array([1.2500, 1.5345])*1e-7   # (1392.3 ± 284.5/2) nm → microns

uvis_bands = [F625W,F475W,F410M,F467M,F547M]
ir_bands = [F139M,F160W,F105W,F125W]
lam_arr = []; del_arr = []
#lam_arr
#[6.242e-08, 4.773e-08, 4.1095e-08, 4.6829999999999996e-08, 5.4479999999999996e-08, 1.38325e-07, 1.5368500000000002e-07, 1.0552e-07, 1.3922499999999998e-07]
for i in range(len(uvis_bands)):
    lam_arr.append((uvis_bands[i][0]+uvis_bands[i][1])/2)
    del_arr.append(uvis_bands[i][1]-uvis_bands[i][0])
#for i in range(len(ir_bands)):
    #lam_arr.append((ir_bands[i][0]+ir_bands[i][1])/2)
    #del_arr.append(ir_bands[i][1]-ir_bands[i][0])
#uvis_thru = [0.28, 0.27, 0.22,0.28,0.27]
ir_thru = [0.54, 0.56, 0.52, 0.56]
uvis_loc_arr = [-7327, -4645, -3148, -4450, -5966]
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
int_time = 7200 #s
#lam_arr = [b_lam, v_lam, r_lam, k_lam, h_lam, j_lam]
#del_arr = [b_del, v_del, r_del, k_del, h_del, j_del]
loc_arr = [-3945, -6085, -7860, 1310, 4275, 7160]

snr = []
for i in range(len(lam_arr)):
    signal = fp[uvis_loc_arr[i]]*(lam_arr[i]*del_arr[i]/(h*c))*effarea*throughput*int_time*np.pi*(rp/sm)**2
    #signal = fp[loc_arr[i]]*(lam_arr[i]*del_arr[i]/(h*c))*24*int_time*np.pi*(rp/sm)**2
    snr_signal = signal/np.sqrt(signal)
    precision = 100/snr_signal
    #snr.append(signal/np.sqrt(signal))
    snr.append(precision)
str_bands_uvis = ['F625W','F475W','F410M','F467M','F547M']
#str_bands_ir = ['F139M','F160W','F105W','F125W']
print("B V R K H J")
#print(str_bands_uvis)
print(snr)
pdb.set_trace() 



'''

f = open("./RefFiles/obsrefspec.txt",'r')
lines = f.readlines()[110:] #data doesn't start til 0.3 microns
wv_obs = []; alb_obs = []
for i in range(len(lines)):
    wv_obs.append(float(lines[i].split()[0]))
    alb_obs.append(float(lines[i].split()[1]))

f = open("background.txt",'r')
lines = f.readlines()[1:] 
wv_back = []; alb_back = []
for i in range(len(lines)):
    wv_back.append(float(lines[i].split(',')[0]))
    alb_back.append(float(lines[i].split(',')[1]))

f = open("nds-2018.txt",'r')
lines = f.readlines()[1:] 
wv_nds = []; alb_nds = []
for i in range(len(lines)):
    wv_nds.append(float(lines[i].split(',')[0]))
    alb_nds.append(float(lines[i].split(',')[1]))



#wno_base, alb_base = pdi.mean_regrid(wno_base, alb_base, R=100)

plt.plot(1e4/wno_base,alb_base,label="1D Cloud Base")
plt.plot(wv_obs,alb_obs, 'k--', label="HST/SpeX Data", alpha = 0.75)
plt.plot(np.array(wv_back)/1000,alb_back,'r-', label="Background Fit ", alpha = 0.8)
plt.plot(wv_nds,alb_nds,'b-', label="NDS Fit", alpha = 0.8)
plt.legend()
plt.show()


