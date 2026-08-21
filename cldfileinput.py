import numpy as np
import os
#replace with path to ref files
os.environ["picaso_refdata"]="/Users/ihuckabee/Documents/picaso-master/reference"
os.environ['PYSYN_CDBS']="/Users/ihuckabee/Documents/picaso-master/grp/redcat/trds"

refdata = os.getenv("picaso_refdata")
import pdb
import copy
from virga import justdoit as vdi
import pandas as pd
import astropy.units as u
from picaso import justplotit as ppi
from picaso import justdoit as pdi

mieff_dir = '/Users/ihuckabee/Documents/picaso-master/picaso/refrind'

numaero = '2'
if numaero == '31' or numaero == '32':
    qext, qscat, cos_qscat, nwave, radius, wave_in = vdi.get_mie("formattedrefrac_3",directory=mieff_dir)
else:
    qext, qscat, cos_qscat, nwave, radius, wave_in = vdi.get_mie("formattedrefrac_"+numaero,directory=mieff_dir)

tau1 = 0.65; tau2 = 1.5; tau3 = 0.04; tau4 = 0.03 #mean taus

meanr13 = 0.05/1e4 #mean radii from irwin paper given in microns
meanr2 = 0.68/1e4 #converted to cm 
meanr4 = 2.5/1e4

mean_indexr13 = 53
mean_indexr2 = 553
mean_indexr4 = 2500



f = open("/Users/ihuckabee/Documents/picaso-master/neptune_stuff/nepcustomradii.txt",'r') #LOCAL
#f = open("/Users/ihuckabee/Documents/picaso-master/nepcustomradii.txt",'r') #ASH
lines = f.readlines()[1:]
i = 0; rsize = []; r13 = []; r2 = []; r4 = []; 
for i in range(len(lines)):
    rsize.append(float(lines[i].split()[0]))
    r13.append(float(lines[i].split()[1]))
    r2.append(float(lines[i].split()[2]))
    r4.append(float(lines[i].split()[4]))
f.close()
r13 = np.array(r13) 
r2 = np.array(r2) 
r4 = np.array(r4) 
rsize = np.array(rsize) 


#r is the grid of radius bins, n is the distribution per particle per radius bin 
r_arr = np.array([rsize.T,r13.T,r2.T,r4.T]).T
df = pd.DataFrame(r_arr, columns = ['size','146K','84K', '54.7K'])
non_nans_13 = np.where(np.isfinite(df['146K']))
r13 = df['size'].iloc[non_nans_13]
n13 = df['146K'].iloc[non_nans_13]

non_nans_2 = np.where(np.isfinite(df['84K']))
r2 = df['size'].iloc[non_nans_2]
n2 = df['84K'].iloc[non_nans_2]

non_nans_4 = np.where(np.isfinite(df['54.7K']))
r4 = df['size'].iloc[non_nans_4]
n4 = df['54.7K'].iloc[non_nans_4]

r13 = np.array(r13)
r2 = np.array(r2)  
r4 = np.array(r4)     


comp_file = pdi.pd.read_csv("/Users/ihuckabee/Documents/picaso-master/neptune_stuff/ToShare/RefFiles/nepforpicaso_salr_abridged.txt", delim_whitespace=True)
# ---------------------------------------------------------------
# INPUTS -- replace with your actual data
# ---------------------------------------------------------------

# Pressure grid (bar), layer midpoints, increasing with depth (or however
# your PICASO grid is ordered -- just be consistent)
pressure_edges = comp_file['pressure']
pres = pressure_edges.values   # or .values
pressure_mid = np.sqrt(pres[:-1] * pres[1:]) 


# Atmosphere profile needed for scale height -- pull from your PICASO
# atmosphere object (same grid as pressure_mid)
temperature = comp_file['temperature']          # K, shape (nlayer,)
temp = temperature.values

base_pressure_1 = 10. 
haze_top_pressure_1 = 2.05 
base_pressure_2 = 2.05
haze_top_pressure_2 = 1.6
base_pressure_31 = 1.6
haze_top_pressure_31 = 0.2
base_pressure_4 = 0.2
haze_top_pressure_4 = 0.08
base_pressure_32 = 0.08
haze_top_pressure_32 = 0.01

# Normalized extinction cross section vs wavelength (dimensionless,
# sigma(lambda)/sigma(0.8um)); this becomes your OPD wavelength shape
wavelength = wavelength = np.array([0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4])       # um, shape (nwave,)
if numaero == '1':
    tau_reference = tau1
    meanr = meanr13
    mean_indexr = mean_indexr13
    rad = r13
    numr = n13
    dtau_dP_08 = 10.     # bar^-1, at 0.8um, interpolated onto pressure_mid
    base_pressure = 10. 
    haze_top_pressure = 2.05 
elif numaero == '2':
    tau_reference = tau2
    meanr = meanr2
    mean_indexr = mean_indexr2
    numr = n2
    rad = r2
    dtau_dP_08 = 4.
    base_pressure = 3.
    haze_top_pressure = 1.  

elif numaero == '31':
    tau_reference = tau3
    meanr = meanr13
    mean_indexr = mean_indexr13
    numr = n13
    rad = r13
    dtau_dP_08 = 0.05
    base_pressure = 1.6
    haze_top_pressure = 0.2

elif numaero == '32':
    tau_reference = tau3
    meanr = meanr13
    mean_indexr = mean_indexr13
    numr = n13
    rad = r13
    dtau_dP_08 = 0.05
    base_pressure = 0.08
    haze_top_pressure = 0.01

           
elif numaero == '4':
    tau_reference = tau4
    meanr = meanr4
    mean_indexr = mean_indexr4
    numr = n4
    rad = r4
    dtau_dP_08 = 0.5
    base_pressure = 0.2
    haze_top_pressure = 0.08     


                                
delta_P = pres[1:] - pres[:-1]    # bar, shape (nlayer,)
tau_layer_08 = np.zeros(len(pressure_mid))
p_range = (haze_top_pressure, base_pressure)# bar, bounds over which tau_reference applies
# Known reference optical depth and the pressure range it applies to



in_range = (pressure_mid >= p_range[0]) & (pressure_mid <= p_range[1])
tau_layer_08[in_range] = dtau_dP_08 * delta_P[in_range]  
tau_computed = tau_layer_08.sum()
calibration_factor = tau_reference / tau_computed
print(f"tau_computed (uncalibrated) over range: {tau_computed:.4f}")
print(f"calibration factor f = {calibration_factor:.4f}")

tau_layer_ref_calibrated = tau_layer_08 * calibration_factor

Q_ext_08 = qext[-7][mean_indexr]
sig_ext_08 = np.pi*meanr**2*Q_ext_08

N_column = tau_layer_ref_calibrated / sig_ext_08 
ind_ndz = np.where(N_column != 0)[0][0]
#import pdb; pdb.set_trace()
opd_1,w0_1,g0_1,wavenumber_grid_1=vdi.calc_optics_user_r_dist(wave_in, N_column[ind_ndz] ,rad, u.nm, numr/100, qext, qscat, cos_qscat, )



df = vdi.picaso_format_slab(base_pressure,opd_1[:-1], w0_1[:-1], g0_1[:-1], wavenumber_grid_1[:-1], pressure_edges[:-1],p_top=haze_top_pressure)
import pdb; pdb.set_trace()
df.to_csv('df_'+numaero+'_reviewer.csv')


