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
from astropy.io import fits
comp_file = pdi.pd.read_csv("/Users/ihuckabee/Documents/picaso-master/neptune_stuff/ToShare/RefFiles/nepforpicaso_salr_abridged.txt", delim_whitespace=True)
nlayer = len(comp_file['pressure'])-1 #one less than the number of PT points in your input
df_1 = pd.read_csv('/Users/ihuckabee/Documents/picaso-master/neptune_stuff/ToShare/RefFiles/cloud1_trunc.csv')
df_2 = pd.read_csv('/Users/ihuckabee/Documents/picaso-master/neptune_stuff/ToShare/RefFiles/cloud2_trunc.csv')
df_3 = pd.read_csv('/Users/ihuckabee/Documents/picaso-master/neptune_stuff/ToShare/RefFiles/cloud3_trunc.csv')

df_base = copy.deepcopy(df_1[["pressure", "wavenumber"]])
df_base['g0'] = df_1['g0']+df_2['g0']+df_3['g0']#+df_4['g0']
df_base['w0'] = df_1['w0']+df_2['w0']+df_3['w0']#+df_4['w0']
df_base['opd'] = df_1['opd']+df_2['opd']+df_3['opd']#+df_4['opd']
df_no1 = copy.deepcopy(df_1[["pressure", "wavenumber"]])
df_no1['g0'] = df_1['g0']
df_no1['w0'] = df_1['w0']
df_no1['opd'] = df_1['opd']
df_darkspot_old = copy.deepcopy(df_1[["pressure", "wavenumber"]])
df_darkspot_old['g0'] = df_2['g0']+df_3['g0']#+df_4['g0']
df_darkspot_old['w0'] = df_2['w0']+df_3['w0']#+df_4['w0']
df_darkspot_old['opd'] = df_2['opd']+df_3['opd']#+df_4['opd']



#SETTING UP THE PLANET
opacity = pdi.opannection(wave_range=[0.3,1.0])
neptune = pdi.inputs()
neptune.phase_angle(0)
neptune.gravity(radius=0.3463884070945,radius_unit=pdi.u.Unit('R_jup'), mass=0.053740779768177, mass_unit=pdi.u.Unit('M_jup'))
#neptune.star(opacity, temp=5778, metal=0.00, logg=4.4374, semi_major = 50., semi_major_unit = u.Unit("au"), radius=1.0, radius_unit=pdi.u.R_sun,)
neptune.star(opacity, temp=5778, metal=0.00, logg=4.4374, semi_major = 30.178, semi_major_unit = u.Unit("au"), radius=1.0, radius_unit=pdi.u.R_sun,)
neptune.atmosphere(filename="/Users/ihuckabee/Documents/picaso-master/neptune_stuff/ToShare/RefFiles/nepforpicaso_salr_abridged.txt", delim_whitespace=True)



base_pressure_1 = 10. 
haze_top_pressure_1 = 2.05 
pressure = comp_file['pressure']
tau_ref = 0.65
wavelength = np.array([0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
wavenumber_grid = 1e4/wavelength[::-1]
Q_ext = np.array([7.49127, 3.50126, 1.88824, 1.00575, 0.791879, 1.0000, 0.299883, 0.19146])
w0 = np.array([0.889673, 0.955005, 0.93751, 0.990416, 0.719135, 0.362478, 0.958898, 0.948567])
g0 = np.array([0.62539445, 0.51918708, 0.38986472, 0.288171, 0.21448347, 0.16618509, 0.13258767, 0.10693117])

Q_ref = Q_ext[-3]

tau_lambda = tau_ref*(Q_ext/Q_ref)

in_cloud = (pressure <= base_pressure_1) & (pressure >= haze_top_pressure_1)
N = np.sum(in_cloud)

opd = np.zeros(len(wavelength))


opd= tau_lambda / N


df_dark_clouds= vdi.picaso_format_slab(base_pressure_1,opd, w0, g0, wavenumber_grid, pressure[:-1],p_top=haze_top_pressure_1)

df_darkspot_new = copy.deepcopy(df_1[["pressure", "wavenumber"]])
df_darkspot_new['g0'] = df_dark_clouds['g0'] + df_2['g0']+df_3['g0']#+df_4['g0']
df_darkspot_new['w0'] = df_dark_clouds['w0'] + df_2['w0']+df_3['w0']#+df_4['w0']
df_darkspot_new['opd'] = df_dark_clouds['opd'] +  df_2['opd']+df_3['opd']#+df_4['opd']

neptune_dark_new = copy.deepcopy(neptune)
neptune_dark_new.clouds(df=df_darkspot_new.astype(float))
df_cf = neptune.spectrum(opacity)
df_darkspot_new = neptune_dark_new.spectrum(opacity)

wno_cf, alb_cf, fpfs_cf = df_cf['wavenumber'] , df_cf['albedo'] , df_cf['fpfs_reflected']
wno_cf, alb_cf = pdi.mean_regrid(wno_cf, alb_cf, R=100)
wno_dark, alb_dark, fpfs_dark = df_darkspot_new['wavenumber'] , df_darkspot_new['albedo'] , df_darkspot_new['fpfs_reflected']
wno_dark, alb_dark = pdi.mean_regrid(wno_dark, alb_dark, R=100)

#ppi.show(ppi.spectrum([wno_dark,wno_cf], [alb_dark, alb_cf],legend = ["Darkened Aerosol 1", "Cloud-Free"], plot_width=750))


neptune_background = copy.deepcopy(neptune)
neptune_no1 = copy.deepcopy(neptune)
neptune_background.clouds(df=df_base.astype(float))
neptune_no1.clouds(df=df_no1.astype(float))
df_base= neptune_background.spectrum(opacity)
df_no1 = neptune_no1.spectrum(opacity)

neptune_darkspot_old = copy.deepcopy(neptune)
neptune_darkspot_old.clouds(df=df_darkspot_old.astype(float))
df_darkspot_old = neptune_darkspot_old.spectrum(opacity)

wno_base, alb_base, fpfs_base = df_base['wavenumber'] , df_base['albedo'] , df_base['fpfs_reflected']
wno_base, alb_base = pdi.mean_regrid(wno_base, alb_base, R=100)
wno_no1, alb_no1, fpfs_no1 = df_no1['wavenumber'] , df_no1['albedo'] , df_no1['fpfs_reflected']
wno_no1, alb_no1 = pdi.mean_regrid(wno_no1, alb_no1, R=100)
wno_darkspot_old, alb_darkspot_old, fpfs_darkspot_old = df_darkspot_old['wavenumber'] , df_darkspot_old['albedo'] , df_darkspot_old['fpfs_reflected']
wno_darkspot_old, alb_darkspot_old = pdi.mean_regrid(wno_darkspot_old, alb_darkspot_old, R=100)


#loading in vlt data. need solar spectrum first 
hdul_sun = fits.open(name='/Users/ihuckabee/Downloads/solar_spec.fits', memmap=False, cache=False, lazy_load_hdus=False)
mast_data_sun = hdul_sun[1].data
mast_wavelength_sun = mast_data_sun.WAVELENGTH #angstrom
mast_flux_sun = mast_data_sun.FLUX #erg s-1 cm-2 A-1 
#convert to W cm-2 um-1 

wno_ref, fpfs_ref, fs = np.load("./ToShare/RefFiles/spec_thrutest.npy")#fs in #erg/cm2/s/cm
wno_ref, fs = pdi.mean_regrid(wno_ref, fs, R=100)
sm = 4.514565e12  # m semi major axis of planet
rp = 2.47639e7  # m radius of planet
d_to_nep_close = 29*1.496e+11 #m 
d_to_nep_far = 31*1.496e+11 #m 
d_avg = (d_to_nep_far+d_to_nep_close)/2 #distance from observer to planet
rs = 6.96e8 #m radius of star
F_sun_erg = fs[::-1][:-93]*(rs/(d_avg))**2 #scaled to what we observe at Neptune
h = 6.6261e-27  # erg s
c = 3e10  # cm s-1


from scipy.interpolate import interp1d

# ── File path ────────────────────────────────────────────────────────────────
SPEC_FILE = "/Users/ihuckabee/Downloads/compare_spectra_nds.txt"


# Input: wavenumber grid in cm^-1, flux in erg/s/cm^2/cm
wave_solar_um = 1e4 / wno_ref[::-1][:-93]     # wavenumber (cm^-1) -> wavelength (um)
wave_solar_cm = wave_solar_um * 1e-4
# Convert erg/s/cm^2/cm -> W/cm^2/um
# F_lam = F_wn * (1/lambda_cm^2) with unit adjustment
# lambda_cm = lambda_um * 1e-4, so 1/lambda_cm^2 = 1/(lambda_um^2 * 1e-8)
# 1 erg/s = 1e-7 W
F_sun_W = F_sun_erg * 1e-7  #in units of W/cm2/cm 
F_sun_W = F_sun_W * 1e-4 #W/cm2/um 
F_sun_W = F_sun_W#/ (np.pi *(rp/sm)**2) #W/cm2/sr-1/um
#/ (wave_solar_um)**2 #* 1e-4

# ── User-supplied values ─────────────────────────────────────────────────────
# Neptune's heliocentric distance at time of observation (AU)

# ── Parse the file ───────────────────────────────────────────────────────────
bg_rows  = []   # background spectrum rows
nds_rows = []   # NDS-2018 spectrum rows

section = None
with open(SPEC_FILE) as fh:
    for line in fh:
        line = line.strip()
        if not line:
            continue
        if line.startswith("Fit to background"):
            section = "bg"
            continue
        if line.startswith("Fits to NDS"):
            section = "nds"
            continue
        if line.startswith("Wavelength") or line.startswith("Radiances") or line[0].isalpha():
            continue                        # skip all other header lines
        cols = line.split()
        if section == "bg" and len(cols) >= 4:
            bg_rows.append([float(c) for c in cols[:4]])
        elif section == "nds" and len(cols) >= 4:
            nds_rows.append([float(c) for c in cols[:5]])

bg  = np.array(bg_rows)   # columns: wave_nm, radiance, error, fit_radiance
nds = np.array(nds_rows)  # columns: wave_nm, radiance, error, fit1, fit2

# ── Wavelength arrays (A -> um) ─────────────────────────────────────────────
wave_bg_um  = bg[:, 0]  * 1e-3
wave_nds_um = nds[:, 0] * 1e-3

# ── Radiance arrays (already W cm^-2 sr^-1 um^-1) ───────────────────────────
rad_bg       = bg[:, 1]   # measured background radiance
fit_bg       = bg[:, 3]   # fitted background radiance

rad_nds      = nds[:, 1]  # measured NDS-2018 radiance
fit_nds      = nds[:, 3]  # fitted NDS-2018 radiance (model 1)

# ── Interpolate F_sun onto each spectrum's wavelength grid ───────────────────
# wave_solar_nm must be in nm to match; convert to um for the interpolator

interp_solar = interp1d(wave_solar_um, F_sun_W,
                        kind='linear', bounds_error=False,
                        fill_value=np.nan)

Fsun_on_bg  = interp_solar(wave_bg_um)
Fsun_on_nds = interp_solar(wave_nds_um)

# ── Compute I/F (albedo) ─────────────────────────────────────────────────────
# I/F = pi * L / F_sun_at_Neptune
# (pi converts radiance per sr to equivalent flux for a Lambertian reflector)
IF_bg_measured  = np.pi * rad_bg  /Fsun_on_bg
IF_bg_fit       = np.pi * fit_bg  / Fsun_on_bg
IF_nds_measured = np.pi * rad_nds / Fsun_on_nds
IF_nds_fit      = np.pi * fit_nds / Fsun_on_nds
import pdb; pdb.set_trace()
# ── Ready to plot ─────────────────────────────────────────────────────────────
# x-axis:  wave_bg_um  or  wave_nds_um  (wavelength in microns)
# y-axis:  IF_bg_fit, IF_nds_fit  (or the _measured variants)
#
# Example:
# import matplotlib.pyplot as plt
# fig, ax = plt.subplots()
# ax.plot(wave_bg_um,  IF_bg_fit,  label="Background fit")
# ax.plot(wave_nds_um, IF_nds_fit, label="NDS-2018 fit")
# ax.set_xlabel("Wavelength (μm)")
# ax.set_ylabel("I/F")
# ax.legend()
# plt.show()
#import pdb; pdb.set_trace()

f = open("./ToShare/RefFiles/obsrefspec.txt",'r')
lines = f.readlines()[110:] #data doesn't start til 0.3 microns
wv_obs = []; alb_obs = []
for i in range(len(lines)):
    wv_obs.append(float(lines[i].split()[0]))
    alb_obs.append(float(lines[i].split()[1]))


import matplotlib.pyplot as plt
#plt.plot(1e4/wno_dark, alb_dark, label = 'Darkened Aerosol 1')
#plt.plot(1e4/wno_no1, alb_no1, label = 'Aerosol 1')
#plt.plot(1e4/wno_cf, alb_cf, label = 'Cloud-Free')
plt.plot(1e4/wno_base, alb_base, '-.',color = 'b', alpha = 0.5, label = 'Cloud Background - PICASO')
#plt.plot(1e4/wno_darkspot_old, alb_darkspot_old, '-.',color = 'r', alpha = 0.5, label = 'Dark Spot - Old')
#plt.plot(1e4/wno_dark, alb_dark, '-',color = 'k', alpha = 0.75, label = 'Dark Spot - Updated')
plt.plot(wave_bg_um,  IF_bg_fit,color = 'b', label="Background fit")
plt.plot(wave_nds_um, IF_nds_fit,color = 'r', label="NDS-2018 fit")
plt.plot(wv_obs,alb_obs, color = 'k', alpha = 0.5, label = 'Data')
plt.legend()
plt.xlabel("Wavelength (um)")
plt.ylabel("Albedo")
plt.show()
import pdb; pdb.set_trace()