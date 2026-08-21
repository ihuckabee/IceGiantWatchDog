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
"""
Compute per-layer cloud optical depth (OPD) from a digitized
extinction-coefficient-vs-pressure plot, extend to a wavelength grid
using a normalized extinction cross section, and calibrate against a
known reference tau value.

Fill in the INPUTS section with your actual arrays/values.
"""

comp_file = pdi.pd.read_csv("/Users/ihuckabee/Documents/picaso-master/neptune_stuff/ToShare/RefFiles/nepforpicaso_salr_abridged.txt", delim_whitespace=True)
# ---------------------------------------------------------------
# INPUTS -- replace with your actual data
# ---------------------------------------------------------------

# Pressure grid (bar), layer midpoints, increasing with depth (or however
# your PICASO grid is ordered -- just be consistent)
pressure_edges = comp_file['pressure']
pres = pressure_edges.values   # or .values
pressure_mid = np.sqrt(pres[:-1] * pres[1:]) 

# beta_ext read off your plot (already divided by 100), km^-1, at 0.8 um
# one value per layer (or interpolate your digitized curve onto pressure_mid)
beta_ext = np.zeros(len(pressure_mid))              # shape (nlayer,)
beta_ext[-14:-12] = 11/100; beta_ext[-16:-14] = 10/100; 
# Atmosphere profile needed for scale height -- pull from your PICASO
# atmosphere object (same grid as pressure_mid)
temperature = comp_file['temperature']          # K, shape (nlayer,)
temp = temperature.values
temperature_mid = np.sqrt(temp[:-1] * temp[1:]) 
molar_mass = {
    'H2': 2.016, 'He': 4.003, 
    'CH4': 16.043, 'H2S': 34.0809
    
    # add whatever species you have
}

mixingratios=comp_file[['CH4', 'H2S', 'H2', 'He']].copy()
mmw = sum(mixingratios[gas] * molar_mass[gas] for gas in mixingratios)   # mean molecular weight, shape (nlayer,)
mmw_vals = mmw.values
mmw_mid =  np.sqrt(mmw_vals[:-1] * mmw_vals[1:]) 
g_planet = 11.1*1e3 #m/s2                           # surface/local gravity, cm/s^2 (or m/s^2 -- pick units and stay consistent)

# Normalized extinction cross section vs wavelength (dimensionless,
# sigma(lambda)/sigma(0.8um)); this becomes your OPD wavelength shape
wavelength = wavelength = np.array([0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4])       # um, shape (nwave,)
sigma_norm = np.array([0.0, 
    7.49127, 3.50126, 1.88824, 1.00575,
    0.791879, 1.00000, 0.299883, 0.19146, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
])           # shape (nwave,), sigma_norm(0.8um) == 1

# Known reference optical depth and the pressure range it applies to
tau_reference = 0.65
p_range = (2.05, 10)               # bar, bounds over which tau_reference applies

# ---------------------------------------------------------------
# CONSTANTS
# ---------------------------------------------------------------
k_B = 1.380649e-16   # erg/K  (cgs, pairs with g in cm/s^2, mmw in g/mol -> need /N_A)
N_A = 6.02214076e23  # 1/mol

# ---------------------------------------------------------------
# STEP 1: scale height + layer thickness (hydrostatic, cgs -> convert to km)
# ---------------------------------------------------------------
def scale_height_km(T, mu, g):
    """T in K, mu in g/mol, g in cm/s^2 -> H in km"""
    H_cm = (k_B * T) / ((mu / N_A) * g)
    return H_cm * 1e-5  # cm -> km

H = scale_height_km(temperature_mid, mmw_mid, g_planet)   # shape (nlayer,)

# dz per layer using ln(P_bot/P_top); make sure pressure_edges is ordered
# top-of-atmosphere -> deep, so P_bot > P_top for each layer
dz_km = H * np.log(pres[1:] / pres[:-1])  # shape (nlayer,)
dz_km = np.abs(dz_km)  # guard against ordering sign flips

# ---------------------------------------------------------------
# STEP 2: optical depth per layer at 0.8 um
# ---------------------------------------------------------------
tau_layer_ref = beta_ext * dz_km   # shape (nlayer,)

# ---------------------------------------------------------------
# STEP 3: calibrate against known tau_reference
# ---------------------------------------------------------------
in_range = (pressure_mid >= p_range[0]) & (pressure_mid <= p_range[1])
tau_computed = tau_layer_ref[in_range].sum()

calibration_factor = tau_reference / tau_computed
print(f"tau_computed (uncalibrated) over range: {tau_computed:.4f}")
print(f"calibration factor f = {calibration_factor:.4f}")

tau_layer_ref_calibrated = tau_layer_ref * calibration_factor

# ---------------------------------------------------------------
# STEP 4: extend to full wavelength grid -> OPD(nlayer, nwave)
# ---------------------------------------------------------------
opd = tau_layer_ref_calibrated[:, None] * sigma_norm[None, :-1]  # shape (nlayer, nwave)

# ---------------------------------------------------------------
# sanity checks
# ---------------------------------------------------------------
#assert opd.shape == (len(pressure_mid), len(wavelength))
print("OPD array shape:", opd.shape)
print("Total column tau at 0.8um (calibrated):", tau_layer_ref_calibrated.sum())



# opd, along with your existing w0 (nlayer, nwave) and g0 (nlayer, nwave),
# is now ready to feed into PICASO's cloud structure, e.g.:
#   df = pd.DataFrame({'lvl':..., 'opd':opd.flatten(), 'w0':w0.flatten(),
#                       'g0':g0.flatten(), 'wave':...})
#   case.clouds(df=df)   # or however your PICASO version expects it

w0 = np.array([0.0,
    0.889673, 0.955005, 0.93751, 0.990416,
    0.719135, 0.362478, 0.958898, 0.948567,  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
])  # single-scattering albedo

g1 = np.array([
    0.0, 0.625467, 0.519249, 0.396603, 0.345827,
    0.32571, 0.316067, 0.309884, 0.306829,  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
]) 

w0grid = tau_layer_ref_calibrated.copy(); w0grid[in_range] = 1; w0grid = w0grid[:, None] * w0[None,:]
g0grid = tau_layer_ref_calibrated.copy(); g0grid[in_range] = 1; g0grid = g0grid[:, None] * g1[None,:]



wavenumber = np.array([ 4166.66666667,  4347.82608696,  4545.45454545,  4761.9047619 ,
        5000.        ,  5263.15789474,  5555.55555556,  5882.35294118,
        6250.        ,  6666.66666667,  7142.85714286,  7692.30769231,
        8333.33333333,  9090.90909091, 10000.        , 11111.11111111,
       12500.        , 14285.71428571, 16666.66666667, 20000.        ,
       25000.        , 33333.33333333, 50000.        ])



df_1 = pd.read_csv('df_2_reviewer.csv')
df_1mod = copy.deepcopy(df_1[["pressure", "wavenumber"]])

df_1mod['opd'] = opd.flatten()
df_1mod['w0'] = w0grid.flatten()
df_1mod['g0'] = g0grid.flatten()

df_1mod.to_csv('df_1mod.csv')
