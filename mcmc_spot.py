import numpy as np
import emcee
import matplotlib.pyplot as plt
import os
#os.environ["picaso_refdata"] = "/home/izzyh/Documents/picaso-master/reference"
#os.environ["PYSYN_CDBS"] = "/home/izzyh/Documents/picaso-master/grp/redcat/trds"
os.environ["picaso_refdata"]="/Users/ihuckabee/Documents/picaso-master/reference"
os.environ['PYSYN_CDBS']="/Users/ihuckabee/Documents/picaso-master/grp/redcat/trds"
from scipy.special import expit  # logistic function to smoothly convert real values to [0,1]
refdata = os.getenv("picaso_refdata")
import sys
sys.path.append("..")
import pdb
from virga import justdoit as vdi
import pandas as pd
import astropy.units as u
from picaso import justplotit as ppi
from picaso import justdoit as pdi
import pandas as pd
import xarray as xr
import time
from bokeh.plotting import show, figure
import pickle
import math
import copy
import matplotlib.pyplot as plt
from disco_mcmc import lc_forMCMC


def log_prior(theta):
    num_spots = int(theta[0])
    if not (0 <= num_spots <= 5): return -np.inf

    spot_params = theta[1:]
    for i in range(num_spots):
        lat, lon, dlat_frac, dlon_frac, spot_type_logit = spot_params[i*5:(i+1)*5]
        if not (-90 <= lat <= 90): return -np.inf
        if not (0 <= lon <= 360): return -np.inf
        if not (0.01 <= dlat_frac <= 0.5): return -np.inf
        if not (0.01 <= dlon_frac <= 0.5): return -np.inf
        # spot_type is modeled as a logit transformed real number
    return 0.0

def log_likelihood(theta, time_array, observed_lc, lc_err):
    num_spots = int(theta[0])
    spot_params = theta[1:]
    
    spots = []
    for i in range(num_spots):
        lat, lon, dlat_frac, dlon_frac, spot_type_logit = spot_params[i*5:(i+1)*5]
        dlat = dlat_frac * 180
        dlon = dlon_frac * 360
        spot_type = int(expit(spot_type_logit) > 0.5)
        spots.append([lat, lon, dlat, dlon, spot_type])
    
    spot_map = #generate_spot_map(spots)  num spots, all that
    model_lc = lc_forMCMC(spots_dict) #simulate_light_curve(spot_map, time_array)
    
    return -0.5 * np.sum(((observed_lc - model_lc) / lc_err)**2)

def log_posterior(theta, time_array, observed_lc, lc_err):
    lp = log_prior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + log_likelihood(theta, time_array, observed_lc, lc_err)

def run_mcmc(time_array, observed_lc, lc_err, n_walkers=50, n_steps=1000):
    ndim = 1 + 5*5  # max: 1 (number of spots) + 5 params * 5 spots
    p0 = []
    
    for i in range(n_walkers):
        num_spots = np.random.randint(0, 6)
        spot_params = []
        for j in range(5):
            lat = np.random.uniform(-90, 90)
            lon = np.random.uniform(0, 360)
            dlat_frac = np.random.uniform(0.01, 0.5)
            dlon_frac = np.random.uniform(0.01, 0.5)
            spot_type_logit = np.random.normal()  # logit ~ real number
            spot_params += [lat, lon, dlat_frac, dlon_frac, spot_type_logit]
        theta = [num_spots] + spot_params
        p0.append(theta)
    
    sampler = emcee.EnsembleSampler(n_walkers, ndim, log_posterior, args=(time_array, observed_lc, lc_err))
    sampler.run_mcmc(p0, n_steps, progress=True)
    return sampler

def extract_results(sampler):
    flat_samples = sampler.get_chain(discard=200, thin=10, flat=True)
    best_index = np.argmax(sampler.get_log_prob(discard=200, thin=10, flat=True))
    best_params = flat_samples[best_index]
    
    num_spots = int(best_params[0])
    spots = []
    for i in range(num_spots):
        lat, lon, dlat_frac, dlon_frac, spot_type_logit = best_params[1+i*5:1+(i+1)*5]
        dlat = dlat_frac * 180
        dlon = dlon_frac * 360
        spot_type = int(expit(spot_type_logit) > 0.5)
        spots.append({
            "latitude": lat,
            "longitude": lon,
            "length_deg": dlat,
            "width_deg": dlon,
            "type": "bright" if spot_type else "dark"
        })
    return num_spots, spots


wave_range=[0.3,2.5]
opacity = pdi.opannection(wave_range=[0.3,2.5])
neptune = pdi.inputs()
nep_ptprofile = "/Users/ihuckabee/Documents/picaso-master/ToShare/RefFiles/nepforpicaso_salr_abridged.txt"
neptune = pdi.inputs()
grid_size = 2
neptune.phase_angle(0, num_gangle=grid_size, num_tangle=grid_size) #radians
neptune.gravity(radius=0.3463884070945,radius_unit=pdi.u.Unit('R_jup'), mass=0.053740779768177, mass_unit=pdi.u.Unit('M_jup'))
neptune.star(opacity, temp=5778, metal=0.00, logg=4.4374, semi_major = 30.178, semi_major_unit = u.Unit("au"), radius=1.0, radius_unit=pdi.u.R_sun,)

f = open(nep_ptprofile,'r')
lines = f.readlines()[1:]
pressure = []
temp = []
ch4 = []
h2s = []
h2 = []
he = []
for i in range(len(lines)):
        pressure.append(float(lines[i].split()[0]))
        temp.append(float(lines[i].split()[1]))
        ch4.append(float(lines[i].split()[2]))
        h2s.append(float(lines[i].split()[3]))
        h2.append(float(lines[i].split()[4]))
        he.append(float(lines[i].split()[5]))

#creating the lat and long grids for each species
pres = np.array(pressure) 
lon = np.linspace(-180,180,512)
lat = np.linspace(-90,90,128)
nep_data = np.tile(np.array(temp), (len(lon), len(lat), 1))
chem_ch4 = np.tile(np.array(ch4), (len(lon), len(lat), 1))
chem_h2s = np.tile(np.array(h2s), (len(lon), len(lat), 1))
chem_h2 = np.tile(np.array(h2), (len(lon), len(lat), 1))
chem_he = np.tile(np.array(he), (len(lon), len(lat), 1))
nwno = 23
wno_grid =np.array([ 4000.        ,  4166.66666667,  4347.82608696,  4545.45454545,
4761.9047619 ,  5000.        ,  5263.15789474,  5555.55555556,
5882.35294118,  6250.        ,  6666.66666667,  7142.85714286,
7692.30769231,  8333.33333333,  9090.90909091, 10000.        ,
11111.11111111, 12500.        , 14285.71428571, 16666.66666667,
20000.        , 25000.        , 33333.33333333])#np.linspace(1e4/wave_range[1],1e4/wave_range[0],nwno)

#reading in cloud info calculated previously 
df_1 = pd.read_csv('/Users/ihuckabee/Documents/picaso-master/ToShare/RefFiles/cloud1_v3.csv')
df_2 = pd.read_csv('/Users/ihuckabee/Documents/picaso-master/ToShare/RefFiles/cloud2_v3.csv')
df_3 = pd.read_csv('/Users/ihuckabee/Documents/picaso-master/ToShare/RefFiles/cloud3_v3.csv')
df_4 = pd.read_csv('/Users/ihuckabee/Documents/picaso-master/ToShare/RefFiles/cloud4_v3.csv')

df_base = copy.deepcopy(df_1[["pressure", "wavenumber"]])
df_base['g0'] = df_1['g0']+df_2['g0']+df_3['g0']#+df_32['g0']
df_base['w0'] = df_1['w0']+df_2['w0']+df_3['w0']#+df_32['w0']
df_base['opd'] = df_1['opd']+df_2['opd']+df_3['opd']#+df_32['opd']
opd_array = np.zeros((len(lon), len(lat), len(pres)-1, nwno)); g0_array = np.zeros((len(lon), len(lat), len(pres)-1, nwno)); w0_array = np.zeros((len(lon), len(lat), len(pres)-1, nwno))
test = np.zeros((len(lon), len(lat), len(pres)-1, nwno))
#iterate over each latitude and longitude point
opd_base = np.array(df_base['opd']).reshape(len(pres)-1,nwno)
g0_base= np.array(df_base['g0']).reshape(len(pres)-1,nwno)
w0_base = np.array(df_base['w0']).reshape(len(pres)-1,nwno)

opd_4 = np.array(df_4['opd']).reshape(len(pres)-1,nwno)
g0_4 = np.array(df_4['g0']).reshape(len(pres)-1,nwno)
w0_4 = np.array(df_4['w0']).reshape(len(pres)-1,nwno)

p_rot = 16.11  # hours
timestep = 1  # hours
num_rot = 10 
timedur = num_rot * p_rot  # hours
rnep = 0.3463884070945 * 7.1492e07  # meters
lon_len = len(lon)
lat_len = len(lat)

with open('windprofile.npy', 'rb') as f:
    xspeed, ylat = np.load(f)
lat_base = [int(lat_len/8), int(lat_len - lat_len/8)]
lon_base = [0, lon_len-1]

#defining spot dict
spots_dict = {}

#picaso object params
spots_dict['picaso_object']['planet'] = neptune
spots_dict['picaso_object']['opacity'] = opacity 

#rotation parameters
spots_dict['rot_params']['p_rot'] = p_rot 
spots_dict['rot_params']['timestep'] = timestep 
spots_dict['rot_params']['num_rot'] = num_rot 
spots_dict['rot_params']['rnep'] = rnep 

#spatial grid
spots_dict['spatial_grid']['lon_array'] = lon
spots_dict['spatial_grid']['lat_array'] = lat


#vertical profile
spots_dict['vertprofile']['pressure'] = pres
spots_dict['vertprofile']['temp_tile'] = nep_data
spots_dict['vertprofile']['CH4'] = chem_ch4 
spots_dict['vertprofile']['H2S'] = chem_h2s
spots_dict['vertprofile']['H2'] = chem_h2
spots_dict['vertprofile']['He'] = chem_he
spots_dict['vertprofile']['nwno'] = nwno

#wind profile
spots_dict['windprofile']['xspeed'] = xspeed
spots_dict['windprofile']['ylat'] = ylat

#cloud_info 
spots_dict['cloud_info']['cloud_base'] = df_base
spots_dict['cloud_info']['cloud_4'] = df_4
spots_dict['cloud_info']['wnogrid'] = wno_grid

time_array = np.linspace(0, 1, 100)
true_spot_info = #this is the spot info for what we're trying to fit 
true_lc = #this is what we're trying to fit 

lc_err = 0.01 * np.ones_like(true_lc)
noisy_lc = true_lc + np.random.normal(0, lc_err[0], size=true_lc.shape)

sampler = run_mcmc(time_array, noisy_lc, lc_err)
num_spots, spot_info = extract_results(sampler)

print(f"Estimated number of spots: {num_spots}")
for i, s in enumerate(spot_info):
    print(f"Spot {i+1}: {s}")

pdb.set_trace()