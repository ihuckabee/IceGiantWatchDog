import numpy as np
import os

#replace with path to ref files
os.environ["picaso_refdata"]="/Users/ihuckabee/Documents/picaso-master/reference"
os.environ['PYSYN_CDBS']="/Users/ihuckabee/Documents/picaso-master/grp/redcat/trds"

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
import copy
timestart = time.perf_counter()
wave_range=[0.3,2.5]
opacity = pdi.opannection(wave_range=[0.3,2.5])
neptune = pdi.inputs()
nep_ptprofile = "./RefFiles/nepforpicaso_salr_abridged.txt"
nep_clouds = "./RefFiles/cloudall_v2.csv"
hstspexdata = "./RefFiles/obsrefspec.txt"
cloudfree3d = "./RefFiles/3x3_cldfree.pic"
cloudfree1d = "./RefFiles/cloudfreespec.npy"
cloud1d = "./RefFiles/1Dcloudbase.npy"

neptune = pdi.inputs()

#define spatial grid
grid_size = 2
neptune.phase_angle(0, num_gangle=grid_size, num_tangle=grid_size) #radians
neptune.gravity(radius=0.3463884070945,radius_unit=pdi.u.Unit('R_jup'), mass=0.053740779768177, mass_unit=pdi.u.Unit('M_jup'))
neptune.star(opacity, temp=5778, metal=0.00, logg=4.4374, semi_major = 30.178, semi_major_unit = u.Unit("au"), radius=1.0, radius_unit=pdi.u.R_sun,)

#PT profile
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
lon = np.linspace(-180,180,128)
lat = np.linspace(-90,90,64)
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
#df_31 = pd.read_csv('/Users/ihuckabee/Documents/picaso-master/ToShare/RefFiles/cloud31.csv')
#df_32 = pd.read_csv('/Users/ihuckabee/Documents/picaso-master/ToShare/RefFiles/cloud32.csv')
pdb.set_trace()

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

#opd_4 = np.array(df_4['opd']).reshape(len(pres)-1,nwno)
#g0_4 = np.array(df_4['g0']).reshape(len(pres)-1,nwno)
#w0_4 = np.array(df_4['w0']).reshape(len(pres)-1,nwno)

#lat range: -90 to 90, 64 steps
#lon range: -180 to 180, 128 steps 
#define locations of clouds
lat_base = [0,64]
lon_base = [0,128]
lon_spot = [32,96] 
lat_spot = [30,34]
#lon_dspot = [16,48]
#lat_dspot =[80,90] 
#populate cloud param arrays to feed into xarray


for press in range(len(pres)-1):
    for wave in range(nwno):
        opd_array[lon_base[0]:lon_base[1], lat_base[0]:lat_base[1],press,wave] = opd_base[press,wave]
        g0_array[lon_base[0]:lon_base[1], lat_base[0]:lat_base[1],press,wave] = g0_base[press,wave]
        w0_array[lon_base[0]:lon_base[1], lat_base[0]:lat_base[1],press,wave] = w0_base[press,wave]
        test[lon_base[0]:lon_base[1], lat_base[0]:lat_base[1],press,wave] += 1

        opd_array[lon_spot[0]:lon_spot[1], lat_spot[0]:lat_spot[1],press,wave] += 0.0 #+= opd_4[press,wave]
        g0_array[lon_spot[0]:lon_spot[1], lat_spot[0]:lat_spot[1],press,wave] += 0.0 #+= g0_4[press,wave]
        w0_array[lon_spot[0]:lon_spot[1], lat_spot[0]:lat_spot[1],press,wave] += 0.0#+= w0_4[press,wave]
        test[lon_spot[0]:lon_spot[1], lat_spot[0]:lat_spot[1],press,wave] = 0 



import matplotlib.pyplot as plt
plt.imshow(test[:,:,0,0].T,cmap="bone")
plt.colorbar(label="arbitrary color value")
plt.xlabel("Lon")
plt.ylabel("Lat")
#plt.titl("Grid Visualization")
#plt.show()
plt.close()
#import pdb; pdb.set_trace()
#setting up xarrays
ds = xr.Dataset(
    data_vars=dict(
        temperature=(["lon", "lat","pressure"], nep_data,{'units': 'Kelvin'})#, required
        #kzz = (["x", "y","z"], gcm_out['kzz'])#could add other data components if wanted
    ),
    coords=dict(
        lon=(["lon"], lon,{'units': 'degrees'}),#required
        lat=(["lat"], lat,{'units': 'degrees'}),#required
        pressure=(["pressure"], pres,{'units': 'bar'})#required*
    ),
    attrs=dict(description="coords with vectors"),
)

ds_chem = xr.Dataset(
    data_vars=dict(
        CH4=(["lon", "lat","pressure"], chem_ch4,{'units': 'v/v'}),
        H2S=(["lon", "lat","pressure"], chem_h2s,{'units': 'v/v'}),
        H2=(["lon", "lat","pressure"], chem_h2,{'units': 'v/v'}),
        He=(["lon", "lat","pressure"], chem_he,{'units': 'v/v'}),

    ),
    coords=dict(
        lon=(["lon"], lon,{'units': 'degrees'}),#required
        lat=(["lat"], lat,{'units': 'degrees'}),#required
        pressure=(["pressure"], pres,{'units': 'bar'})#required*
    ),
    attrs=dict(description="coords with vectors"),
)

ds_withchem= ds.update(ds_chem)

#adding in clouds

ds_cld= xr.Dataset(
    data_vars=dict(
        opd=(["lon", "lat","pressure","wno"], opd_array,{'units': 'depth per layer'}),
        g0=(["lon", "lat","pressure","wno"], g0_array,{'units': 'none'}),
        w0=(["lon", "lat","pressure","wno"], w0_array,{'units': 'none'}),
    ),
    coords=dict(
        lon=(["lon"], lon,{'units': 'degrees'}),#required
        lat=(["lat"], lat,{'units': 'degrees'}),#required
        pressure=(["pressure"], pres[:-1],{'units': 'bar'}),#required
        wno=(["wno"], wno_grid,{'units': 'cm^(-1)'})#required for clouds NEEDS CHEWCKING FOR UNITS AND WHATENOT. what do the 24 wv points map to
    ),
    attrs=dict(description="coords with vectors"),
)

neptune.atmosphere_3d(ds_withchem,regrid=True,plot=False,verbose=True)
neptune.clouds_3d(ds_cld)
#pdb.set_trace()
#ds_cld['opd'].isel(pressure=50,wno=0).plot(x='lon',y='lat'); plt.show() ->> TO SHOW CLOUD AT GIVEN PRcES/WNO
#next line is what takes a long time
out3d = neptune.spectrum(opacity,calculation='reflected',dimension='3d',full_output=True)
wno,fpfs = out3d['full_output']['wavenumber'], out3d['fpfs_reflected'] #wno in cm^-1

pdb.set_trace()
outfile = open("./RefFiles/2x2_cldv3.pic","wb")
pickle.dump(np.array([wno,fpfs]),outfile)
#ppi.show(ppi.disco(out3d['full_output'] ,wavelength=[1.63], calculation='reflected'))
#wno, alb = pdi.mean_regrid(wno,alb,R=100)
plt.close()
#taucld = out3d['full_output']['taucld'][52].T[0][0][0]
#plt.plot(1e4/wno[::-1][:-1],np.log10(taucld[::-1][:-1]))
plt.plot(1e4/wno, fpfs)
plt.show()
#pdb.set_trace()

