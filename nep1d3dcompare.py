import pandas as pd
import numpy as np
import os
os.environ["picaso_refdata"]="/Users/ihuckabee/Documents/picaso-master/reference"
os.environ['PYSYN_CDBS']="/Users/ihuckabee/Documents/picaso-master/grp/redcat/trds"
refdata = os.getenv("picaso_refdata")
from picaso import justdoit as pdi
from picaso import justplotit as ppi
from virga import justdoit as vdi
import pdb
import xarray as xr
import astropy.units as u
import matplotlib.pyplot as plt
wave_range=[0.3,1.]
opacity = pdi.opannection(wave_range)
opa = pdi.opannection(wave_range=wave_range)
neptune3d = pdi.inputs()
neptune1d = pdi.inputs()
#define spatial grid
grid_size = 4
neptune3d.phase_angle(0, num_gangle=grid_size, num_tangle=grid_size) #radians
neptune3d.gravity(radius=0.3463884070945,radius_unit=pdi.u.Unit('R_jup'), mass=0.053740779768177, mass_unit=pdi.u.Unit('M_jup'))
neptune3d.star(opacity, temp=5778, metal=0.00, logg=4.4374, semi_major = 30.178, semi_major_unit = u.Unit("au"), radius=1.0, radius_unit=pdi.u.R_sun,)
#litrly exact copy
neptune1d.phase_angle(0, num_gangle=grid_size, num_tangle=grid_size) #radians
neptune1d.gravity(radius=0.3463884070945,radius_unit=pdi.u.Unit('R_jup'), mass=0.053740779768177, mass_unit=pdi.u.Unit('M_jup'))
neptune1d.star(opacity, temp=5778, metal=0.00, logg=4.4374, semi_major = 30.178, semi_major_unit = u.Unit("au"), radius=1.0, radius_unit=pdi.u.R_sun,)
neptune1d.atmosphere(filename="/Users/ihuckabee/Documents/picaso-master/ToShare/RefFiles/nepforpicaso_salr_abridged.txt", delim_whitespace=True)
#setting up chemistry (now with 67 pressure layers instead of 999, just took the temp/chem info from every nth layer)
f = open("/Users/ihuckabee/Documents/picaso-master/ToShare/RefFiles/nepforpicaso_salr_abridged.txt",'r')
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

#works out so h2 and he are in that 80% and 19% ratio ish
chem_other = 1 - np.array(ch4)+np.array(h2s)
#h2 = 0.81*chem_other
#he = 0.19*chem_other

#creating the lat and long grids for each species
pres = np.array(pressure) 
lon = np.linspace(-180,180,128)
lat = np.linspace(-90,90,64)
nep_data = np.tile(np.array(temp), (len(lon), len(lat), 1))
chem_ch4 = np.tile(np.array(ch4), (len(lon), len(lat), 1))
chem_h2s = np.tile(np.array(h2s), (len(lon), len(lat), 1))
chem_h2 = np.tile(np.array(h2), (len(lon), len(lat), 1))
chem_he = np.tile(np.array(he), (len(lon), len(lat), 1))
wno_grid = 1e4/np.arange(0.2,2.6,0.1)
nwno = 24 

#reading in cloud info calculated previously 
df_1 = pd.read_csv('/Users/ihuckabee/Documents/picaso-master/ToShare/RefFiles/cloud1.csv')
df_2 = pd.read_csv('/Users/ihuckabee/Documents/picaso-master/ToShare/RefFiles/cloud2.csv')
df_31 = pd.read_csv('/Users/ihuckabee/Documents/picaso-master/ToShare/RefFiles/cloud31.csv')
df_4 = pd.read_csv('/Users/ihuckabee/Documents/picaso-master/ToShare/RefFiles/cloud4.csv')
df_32 = pd.read_csv('/Users/ihuckabee/Documents/picaso-master/ToShare/RefFiles/cloud32.csv')

df_base = df_1[["pressure", "wavenumber"]]
df_base['opd'] = df_1['opd']+df_2['opd']+df_31['opd']+df_32['opd']
df_base['g0'] = df_1['g0']+df_2['g0']+df_31['g0']+df_32['g0']
df_base['w0'] = df_1['w0']+df_2['w0']+df_31['w0']+df_32['w0']

opd_array = np.zeros((len(lon), len(lat), len(pres)-1, nwno)); g0_array = np.zeros((len(lon), len(lat), len(pres)-1, nwno)); w0_array = np.zeros((len(lon), len(lat), len(pres)-1, nwno))
#iterate over each latitude and longitude point
opd_base = np.array(df_base['opd']).reshape(len(pres)-1,nwno)
g0_base= np.array(df_base['g0']).reshape(len(pres)-1,nwno)
w0_base = np.array(df_base['w0']).reshape(len(pres)-1,nwno)

opd_4 = np.array(df_4['opd']).reshape(len(pres)-1,nwno)
g0_4 = np.array(df_4['g0']).reshape(len(pres)-1,nwno)
w0_4 = np.array(df_4['w0']).reshape(len(pres)-1,nwno)

#define locations of clouds
lat_base = [16,48]
lon_spot = [48,80]
lat_spot = [24,40]
#populate cloud param arrays to feed into xarray
for lo in range(len(lon)):
    for la in range(len(lat)):
        for press in range(len(pres)-1):
            for wave in range(nwno):
                opd_array[lo, la, press, wave] = opd_base[press][wave]
                g0_array[lo, la, press, wave] = g0_base[press][wave]
                w0_array[lo, la, press, wave] = w0_base[press][wave]
                '''
                 if (la >= lat_base[0] and la <= lat_base[1]):
                        opd_array[lo, la, press, wave] = opd_base[press][wave]
                        g0_array[lo, la, press, wave] = g0_base[press][wave]
                        w0_array[lo, la, press, wave] = w0_base[press][wave]
                 if (lo >= lon_spot[0]  and lo <= lon_spot[1]) and (la >= lat_spot[0] and la <= lat_spot[1]):
                        opd_array[lo, la, press, wave] += opd_4[press][wave]
                        g0_array[lo, la, press, wave] += g0_4[press][wave]
                        w0_array[lo, la, press, wave] += w0_4[press][wave]
                 elif not (lat_base[0] <= la <= lat_base[1]) and not (lon_spot[0] <= lo <= lon_spot[1] and lat_spot[0] <= la <= lat_spot[1]):
                        opd_array[lo, la, press, wave] = 0.0 #opd_4[press][wave]
                        g0_array[lo, la, press, wave] = 0.0 #g0_4[press][wave]
                        w0_array[lo, la, press, wave] = 0.0 #w0_4[press][wave]
                '''
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
#1D
neptune1d
#neptune1d.clouds(df=df_base.astype(float))
df_spec = neptune1d.spectrum(opacity)
wno_base, alb_base, fpfs_base = df_spec['wavenumber'] , df_spec['albedo'] , df_spec['fpfs_reflected']
wno1d, alb1d = pdi.mean_regrid(wno_base, alb_base, R=100)

ds_withchem= ds.update(ds_chem)
neptune3d.atmosphere_3d(ds_withchem,regrid=True,plot=False,verbose=True)
#neptune3d.clouds_3d(ds_cld,plot=True)
pdb.set_trace()
out3d = neptune3d.spectrum(opacity,calculation='reflected',dimension='3d',full_output=True)
wno3d,alb3d = pdi.mean_regrid(out3d['wavenumber'],out3d['albedo'],R=100)

plt.close()
plt.plot(1e4/wno3d,alb3d, label = '3d')
plt.plot(1e4/wno1d,alb1d, label = '1d')
plt.title("Neptune 1D vs 3D")
plt.xlabel("Wavelength (microns)")
plt.ylabel("Albedo")
plt.legend()
plt.show()
pdb.set_trace()