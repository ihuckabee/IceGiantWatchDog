import pandas as pd
import numpy as np
import os
os.environ["picaso_refdata"]="/Users/ihuckabee/Documents/picaso-master/reference"
os.environ['PYSYN_CDBS']="/Users/ihuckabee/Documents/picaso-master/grp/redcat/trds"
refdata = os.getenv("picaso_refdata")
from picaso import justdoit as jdi
from picaso import justplotit as jpi
from virga import justdoit as vdi
import pdb
import xarray as xr
import astropy.units as u
import matplotlib.pyplot as plt
opacity = jdi.opannection(wave_range=[0.3,2.5])

#1D SET UP W/ CLOUDS
case1d = jdi.inputs()
case1d.phase_angle(0) #radians
case1d.gravity(radius=1,radius_unit=jdi.u.Unit('R_jup'), mass=1, mass_unit=jdi.u.Unit('M_jup')) #any astropy units available
case1d.star(opacity, temp=5778, metal=0.00, logg=4.4374, semi_major = 30.178, semi_major_unit = u.Unit("au"), radius=1.0, radius_unit=jdi.u.R_sun,)
#atmosphere
case1d.atmosphere(filename= jdi.HJ_pt(), delim_whitespace=True)
#case1d.clouds(filename= jdi.HJ_cld(), delim_whitespace=True)
#jpi.output_notebook()
#gcm_out = jdi.HJ_pt_3d(as_xarray=True)
pt_1d = pd.read_csv(jdi.HJ_pt(), delim_whitespace=True)
#pt_3d = pd.read_csv(jdi.HJ_pt_3d(as_xarray=True), delim_whitespace=True)
# create coords
lon = np.linspace(-180,180,128)
lat = np.linspace(-90,90,64)
pres_1d = np.array(pt_1d['pressure'])
temp = np.array(pt_1d['temperature'])
#pdb.set_trace()
h2_1D = pt_1d['H2']; h_1D = pt_1d['H']; he_1D = pt_1d['He']; h2o_1D = pt_1d['H2O']; ch4_1D = pt_1d['CH4']; co_1D = pt_1d['CO']; nh3_1D = pt_1d['NH3']; n2_1D = pt_1d['N2']; h2s_1D = pt_1d['H2S'];fe_1D = pt_1d['Fe']; na_1D = pt_1d['Na']; k_1D = pt_1d['K']; co2_1D = pt_1d['CO2']; sio_1D = pt_1d['SiO'];  
for i in range(len(temp)):
        chemsum = h2_1D[i] + h_1D[i]+ he_1D[i] + h2o_1D[i]+ co_1D[i]+ h2s_1D[i]+ ch4_1D[i] + fe_1D[i]+ na_1D[i]+ k_1D[i]+ co2_1D[i] +  nh3_1D[i]
        remainder = 1-chemsum
        h2_1D[i] += remainder
h2_3D = np.tile(np.array(h2_1D), (len(lon), len(lat), 1))
h_3D = np.tile(np.array(h_1D), (len(lon), len(lat), 1))
he_3D = np.tile(np.array(he_1D), (len(lon), len(lat), 1))
h2o_3D = np.tile(np.array(h2o_1D), (len(lon), len(lat), 1))
co_3D = np.tile(np.array(co_1D), (len(lon), len(lat), 1))
h2s_3D = np.tile(np.array(h2s_1D), (len(lon), len(lat), 1))
ch4_3D = np.tile(np.array(ch4_1D), (len(lon), len(lat), 1))
fe_3D = np.tile(np.array(fe_1D), (len(lon), len(lat), 1))
na_3D = np.tile(np.array(na_1D), (len(lon), len(lat), 1))
k_3D = np.tile(np.array(k_1D), (len(lon), len(lat), 1))
co2_3D = np.tile(np.array(co2_1D), (len(lon), len(lat), 1))
temp_data = np.tile(np.array(temp), (len(lon), len(lat), 1))
#fake_chem_H2O = np.random.rand(len(lon), len(lat),len(pres_1d))*0.1+0.1 # create fake data
#fake_chem_H2 = 1-fake_chem_H2O # create data
#pdb.set_trace()
# put data into a dataset
ds = xr.Dataset(
    data_vars=dict(
        temperature=(["lon", "lat","pressure"], temp_data,{'units': 'Kelvin'})#, required
        #kzz = (["x", "y","z"], gcm_out['kzz'])#could add other data components if wanted
    ),
    coords=dict(
        lon=(["lon"], lon,{'units': 'degrees'}),#required
        lat=(["lat"], lat,{'units': 'degrees'}),#required
        pressure=(["pressure"], pres_1d,{'units': 'bar'})#required*
    ),
    attrs=dict(description="coords with vectors"),
)
ds_chem = xr.Dataset(
    data_vars=dict(
        H2O=(["lon", "lat","pressure"], h2o_3D,{'units': 'v/v'}),
        H2=(["lon", "lat","pressure"], h2_3D,{'units': 'v/v'}),
        H=(["lon", "lat","pressure"], h_3D,{'units': 'v/v'}),
        He=(["lon", "lat","pressure"], he_3D,{'units': 'v/v'}),
        CO=(["lon", "lat","pressure"], co_3D,{'units': 'v/v'}),
        H2S=(["lon", "lat","pressure"], h2s_3D,{'units': 'v/v'}),
        Na=(["lon", "lat","pressure"], na_3D,{'units': 'v/v'}),
        K=(["lon", "lat","pressure"], k_3D,{'units': 'v/v'}),
        Fe=(["lon", "lat","pressure"], fe_3D,{'units': 'v/v'}),
        CH4=(["lon", "lat","pressure"], ch4_3D,{'units': 'v/v'}),
        CO2=(["lon", "lat","pressure"], co2_3D,{'units': 'v/v'}),
    ),
    coords=dict(
        lon=(["lon"], lon,{'units': 'degrees'}),#required
        lat=(["lat"], lat,{'units': 'degrees'}),#required
        pressure=(["pressure"], pres_1d,{'units': 'bar'})#required*
    ),
    attrs=dict(description="coords with vectors"),
)

all_gcm = ds.update(ds_chem)
# create coords
wno_grid = np.linspace(1e4/0.3,1e4/2.5/1,10)#cloud properties are defined on a wavenumber grid

#create box-band cloud model
fake_opd = np.zeros((len(lon), len(lat),len(pres_1d)-1, len(wno_grid))) # create fake data
fake_asymmetry_g0 = np.zeros((len(lon), len(lat),len(pres_1d)-1, len(wno_grid)))
fake_ssa_w0 =  np.zeros((len(lon), len(lat),len(pres_1d)-1, len(wno_grid)))
where_pres = np.where(((pres_1d<0.01) & (pres_1d>.001)))#creating a grey cloud band
for ip in where_pres[0]:
    fake_opd[:,:,ip,:]= 10 #optical depth of 10 (>>1)
    fake_asymmetry_g0[:,:,ip,:]=0.8
    fake_ssa_w0[:,:,ip,:]=0.9
#make up asymmetry and single scattering properties
#this is what im playing around with rn. making the 1d and 3d example clouds the same
df_cld= vdi.picaso_format_slab(0.01, 10+np.zeros((len(wno_grid))), 0.8+np.zeros((len(wno_grid))), 0.9+np.zeros((len(wno_grid))), wno_grid, pres_1d[:-1],p_top=0.001)
case1d.clouds(df=df_cld.astype(float))
#pdb.set_trace()
# put data into a dataset
ds_cld= xr.Dataset(
    data_vars=dict(
        opd=(["lon", "lat","pressure","wno"], fake_opd,{'units': 'depth per layer'}),
        g0=(["lon", "lat","pressure","wno"], fake_asymmetry_g0,{'units': 'none'}),
        w0=(["lon", "lat","pressure","wno"], fake_ssa_w0,{'units': 'none'}),
    ),
    coords=dict(
        lon=(["lon"], lon,{'units': 'degrees'}),#required
        lat=(["lat"], lat,{'units': 'degrees'}),#required
        pressure=(["pressure"], pres_1d[:-1],{'units': 'bar'}),#required
        wno=(["wno"], wno_grid,{'units': 'cm^(-1)'})#required for clouds
    ),
    attrs=dict(description="coords with vectors"),
)
#first step is identical to what's been done in the past
case_3d = jdi.inputs()
case_3d.phase_angle(0, num_tangle=3, num_gangle=3)
#turning off alerts since we already went through this
case_3d.atmosphere_3d(all_gcm, regrid=True, plot=False,verbose=False)

case_3d.clouds_3d(ds_cld,plot=True) #investigate this 
#pdb.set_trace()
case_3d.gravity(radius=1,radius_unit=jdi.u.Unit('R_jup'),
                mass=1, mass_unit=jdi.u.Unit('M_jup')) #any astropy units available
case_3d.star(opacity, temp=5778, metal=0.00, logg=4.4374, semi_major = 0.1, semi_major_unit = u.Unit("au"), radius=1.0, radius_unit=jdi.u.R_sun,)
out3d_cld = case_3d.spectrum(opacity,calculation='reflected',dimension='3d',full_output=True)
out1d_cld = case1d.spectrum(opacity, calculation='reflected')
#pdb.set_trace()
wno_3d,alb_3d = jdi.mean_regrid(out3d_cld['wavenumber'],out3d_cld['albedo'],R=100)
wno_1d, alb_1d = jdi.mean_regrid(out1d_cld['wavenumber'] , out1d_cld['albedo'],R=100)
#pdb.set_trace()
plt.close()
plt.plot(1e4/wno_3d,alb_3d, label = '3d')
plt.plot(1e4/wno_1d,alb_1d, label = '1d')
plt.title("PICASO Hot Jupiter w/ cloud (opd = 10)")
plt.xlabel("Wavelength (microns)")
plt.ylabel("Albedo")
plt.legend()
plt.show()
pdb.set_trace()