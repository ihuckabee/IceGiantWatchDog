import numpy as np
import numpy as np
import os
os.environ["picaso_refdata"]="/Users/ihuckabee/Documents/picaso-master/reference"
os.environ['PYSYN_CDBS']="/Users/ihuckabee/Documents/picaso-master/grp/redcat/trds"
refdata = os.getenv("picaso_refdata")
import pdb
from virga import justdoit as vdi
import pandas as pd
import astropy.units as u
from picaso import justplotit as ppi
from picaso import justdoit as pdi
import pandas as pd
#import xesmf as xe
import xarray as xr
from bokeh.plotting import show, figure
wave_range=[0.3,2.5]
opacity = pdi.opannection(wave_range)
opa = pdi.opannection(wave_range=wave_range)
gcm_out =  pdi.HJ_pt_3d()
#pdb.set_trace()
neptune = pdi.inputs()
neptune.phase_angle(0, num_gangle=3, num_tangle=3) #radians
#neptune.gravity(gravity=11.15, gravity_unit=pdi.u.Unit('m/(s**2)'), ) #any astropy units available
neptune.gravity(radius=0.3463884070945,radius_unit=pdi.u.Unit('R_jup'), mass=0.053740779768177, mass_unit=pdi.u.Unit('M_jup'))
neptune.star(opacity, temp=5778, metal=0.00, logg=4.4374, semi_major = 30.178, semi_major_unit = u.Unit("au"), radius=1.0, radius_unit=pdi.u.R_sun,)#opacity db, pysynphot database, temp, metallicity, logg
#need to create the xarray thing then feed it into atm3d 
#out3d = neptune.spectrum(opacity,calculation='reflected',dimension='3d',full_output=True)
#pdb.set_trace()
f = open("/Users/ihuckabee/Documents/picaso-master/nepforpicaso_salr.txt",'r')
lines = f.readlines()[1:]
pressure = []
temp = []
ch4 = []
h2s = []
for i in range(len(lines)):
        pressure.append(float(lines[i].split()[0]))
        temp.append(float(lines[i].split()[1]))
        ch4.append(float(lines[i].split()[2]))
        h2s.append(float(lines[i].split()[3]))
chem_other = 1 - np.array(ch4)+np.array(h2s)
h2 = 0.81*chem_other
he = 0.19*chem_other
pdb.set_trace()
pres = np.array(pressure) 
lon = np.linspace(-180,180,32)
lat = np.linspace(-90,90,16)
nwno = 24 
nep_data = np.tile(np.array(temp), (len(lon), len(lat), 1))
chem_ch4 = np.tile(np.array(ch4), (len(lon), len(lat), 1))
chem_h2s = np.tile(np.array(h2s), (len(lon), len(lat), 1))
chem_h2 = np.tile(np.array(h2), (len(lon), len(lat), 1))
chem_he = np.tile(np.array(he), (len(lon), len(lat), 1))
wno_grid = np.arange(0.2,2.6,0.1)


df_all = pd.read_csv('/Users/ihuckabee/Documents/picaso-master/cloudall.csv')
df_1 = pd.read_csv('/Users/ihuckabee/Documents/picaso-master/cloud1.csv')
df_2 = pd.read_csv('/Users/ihuckabee/Documents/picaso-master/cloud2.csv')
df_31 = pd.read_csv('/Users/ihuckabee/Documents/picaso-master/cloud31.csv')
df_4 = pd.read_csv('/Users/ihuckabee/Documents/picaso-master/cloud4.csv')
df_32 = pd.read_csv('/Users/ihuckabee/Documents/picaso-master/cloud32.csv')
df_base = df_1[["pressure", "wavenumber"]]
df_base['opd'] = df_1['opd']+df_2['opd']+df_31['opd']+df_32['opd']
df_base['g0'] = df_1['g0']+df_2['g0']+df_31['g0']+df_32['g0']
df_base['w0'] = df_1['w0']+df_2['w0']+df_31['w0']+df_32['w0']
opd_array = np.zeros((len(lon), len(lat), len(pres)-1, nwno)); g0_array = np.zeros((len(lon), len(lat), len(pres)-1, nwno)); w0_array = np.zeros((len(lon), len(lat), len(pres)-1, nwno))
# Iterate over each latitude and longitude point
opd_base = np.array(df_base['opd']).reshape(998,24)
g0_base= np.array(df_base['g0']).reshape(998,24)
w0_base = np.array(df_base['w0']).reshape(998,24)

opd_4 = np.array(df_4['opd']).reshape(998,24)
g0_4 = np.array(df_4['g0']).reshape(998,24)
w0_4 = np.array(df_4['w0']).reshape(998,24)

for lo in range(len(lon)):
    for la in range(len(lat)):
        for press in range(len(pres)-1):
            for wave in range(nwno):
                opd_array[lo, la, press, wave] = opd_base[press][wave]
                g0_array[lo, la, press, wave] = g0_base[press][wave]
                w0_array[lo, la, press, wave] = w0_base[press][wave]
            '''    if (lo >= 18) and (la >= 6 and la <= 10):
                    opd_array[lo, la, press, wave] = opd_4[press][wave]
                    g0_array[lo, la, press, wave] = g0_4[press][wave]
                    w0_array[lo, la, press, wave] = w0_4[press][wave]
                else:
                    opd_array[lo, la, press, wave] = opd_base[press][wave]
                    g0_array[lo, la, press, wave] = g0_base[press][wave]
                    w0_array[lo, la, press, wave] = w0_base[press][wave]
                    '''        
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
        wno=(["wno"], 1/(wno_grid*1e-4),{'units': 'cm^(-1)'})#required for clouds NEEDS CHEWCKING FOR UNITS AND WHATENOT. what do the 24 wv points map to
    ),
    attrs=dict(description="coords with vectors"),
)
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
all_gcm = ds.update(ds_chem)
#pdb.set_trace()
#n_gauss_angles =3
#n_chebychev_angles=3

#gangle, gweight, tangle, tweight = pdi.get_angles_3d(n_gauss_angles, n_chebychev_angles)
#ubar0, ubar1, cos_theta, latitude, longitude = pdi.compute_disco(n_gauss_angles, n_chebychev_angles, gangle, tangle, phase_angle=0)
neptune.atmosphere_3d(all_gcm,regrid=True,plot=False,verbose=True)
neptune.chemeq_3d(c_o = 0.55, n_cpu=5)
neptune.clouds_3d(ds_cld)
#problems start below
#pdb.set_trace()

out3d = neptune.spectrum(opacity,calculation='reflected',dimension='3d',full_output=True)
#pdb.set_trace()
wno,fpfs = pdi.mean_regrid(out3d['wavenumber'],out3d['fpfs_reflected'],R=100)
sm = 4.514565e12 #m
rp = 2.47639e7 #m
alb = fpfs*(sm/rp)**2
np.save("3dcloudbase_spec.npy", np.array([wno,alb,fpfs]))
#ppi.show(ppi.spectrum(wno, alb, plot_width=500,y_axis_type='log'))
pdb.set_trace()
ppi.show(ppi.disco(out3d['full_output'] ,wavelength=[1.1, 1.4], calculation='reflected'))
pdb.set_trace()
