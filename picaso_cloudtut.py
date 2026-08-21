import os
os.environ["picaso_refdata"]="/Users/ihuckabee/Documents/picaso-master/reference"
os.environ['PYSYN_CDBS']="/Users/ihuckabee/Documents/picaso-master/grp/redcat/trds"
import numpy as np
import os
refdata = os.getenv("picaso_refdata")
import numpy as np
import pandas as pd
import astropy.units as u

#picaso & virga
from picaso import justdoit as pj
from virga import justdoit as vj
#plot tools
from picaso import justplotit as picplt
from virga import justplotit as cldplt

from bokeh.plotting import show, figure
from bokeh.io import output_notebook
#output_notebook()

opacity = pj.opannection(wave_range=[0.3,1])
sum_planet = pj.inputs()
sum_planet.phase_angle(0) #radians
sum_planet.gravity(gravity=25, gravity_unit=u.Unit('m/(s**2)')) #any astropy units available
sum_planet.star(opacity, 5000,0,4.0) #opacity db, pysynphot database, temp, metallicity, logg

df_atmo = pd.read_csv(pj.jupiter_pt(), sep='\s+')
#you will have to add kz to the picaso profile
df_atmo['kz'] = [1e9]*df_atmo.shape[0]

#business as usual
sum_planet.atmosphere(df=df_atmo)

#let's get the cloud free spectrum for reference
cloud_free = sum_planet.spectrum(opacity)

x_cld_free, y_cld_free = pj.mean_regrid(cloud_free['wavenumber'], cloud_free['albedo'], R=150)

metallicity = 1 #atmospheric metallicity relative to Solar
mean_molecular_weight = 2.2 # atmospheric mean molecular weight
directory ='/Users/ihuckabee/Documents/picaso-master/data/virga/'
#import pdb; pdb.set_trace()
#we can get the same full output from the virga run
cld_out = sum_planet.virga(['H2O'],directory, fsed=1,mh=metallicity,
                 mmw = mean_molecular_weight)
out = sum_planet.spectrum(opacity, full_output=True)

x_cldy, y_cldy = pj.mean_regrid(out['wavenumber'], out['albedo'], R=150)

show(picplt.spectrum([x_cld_free, x_cldy],
                     [y_cld_free, y_cldy],plot_width=500, plot_height=300,
                  legend=['Cloud Free','Cloudy']))

#show(picplt.photon_attenuation(out['full_output'],
                            #plot_width=500, plot_height=300))

#fig, dndr = cldplt.radii(cld_out,at_pressure=0.1)
#show(fig)                            
import pdb; pdb.set_trace()
#grab your mie parameters
qext, qscat, g_qscat, nwave,radii,wave = vj.get_mie('H2O',directory)

from bokeh.layouts import row,column
ind = cldplt.find_nearest_1d(radii,30e-4) #remember the radii are in cm

qfig = figure(width=300, height=300,
              x_axis_type='log',y_axis_label ='Asymmetry',
              x_axis_label='Wavelength(um)')

wfig = figure(width=300, height=300,
              x_axis_type='log',y_axis_label ='Qscat/Qext',
              x_axis_label='Wavelength(um)')

qfig.line(1e4*wave[:,0], g_qscat[:,ind]/qscat[:,ind])
wfig.line(1e4*wave[:,0], qscat[:,ind]/qext[:,ind])

show(row(qfig, wfig))

hot_atmo = df_atmo
hot_atmo['temperature'] = hot_atmo['temperature'] + 600

#remember we can use recommend_gas function to look at what the condensation curves look like
recommended = vj.recommend_gas(hot_atmo['pressure'], hot_atmo['temperature'], metallicity,mean_molecular_weight,
                #Turn on plotting and add kwargs for bokeh.figure
                 plot=True, y_axis_type='log',y_range=[1e2,1e-3],
                               plot_height=400, plot_width=600,
                  y_axis_label='Pressure(bars)',x_axis_label='Temperature (K)')

#business as usual
sum_planet.atmosphere(df=hot_atmo)

#make sure clouds are turned off
sum_planet.clouds_reset()

#let's get the cloud free spectrum for reference
cloud_free = sum_planet.spectrum(opacity)
x_cld_free, y_cld_free = pj.mean_regrid(cloud_free['wavenumber'], cloud_free['albedo'], R=150)

#now the cloudy runs
cld_out = sum_planet.virga(['Na2S','ZnS'],directory, fsed=1,mh=metallicity,
                 mmw = mean_molecular_weight)

out = sum_planet.spectrum(opacity, full_output=True)
x_cld, y_cld = pj.mean_regrid(out['wavenumber'], out['albedo'], R=150)

w = [x_cld_free, x_cld]
a = [y_cld_free, y_cld]
show(picplt.spectrum(w,a,plot_width=500, plot_height=300,
                  legend=['Cloud Free','Cloudy']))

show(picplt.photon_attenuation(out['full_output'],
                            plot_width=500, plot_height=300))

fig, dndr = cldplt.radii(cld_out,at_pressure=0.5)
show(fig)

#grab your mie parameters
gas_name = 'Na2S' #ZnS
qext, qscat, g_qscat, nwave,radii,wave = vj.get_mie(gas_name,directory)
ind = cldplt.find_nearest_1d(radii,10e-4) #remember the radii are in cm

qfig = figure(width=300, height=300,
              x_axis_type='log',y_axis_label ='Asymmetry',
              x_axis_label='Wavelength(um)')

wfig = figure(width=300, height=300,
              x_axis_type='log',y_axis_label ='Qscat/Qext',
              x_axis_label='Wavelength(um)')

qfig.line(1e4*wave[:,0], g_qscat[:,ind]/qscat[:,ind])
wfig.line(1e4*wave[:,0], qscat[:,ind]/qext[:,ind])

show(row(qfig, wfig))

df_atmo = pd.read_csv(pj.jupiter_pt(), delim_whitespace=True)
df_atmo['kz'] = [1e10]*df_atmo.shape[0]

sum_planet.atmosphere(df = df_atmo)

all_fseds =  [1, 6, 10]
w = []
a = []
all_outs,cld_outs=[],[]
for fs in all_fseds:
    cld = sum_planet.virga(['H2O'],directory, fsed=fs,mh=metallicity,
                 mmw = mean_molecular_weight)
    cld_outs += [cld]
    out = sum_planet.spectrum(opacity,full_output=True)
    x,y = pj.mean_regrid(out['wavenumber'], out['albedo'], R=150)
    w += [x]
    a += [y]
    all_outs += [out['full_output']]

show(picplt.spectrum(w,a,plot_width=500, plot_height=300,
                     legend=['fs= '+str(i) for i in all_fseds]))

show(column(picplt.photon_attenuation(all_outs[0],title='kz=1e6',
                            plot_width=500, plot_height=300),
   picplt.photon_attenuation(all_outs[2],title='kz=1e10',
                            plot_width=500, plot_height=300)))

fig, dndr = cldplt.radii(cld_outs[0],at_pressure=1e-2)
print('kz=1e6')
show(fig)

fig, dndr = cldplt.radii(cld_outs[2],at_pressure=1e-1)
print('kz=1e10')
show(fig)