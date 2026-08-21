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
from bokeh.plotting import show, figure
wave_range=[0.3,2.5]
opacity = pdi.opannection(wave_range)
opa = pdi.opannection(wave_range=wave_range)

neptune = pdi.inputs()
neptune.phase_angle(0) #radians
neptune.gravity(gravity=11.15, gravity_unit=pdi.u.Unit('m/(s**2)')) #any astropy units available
neptune.star(opacity, temp=5778, metal=0.00, logg=4.4374, radius=1.0, radius_unit=pdi.u.R_sun )#opacity db, pysynphot database, temp, metallicity, logg
neptune.atmosphere(filename="/Users/ihuckabee/Documents/picaso-master/nepforpicaso_salr.txt", delim_whitespace=True)


#case1.phase_angle(num_tangle=10, num_gangle=10)
neptune.phase_angle(num_tangle=10, num_gangle=10)
#neptune.phase_angle(num_tangle=1, num_gangle=10) #remember num_tangle=1 automatically insinuates symmetry
df = neptune.spectrum(opacity, full_output = True)
x,y = df['wavenumber'], df['albedo'] #units of erg/cm2/s/cm
xmicron = 1e4/x
flamy = y*1e-8 #per anstrom instead of per cm
sp = pdi.psyn.ArraySpectrum(xmicron, flamy,waveunits='um',fluxunits='FLAM')
sp.convert("um")
sp.convert('Fnu') #erg/cm2/s/Hz
x = sp.wave #micron
y= sp.flux #erg/cm2/s/Hz
df['fluxnu'] = y
x,y = pdi.mean_regrid(1e4/x, y, R=300) #wavenumber, erg/cm2/s/Hz
df['regridy'] =  y
df['regridx'] = x

asdict = df['full_output']
#ppi.show(ppi.disco(asdict, calculation='reflected', wavelength=[0.5]))
#pdb.set_trace()
show(ppi.disco(asdict, calculation='reflected', wavelength=[1.5]))
