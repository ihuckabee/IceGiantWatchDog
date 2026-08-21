import os
os.environ["picaso_refdata"]="/Users/ihuckabee/Documents/picaso-master/reference"
os.environ['PYSYN_CDBS']="/Users/ihuckabee/Documents/picaso-master/grp/redcat/trds"
import numpy as np

import pdb
refdata = os.getenv("picaso_refdata")
from picaso import justdoit as jdi
from picaso import justplotit as jpi
#from bokeh.plotting import show, figure
opacity = jdi.opannection(wave_range=[0.3,2.5])
neptune = jdi.inputs()
#phase angle
neptune.phase_angle(0) #radians
#define gravity
neptune.gravity(gravity=11.15, gravity_unit=jdi.u.Unit('m/(s**2)')) #any astropy units available

#define star
pdb.set_trace()
neptune.star(opacity, temp=5778, metal=0.00, logg=4.4374, radius=1.0, radius_unit=jdi.u.R_sun )#opacity db, pysynphot database, temp, metallicity, logg
neptune.atmosphere(filename="/Users/ihuckabee/Documents/picaso-master/nepforpicaso_salr.txt", delim_whitespace=True)
comp_file = jdi.pd.read_csv("/Users/ihuckabee/Documents/picaso-master/nepforpicaso_salr.txt", delim_whitespace=True)
df = neptune.spectrum(opacity)
wno, alb, fpfs = df['wavenumber'] , df['albedo'] , df['fpfs_reflected']
wno, alb = jdi.mean_regrid(wno, alb , R=100)
fig = jpi.spectrum(wno, alb, plot_width=500)
fig.line(1e4/wno, alb, line_width=2, color='blue')
#jpi.show(fig)
#jpi.close()
'''
neptune.atmosphere(filename="/Users/ihuckabee/Documents/picaso-master/nepforpicaso_salr.txt", delim_whitespace=True)
comp_file = jdi.pd.read_csv("/Users/ihuckabee/Documents/picaso-master/nepforpicaso_salr.txt", delim_whitespace=True)
df = neptune.spectrum(opacity)
wno, alb, fpfs = df['wavenumber'] , df['albedo'] , df['fpfs_reflected']
wno2, alb2 = jdi.mean_regrid(wno, alb , R=100)
fig = jpi.spectrum(wno2, alb2, plot_width=500)
fig.line(1e4/wno2, alb2, line_width=2, color='blue', label = 'salr')
'''
import pandas as pd
nwno = 42 #this is the default number for A&M cloud code (see below if your wave grid is different)
nlayer = len(comp_file['pressure'])-1 #one less than the number of PT points in your input
#cld = pd.read_csv("/Users/ihuckabee/Documents/picaso-master/nepcloud_irwin21.txt", delim_whitespace=True)
pdb.set_trace()
neptune.clouds(filename = "/Users/ihuckabee/Documents/picaso-master/nepcloud_irwin21.cld", delim_whitespace=True)
#neptune.clouds(g0=[0.8], w0=[1.], opd=[1.0], p = [np.log(1.5)], dp=[0.22])
#^find the right dp. rn it seems that p = log(loc of slab (bars)), so = 0.0 = log(1bar))
#nwno = 196 #this is the default number for A&M cloud code (see below if your wave grid is different)
nlayer = len(comp_file['pressure'])-1 #one less than the number of PT points in your input
#jpi.show(jpi.plot_cld_input(nwno, nlayer,df=neptune.inputs['clouds']['profile']))
df = neptune.spectrum(opacity)
wno_c, alb_c, fpfs_c = df['wavenumber'] , df['albedo'] , df['fpfs_reflected']
wno_c, alb_c = jdi.mean_regrid(wno_c, alb_c, R=100)
#pdb.set_trace()
jpi.show(jpi.spectrum([wno, wno_c], [alb,alb_c],legend = ["Cloud-free", "Irwin Cloud"], plot_width=750))
pdb.set_trace()