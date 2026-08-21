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
#loading in custom radii, custom refractive indicies to get compute mie.. 
#first need to read in refrind and add "index", then spit it back out into a refrind file
#its a whole thing
#index, wavelength, real, imaginary

#reading in pt profile
comp_file = pdi.pd.read_csv("/Users/ihuckabee/Documents/picaso-master/neptune_stuff/ToShare/RefFiles/nepforpicaso_salr_abridged.txt", delim_whitespace=True)


#refractive indices files - from elijah 

refrac_files = ['/Users/ihuckabee/Downloads/refindex_methane.txt','/Users/ihuckabee/Downloads/refindex_neptune_sample_aerosol1.txt','/Users/ihuckabee/Downloads/refindex_neptune_sample_aerosol2.txt','/Users/ihuckabee/Downloads/refindex_neptune_sample_aerosol3.txt', '/Users/ihuckabee/Documents/picaso-master/neptune_stuff/refindex_aerosol1mod_irwin.txt']
a_names = ["4","1","2","3", "1mod"]
#reformatting the refractive indices files so picaso likes them 
#(only need to do this step once so i commented it out for not)
'''
for k in range(len(refrac_files)):
    wv = []; real = []; im = []; index = []
    f = open(refrac_files[k],'r')
    lines = f.readlines()
    for i in range(len(lines)-1):
        #wv.append(float(lines[i].split()[0]))
        try:
            wv.append(float(lines[i].split()[0]))
            real.append(lines[i].split()[1])
            im.append(lines[i].split()[2])
            #im.append(float(lines[i].split()[2]))
            index.append(str(i))
        except:
            continue
    f.close()

    #spit back out
    file1 = open("/Users/ihuckabee/Documents/picaso-master/picaso/refrind/test_"+a_names[k]+".refrind", "w")
    #i = 0
    #import pdb; pdb.set_trace()
    for v in range(len(lines)):
        try:
            file1.write(index[v]+" "+str(wv[v])+" "+real[v]+" "+im[v]+" \n")
        except:
            continue
    file1.close()
import pdb; pdb.set_trace()
'''
#read in radii distributions that we previously created
f = open("/Users/ihuckabee/Documents/picaso-master/neptune_stuff/nepcustomradii.txt",'r') #LOCAL
#f = open("/Users/ihuckabee/Documents/picaso-master/nepcustomradii.txt",'r') #ASH
lines = f.readlines()[1:]
i = 0; rsize = []; r13 = []; r2 = []; r4 = []; 
for i in range(len(lines)):
    rsize.append(float(lines[i].split()[0]))
    r13.append(float(lines[i].split()[1]))
    r2.append(float(lines[i].split()[2]))
    r4.append(float(lines[i].split()[4]))
f.close()
r13 = np.array(r13) 
r2 = np.array(r2) 
r4 = np.array(r4) 
rsize = np.array(rsize) 


#r is the grid of radius bins, n is the distribution per particle per radius bin 
r_arr = np.array([rsize.T,r13.T,r2.T,r4.T]).T
df = pd.DataFrame(r_arr, columns = ['size','146K','84K', '54.7K'])
non_nans_13 = np.where(np.isfinite(df['146K']))
r13 = df['size'].iloc[non_nans_13]
n13 = df['146K'].iloc[non_nans_13]

non_nans_2 = np.where(np.isfinite(df['84K']))
r2 = df['size'].iloc[non_nans_2]
n2 = df['84K'].iloc[non_nans_2]

non_nans_4 = np.where(np.isfinite(df['54.7K']))
r4 = df['size'].iloc[non_nans_4]
n4 = df['54.7K'].iloc[non_nans_4]

r13 = np.array(r13)
r2 = np.array(r2)  
r4 = np.array(r4)     




mieff_dir = '/Users/ihuckabee/Documents/picaso-master/picaso/refrind'
# here I am assuming the refractive indices are stored in the same place as I want the output
refrind_dir = output_dir = mieff_dir


#Note that this function wants the radii in centimeters, so we have to convert it from nanometers
#vv this step only needs to be done once

#newmies_1mod=vdi.calc_mie_db('formattedrefrac_'+a_names[-1], refrind_dir, output_dir, rmin = np.min(r13/1e7), rmax=np.max(r13/1e7), nradii = len(r13))
#newmies_1=vdi.calc_mie_db('formattedrefrac_'+a_names[1], refrind_dir, output_dir, rmin = np.min(r13/1e7), rmax=np.max(r13/1e7), nradii = len(r13))
#newmies_3=vdi.calc_mie_db('formattedrefrac_'+a_names[3], refrind_dir, output_dir, rmin = np.min(r13/1e7), rmax=np.max(r13/1e7), nradii = len(r13))
#newmies_4=vdi.calc_mie_db('formattedrefrac_'+a_names[0], refrind_dir, output_dir, rmin = np.min(r4/1e7), rmax=np.max(r4/1e7), nradii = len(r4))
#newmies_2=vdi.calc_mie_db('formattedrefrac_'+a_names[2], refrind_dir, output_dir, rmin = np.min(r2/1e7), rmax=np.max(r2/1e7), nradii = len(r2))
#newmies_2=vdi.calc_mie_db("ASHrefrac_2", refrind_dir, output_dir, rmin = np.min(r2/1e7), rmax=np.max(r2/1e7), nradii = len(r2))
#newmies_1=vdi.calc_mie_db("ASHrefrac_1", refrind_dir, output_dir, rmin = np.min(r13/1e7), rmax=np.max(r13/1e7), nradii = len(r13))

#**** FOR OUT TO 2.5 MICRONS ***
qext_1, qscat_1, cos_qscat_1, nwave_1, radius_1, wave_in_1 = vdi.get_mie("formattedrefrac_"+a_names[1],directory=mieff_dir)
#qext_1, qscat_1, cos_qscat_1, nwave_1, radius_1, wave_in_1 = vdi.get_mie("formattedrefrac_"+a_names[1],directory=mieff_dir)
#qext_1_ash, qscat_1_ash, cos_qscat_1_ash, nwave_1_ash, radius_1_ash, wave_in_1_ash = vdi.get_mie("ASHrefrac_1",directory=mieff_dir)

qext_3, qscat_3, cos_qscat_3, nwave_3, radius_3, wave_in_3 = vdi.get_mie("formattedrefrac_"+a_names[3],directory=mieff_dir)
qext_4, qscat_4, cos_qscat_4, nwave_4, radius_4, wave_in_4 = vdi.get_mie("formattedrefrac_"+a_names[0],directory=mieff_dir)
#qext_2, qscat_2, cos_qscat_2, nwave_2, radius_2, wave_in_2 = vdi.get_mie("ASHrefrac_2",directory=mieff_dir)
qext_2, qscat_2, cos_qscat_2, nwave_2, radius_2, wave_in_2 = vdi.get_mie("formattedrefrac_"+a_names[2],directory=mieff_dir)
qext_1mod, qscat_1mod, cos_qscat_1mod, nwave_1mod, radius_1mod, wave_in_1mod = vdi.get_mie("formattedrefrac_"+a_names[-1],directory=mieff_dir)

#**** FOR TRUNCATED ****
'''
qext_1, qscat_1, cos_qscat_1, nwave_1, radius_1, wave_in_1 = vdi.get_mie("formattedrefrac_"+a_names[1]+'_trunc',directory=mieff_dir)
qext_3, qscat_3, cos_qscat_3, nwave_3, radius_3, wave_in_3 = vdi.get_mie("formattedrefrac_"+a_names[3]+'_trunc',directory=mieff_dir)
qext_4, qscat_4, cos_qscat_4, nwave_4, radius_4, wave_in_4 = vdi.get_mie("formattedrefrac_"+a_names[0]+'_trunc',directory=mieff_dir)
qext_2, qscat_2, cos_qscat_2, nwave_2, radius_2, wave_in_2 = vdi.get_mie("formattedrefrac_"+a_names[2]+'_trunc',directory=mieff_dir)
'''

#now we gotta calculate ndz. the picaso tutorials just make it up 
#but i was given a tau @ 0.8 microns so i tried to back out the ndz 
#but it made me go insane
#so what i did was use the cloud opacity/bar plot in the paper to set ndz arbitrarily
#to get the opd i wanted 

meanr13 = 0.05/1e4 #mean radii from irwin paper given in microns
#meanr2 = 0.55/1e4 #converted to cm 
meanr2 = 0.68/1e4 #converted to cm 
meanr4 = 2.5/1e4
#i dont know what i was trying to do with the mean_index stuff below
#radius is natively given in nanometers (i.e. r13 is in nanometers)
#tau is given at 0.8 micrometers so calculating sig_ext at 0.8 to get the ndz 
#find locations of means in the r13 area. 

mean_indexr13 = 53
mean_indexr2 = 553
mean_indexr4 = 2500

#pressure limits of clouds based on the paper
#note that i will technically have 5 separate clouds because cloud 4 is sandiwiched in cloud 3
base_pressure_1 = 10. 
haze_top_pressure_1 = 2.05 
base_pressure_2 = 2.05
haze_top_pressure_2 = 1.6
base_pressure_31 = 1.6
haze_top_pressure_31 = 0.2
base_pressure_4 = 0.2
haze_top_pressure_4 = 0.08
base_pressure_32 = 0.08
haze_top_pressure_32 = 0.01


#all the stuff in the starred area you technically dont need to do
#it didnt even work i still had to arbitrarily scale stuff
#so just guess and check from the get go 
#*********************

#to calculate scale height 
mu = 2.6
mp = 1.6726219e-24
g_cgs = 11.15*1e2 #cm/s^2
boltz = 1.38065e-16 #erg/K

#estimating temps
temp_cloud1 = 127.5 #K
temp_cloud2 = 93.5 #K
temp_cloud3 = 60. #K
temp_cloud4 = 54.1995 #K 
cloudtemps = [temp_cloud1, temp_cloud2, temp_cloud3, temp_cloud4]
cloudpressure = [[base_pressure_1, haze_top_pressure_1], [base_pressure_2, haze_top_pressure_2], [base_pressure_31, haze_top_pressure_31], [base_pressure_4, haze_top_pressure_4]]
stepsize = [5, 1, 6, 3]
z = []
#calculating scale heights 
for i in range(len(cloudtemps)):
    h = boltz*cloudtemps[i]/(mu*g_cgs*mp) 
    z.append(-h*np.log(cloudpressure[i][1]/cloudpressure[i][0])/stepsize[i])

#************************


tau1 = 0.65; tau2 = 1.5; tau3 = 0.04; tau4 = 0.03 #mean taus
#tau1 = 0.8; tau2 = 2.; tau3 = 0.04; tau4 = 0.03 #max taus 
#tau1 = 0.5; tau2 = 1.; tau3 = 0.04; tau4 = 0.03 #min taus 

#ndz: that's particles/cm^2 but over the entire region of atmosphere where we're sticking aerosol

#ndzs feed into the calc_optics_user_r_dist stuff
#THESE NDZ DONT WORK EITHER
#sig_ext1 = np.pi*meanr13**2*qext_1[5][mean_indexr13]; ndz_1 = 6.92e5*(tau1/sig_ext1)/z[0] 
#sig_ext1 = np.pi*meanr13**2*qext_1[-7][mean_indexr13]; ndz_1 = 7e5*(tau1/sig_ext1)/z[0] 
#sig_ext2 = np.pi*meanr2**2*qext_2[-7][mean_indexr2] ; ndz_2= 4.*(tau2/sig_ext2)/z[1]
#sig_ext3 = np.pi*meanr13**2*qext_3[-7][mean_indexr13]; ndz_31 = 6.75e1*(tau3/sig_ext3)/z[2]; ndz_32 = 6.75e1*(tau3/sig_ext3)/z[2]
#sig_ext4 = np.pi*meanr4**2*qext_4[-7][mean_indexr4]; ndz_4 = 7.15e2*(tau4/sig_ext4)/z[3]
#ndz_1mod = (tau1/sig_ext1)/z[0] 

sig_ext1 = np.pi*meanr13**2*qext_1[-7][mean_indexr13]; ndz_1 = 2.75e15
sig_ext2 = np.pi*meanr2**2*qext_2[-7][mean_indexr2] ; ndz_2= 1.e12
sig_ext3 = np.pi*meanr13**2*qext_3[-7][mean_indexr13]; ndz_31 = 1.75e12; ndz_32 = 1.75e12
sig_ext4 = np.pi*meanr4**2*qext_4[-7][mean_indexr4]; ndz_4 = 7.3e8
ndz_1mod = 2.25e14
'''
#NEW SHIT THAT DIDNT WORK
idx_08um = np.where(wave_in_1.T[0] == 0.8)[0][0]
sig_ext1 = np.pi*meanr13**2*qext_1[idx_08um][mean_indexr13]; ndz_1 = 4.6e3*(tau1/sig_ext1)/z[0] 
sig_ext2 = np.pi*meanr2**2*qext_2[idx_08um][mean_indexr2] ; ndz_2= 7.15*(tau2/sig_ext2)/z[1]
sig_ext3 = np.pi*meanr13**2*qext_3[idx_08um][mean_indexr13]; ndz_3 = 6.75e1*(tau3/sig_ext3)/z[2]
sig_ext4 = np.pi*meanr4**2*qext_4[idx_08um][mean_indexr4]; ndz_4 = 1.9*(tau4/sig_ext4)/z[3]
ndz_31 = 6.1e1*(tau3/sig_ext3)/z[2]; ndz_32 = 3.8*(tau3/sig_ext3)/z[2]

#why index 5? to correlate with 0.8um? But 5 doesn't... 
#what unit is sig in? has to cm^2. 

#im given sig as a function of wavelength, i need to scale the number density to get the right opd. the opd at 0.8um is.. ~11.2
'''
#i calculated the number density that got me the opd i wanted 
opd_1,w0_1,g0_1,wavenumber_grid_1=vdi.calc_optics_user_r_dist(wave_in_1, ndz_1 ,r13, u.nm, n13/100, qext_1, qscat_1, cos_qscat_1, )
opd_4,w0_4,g0_4,wavenumber_grid_4=vdi.calc_optics_user_r_dist(wave_in_4, ndz_4 ,r4, u.nm, n4/100, qext_4, qscat_4, cos_qscat_4, )
opd_2,w0_2,g0_2,wavenumber_grid_2=vdi.calc_optics_user_r_dist(wave_in_2, ndz_2 ,r2, u.nm, n2/100, qext_2, qscat_2, cos_qscat_2, )
opd_31,w0_31,g0_31,wavenumber_grid_3=vdi.calc_optics_user_r_dist(wave_in_3, ndz_31 ,r13, u.nm, n13/100, qext_3, qscat_3, cos_qscat_3, )
opd_32,w0_32,g0_32,wavenumber_grid_3=vdi.calc_optics_user_r_dist(wave_in_3, ndz_32 ,r13, u.nm, n13/100, qext_3, qscat_3, cos_qscat_3, )
opd_1mod,w0_1mod,g0_1mod,wavenumber_grid_1mod=vdi.calc_optics_user_r_dist(wave_in_1mod, ndz_1mod ,r13, u.nm, n13/100, qext_1mod, qscat_1mod, cos_qscat_1mod, )

#opd_1mod = np.concatenate((np.zeros(15, dtype=opd_1mod.dtype), opd_1mod, np.zeros(3, dtype=opd_1mod.dtype)))
#w0_1mod = np.concatenate((np.zeros(15, dtype=opd_1mod.dtype), w0_1mod, np.zeros(3, dtype=opd_1mod.dtype)))
#g0_1mod = np.concatenate((np.zeros(15, dtype=opd_1mod.dtype), g0_1mod, np.zeros(3, dtype=opd_1mod.dtype)))

opd_1mod = np.concatenate((opd_1[0:15], opd_1mod, opd_1[-3:]))
g0_1mod = np.concatenate((g0_1[0:15], g0_1mod, g0_1[-3:]))
w0_1mod = np.concatenate((w0_1[0:15], w0_1mod, w0_1[-3:]))

#print("1 mod: "+str(opd_1mod[-7]))
#print("1: "+str(opd_1[-7]))
#print("2: "+str(opd_2[-7]))
#print("3: "+str(opd_31[-7]))
#print("4: "+str(opd_4[-7]))


wavenumber_grid_1mod = copy.deepcopy(wavenumber_grid_1)
nlayer = len(comp_file['pressure'])-1   #one - pt pressure 

#just set an arbitrary pressure grid from 1000 bars to a nanobar
pressure =  np.logspace(-9,3,nlayer) 

#****GOING TO 2.5 MICRONS****

#for each slab (5 total because we're splitting cloud 3 into two (cloud 4 is sandwiched in the middle)):
df_haze_1 = vdi.picaso_format_slab(base_pressure_1,opd_1[:-1], w0_1[:-1], g0_1[:-1], wavenumber_grid_1[:-1], pressure,p_top=haze_top_pressure_1)
df_haze_2 = vdi.picaso_format_slab(base_pressure_2,opd_2[:-1], w0_2[:-1], g0_2[:-1], wavenumber_grid_2[:-1], pressure,p_top=haze_top_pressure_2)
df_haze_31 = vdi.picaso_format_slab(base_pressure_31,opd_32[:-1], w0_32[:-1], g0_32[:-1], wavenumber_grid_3[:-1], pressure,p_top=haze_top_pressure_31)
df_haze_4 = vdi.picaso_format_slab(base_pressure_4,opd_4[:-1], w0_4[:-1], g0_4[:-1], wavenumber_grid_4[:-1], pressure,p_top=haze_top_pressure_4)
df_haze_32 = vdi.picaso_format_slab(base_pressure_32,opd_31[:-1], w0_31[:-1], g0_31[:-1], wavenumber_grid_3[:-1], pressure,p_top=haze_top_pressure_32)
#df_haze_1mod = vdi.picaso_format_slab(base_pressure_1,opd_1mod[:-1], w0_1mod[:-1], g0_1mod[:-1], wavenumber_grid_1mod[:-1], pressure,p_top=haze_top_pressure_1)
nwno = len(wavenumber_grid_1[:-1]) #this will depend on the refractive indices you have. Khare's tholin have wide coverage but low resolution

#**************CREATING ATMOSPHERE****************
opacity = pdi.opannection(wave_range=[0.3,2.5])
neptune = pdi.inputs()
#phase angle
neptune.phase_angle(0) #radians
#define gravity
neptune.gravity(radius=0.3463884070945,radius_unit=pdi.u.Unit('R_jup'), mass=0.053740779768177, mass_unit=pdi.u.Unit('M_jup'))
neptune.star(opacity, temp=5778, metal=0.00, logg=4.4374, semi_major = 30.178, semi_major_unit = u.Unit("au"), radius=1.0, radius_unit=pdi.u.R_sun,)#opacity db, pysynphot database, temp, metallicity, logg
neptune.atmosphere(filename="/Users/ihuckabee/Documents/picaso-master/neptune_stuff/ToShare/RefFiles/nepforpicaso_salr_abridged.txt", delim_whitespace=True)

#*********************this section is for plotting the different clouds
#CLOUD-FREE FOR REFERENCE
#wno, alb, fpfs = df['wavenumber'] , df['albedo'] , df['fpfs_reflected']

#SETTING UP CLOUD 3 AND TOTAL CLOUD
nlayer = len(comp_file['pressure'])-1 #one less than the number of PT points in your input
df_haze = copy.deepcopy(df_haze_31[["pressure", "wavenumber"]])
df_haze_3 = copy.deepcopy(df_haze_31[["pressure", "wavenumber"]])


df_haze_3['g0'] = df_haze_31['g0']+df_haze_32['g0']
df_haze_3['opd'] =  df_haze_31['opd']+df_haze_32['opd']
df_haze_3['w0'] = df_haze_31['w0']+df_haze_32['w0']
df_haze['g0'] = df_haze_3['g0']+df_haze_2['g0']+df_haze_1['g0']#+df_haze_4['g0']
df_haze['w0'] = df_haze_3['w0']+df_haze_2['w0']+df_haze_1['w0']#+df_haze_4['w0']
df_haze['opd'] = df_haze_3['opd']+df_haze_2['opd']+df_haze_1['opd']#+df_haze_4['opd']

df_1review = pd.read_csv('df_1_reviewer.csv')
df_2review = pd.read_csv('df_2_reviewer.csv')
df_31review = pd.read_csv('df_31_reviewer.csv')
df_32review = pd.read_csv('df_32_reviewer.csv')
df_4review = pd.read_csv('df_4_reviewer.csv')

df_darkreview = copy.deepcopy(df_2review[["pressure", "wavenumber"]])
df_darkreview['opd'] = df_2review['opd'] + df_31review['opd'] + df_32review['opd']
df_darkreview['g0'] = df_2review['g0'] + df_31review['g0'] + df_32review['g0']
df_darkreview['w0'] = df_2review['w0'] + df_31review['w0'] + df_32review['w0']



df_haze_1mod = pd.read_csv('df_1mod.csv')
df_darkmod = copy.deepcopy(df_haze_1mod[["pressure", "wavenumber"]])
#import pdb; pdb.set_trace()
df_darkmod['g0'] = df_haze_1mod['g0']#+df_2review['g0']#+ df_31review['g0']# + df_32review['g0']#+df_haze_3['g0']+df_haze_2['g0']
df_darkmod['w0'] = df_haze_1mod['w0']#+df_2review['w0']#+ df_31review['w0']# + df_32review['w0']#+df_haze_3['w0']+df_haze_2['w0']
df_darkmod['opd'] = df_haze_1mod['opd']#+ df_2review['opd']#+ df_31review['opd']# + df_32review['opd']#+df_haze_3['opd']+df_haze_2['opd']

neptune_1 = copy.deepcopy(neptune)
neptune_2 = copy.deepcopy(neptune)
neptune_3 = copy.deepcopy(neptune)
neptune_4 = copy.deepcopy(neptune)
neptune_tot = copy.deepcopy(neptune)
neptune_1mod = copy.deepcopy(neptune)
neptune_real = copy.deepcopy(neptune)
neptune_darkreview = copy.deepcopy(neptune)

neptune_darkreview.clouds(df=df_darkreview.astype(float))
df_um = neptune_tot.spectrum(opacity)
wno_um, alb_um = df_um['wavenumber'] , df_um['albedo']

neptune_tot.clouds(df=df_haze.astype(float))
nlayer = len(comp_file['pressure'])-1 #one less than the number of PT points in your input
df_tot = neptune_tot.spectrum(opacity)
wno_tot, alb_tot, fpfs_tot = df_tot['wavenumber'] , df_tot['albedo'] , df_tot['fpfs_reflected']
#wno_tot, alb_tot = pdi.mean_regrid(wno_tot, alb_tot, R=100)

neptune_3.clouds(df=df_haze_3.astype(float))
nlayer = len(comp_file['pressure'])-1 #one less than the number of PT points in your input
df_3 = neptune_3.spectrum(opacity)
wno_3, alb_3, fpfs_3 = df_3['wavenumber'] , df_3['albedo'] , df_3['fpfs_reflected']
wno_3, alb_3 = pdi.mean_regrid(wno_3, alb_3, R=100)

neptune_1.clouds(df=df_haze_1.astype(float))
nlayer = len(comp_file['pressure'])-1 #one less than the number of PT points in your input
df_1 = neptune_1.spectrum(opacity)
wno_1, alb_1, fpfs_1 = df_1['wavenumber'] , df_1['albedo'] , df_1['fpfs_reflected']
wno_1, alb_1 = pdi.mean_regrid(wno_1, alb_1, R=100)

neptune_2.clouds(df=df_haze_2.astype(float))
nlayer = len(comp_file['pressure'])-1 #one less than the number of PT points in your input
df_2 = neptune_2.spectrum(opacity)
wno_2, alb_2, fpfs_2 = df_2['wavenumber'] , df_2['albedo'] , df_2['fpfs_reflected']
wno_2, alb_2 = pdi.mean_regrid(wno_2, alb_2, R=100)

neptune_4.clouds(df=df_haze_4.astype(float))
nlayer = len(comp_file['pressure'])-1 #one less than the number of PT points in your input
df_4 = neptune_4.spectrum(opacity)
wno_4, alb_4, fpfs_4 = df_4['wavenumber'] , df_4['albedo'] , df_4['fpfs_reflected']
wno_4, alb_4 = pdi.mean_regrid(wno_4, alb_4, R=100)

neptune_1mod.clouds(df=df_darkmod.astype(float))
df_1mod = neptune_1mod.spectrum(opacity)
wno_1mod, alb_1mod, fpfs_1mod = df_1mod['wavenumber'] , df_1mod['albedo'] , df_1mod['fpfs_reflected']



import matplotlib.pyplot as plt
f = open("./ToShare/RefFiles/obsrefspec.txt",'r')
lines = f.readlines()[110:] #data doesn't start til 0.3 microns
wv_obs = []; alb_obs = []
for i in range(len(lines)):
    wv_obs.append(float(lines[i].split()[0]))
    alb_obs.append(float(lines[i].split()[1]))


#REAL CLOUDS THAT WORK 

df_1 = pd.read_csv('/Users/ihuckabee/Documents/picaso-master/neptune_stuff/ToShare/RefFiles/cloud1_v3.csv')
df_2 = pd.read_csv('/Users/ihuckabee/Documents/picaso-master/neptune_stuff/ToShare/RefFiles/cloud2_v3.csv')
df_3 = pd.read_csv('/Users/ihuckabee/Documents/picaso-master/neptune_stuff/ToShare/RefFiles/cloud3_v3.csv')
df_4 = pd.read_csv('/Users/ihuckabee/Documents/picaso-master/neptune_stuff/ToShare/RefFiles/cloud4_v3.csv')


df_works =  copy.deepcopy(df_1[["pressure", "wavenumber"]])
df_works['opd'] = df_2['opd']+df_3['opd']#+df_1['opd']
df_works['w0'] = df_2['w0']+df_3['w0']#+df_1['w0']
df_works['g0'] = df_2['g0']+df_3['g0']#+df_1['g0']

neptune_real.clouds(df=df_works.astype(float))
nlayer = len(comp_file['pressure'])-1 #one less than the number of PT points in your input
df_1_works = neptune_real.spectrum(opacity)
wno_1_works, alb_1_works, fpfs_1_work = df_1_works['wavenumber'] , df_1_works['albedo'] , df_1_works['fpfs_reflected']
#wno_1_works, alb_1_works = pdi.mean_regrid(wno_1_works, alb_1_works, R=100)



f = open("nds-2018.txt",'r')
lines = f.readlines()[1:] 
wv_nds = []; alb_nds = []
for i in range(len(lines)):
    wv_nds.append(float(lines[i].split()[0]))
    alb_nds.append(float(lines[i].split()[1]))




#plt.plot(1e4/wno_1, alb_1, label = "One I Just Made")
#plt.plot(1e4/wno_2, alb_2, label = "2")
#plt.plot(1e4/wno_3, alb_3, label = "3")
#plt.plot(1e4/wno_4, alb_4, label = "4")
plt.plot(1e4/wno_1_works, alb_1_works, label = "The Original Dark Spot")
#plt.plot(1e4/wno_um, alb_um, label = 'Newly Generated Dark Spot (no Aerosol 1)')
#plt.plot(1e4/wno_tot, alb_tot, label = 'Cloud Base')
plt.plot(1e4/wno_1mod, alb_1mod, 'k--', label = "Newly Generated Dark Spot (with modified Aerosol 1)", alpha = 0.8)
plt.plot(wv_obs,alb_obs, 'k-.', alpha = 0.5, label = 'HST/SpeX data')
plt.plot(wv_nds,alb_nds,'r-', label="NDS Fit by Irwin+23", alpha = 0.8)
plt.xlim(0.45,0.95)
plt.ylabel('Albedo')
plt.xlabel('Wavelength (um)')
plt.legend()
plt.show()
import pdb; pdb.set_trace()

#save the clouds - this feeds into the post processing code 
#df_haze_1.to_csv("./ToShare/RefFiles/cloud1_trunc.csv")
#df_haze_2.to_csv("./ToShare/RefFiles/cloud2_trunc.csv")
#df_haze_3.to_csv("./ToShare/RefFiles/cloud3_trunc.csv")
#df_haze_4.to_csv("./ToShare/RefFiles/cloud4_trunc.csv")

#plot
#ppi.show(ppi.spectrum([wno_4,wno_tot], [alb_4, alb_tot],legend = ["Stacked cloud 4", "Full Cloud"], plot_width=750))

pdb.set_trace()