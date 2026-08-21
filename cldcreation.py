import numpy as np
import pdb
from scipy.interpolate import interp1d
file = "/Users/ihuckabee/Downloads/neptune_cloud_contour_exp.txt"
f=open(file,"r")
pressurelines=f.readlines()[13:20]
f=open(file,"r")
oplines = f.readlines()[22:]
f=open(file,"r")
wvlines=f.readlines()[3:11]
pressure = []
for i in pressurelines:
    pressure.append(np.array(i.split()).astype(float))
pressure = np.concatenate(pressure)
opacity = []
for i in oplines:
    opacity.append(np.array(i.split()).astype(float))
opacity = np.concatenate(opacity)
wv = []
for i in wvlines:
    wv.append(np.array(i.split()).astype(float))
wv = np.concatenate(wv)[:-1]
opacity_matrix = []
new_wv = []
for i in range(len(opacity)):
    new_wv.append(opacity[i])
    if len(new_wv) == len(wv):
        opacity_matrix.append(new_wv)
        new_wv = []
loc_aero = []
aerosol1 = np.round((np.loadtxt("/Users/ihuckabee/Downloads/aerosol_1_properties.txt")).T,2)
aerosol2 = (np.loadtxt("/Users/ihuckabee/Downloads/aerosol_2_properties.txt")).T
aerosol3 = (np.loadtxt("/Users/ihuckabee/Downloads/aerosol_3_properties.txt")).T
aerosol4 = (np.loadtxt("/Users/ihuckabee/Downloads/aerosol_4_properties.txt")).T
#pdb.set_trace()
for i in range(len(wv)):
    if np.where(aerosol1[0] == wv[i])[0]!= []:
        #print("1")
        #pdb.set_trace()
        locs = np.where(aerosol1[0] == wv[i])[0]
        mid = int(np.round(len(locs)/2,0))
        loc_aero.append(locs[mid])

#pdb.set_trace()
#file1 = open("nepcloud_irwin21.txt", "w")
#file1.write("lvl wv opd g0 w0 \n")
count = 0
opd = []
g0 = []
w0 = []

pressure_mark = [183,192,201,210,220,229,238,247,256,265,274,283,293,302,311,320,329,338,348,357,366,375,384,393,403,412,421,430,439,448,457,466,476,485,494,503,512,521,530]

for i in range(998):
    for j in range(len(wv)):
        try:
            if i >=pressure_mark[count+1]:
                count+=1
            if i >= pressure_mark[count] and i <= pressure_mark[count+1]: 
                opd.append(opacity_matrix[count][j])
                if count < 6: 
                    g0.append(aerosol1[1][loc_aero[j]])
                    w0.append(aerosol1[2][loc_aero[j]])
                if count >= 6 and count < 9:
                    g0.append(aerosol2[1][loc_aero[j]])
                    w0.append(aerosol2[2][loc_aero[j]])
                if count >= 9 and count < 15:
                    g0.append(aerosol3[1][loc_aero[j]])
                    w0.append(aerosol3[2][loc_aero[j]])
                if count >= 15 and count < 20:
                    g0.append(aerosol4[1][loc_aero[j]])
                    w0.append(aerosol4[2][loc_aero[j]])
                if count >= 20:
                    g0.append(aerosol3[1][loc_aero[j]])
                    w0.append(aerosol3[2][loc_aero[j]])
            else:
                opd.append(0.)
                g0.append(0.)
                w0.append(0.)
        except:
            opd.append(0.)
            g0.append(0.)
            w0.append(0.)

file_pressure = "/Users/ihuckabee/Documents/picaso-master/nepforpicaso_salr.txt"
f=open(file_pressure,"r")
lines=f.readlines()[1:]
pressure_picaso = []
for i in lines:
    pressure_picaso.append(np.array(i.split()).astype(float))
pressure_picaso = np.array(pressure_picaso).T[0]
file1 = open("nepcloud_irwin21.cld", "w")
file1.write("nlayer nwave pressure wavenumber opd g0 w0  \n")
#file1.write("nlayer nwave wavenumber w0 g0 opd \n")
nlvl = 1
nwv = 1
nnum = 0
wavenum = np.round(1/(wv/1e4),3) #cm
#pdb.set_trace()
for i in range(len(opd)):
    if nwv > len(wv):
        #pdb.set_trace()
        nnum = 0
        nwv = 1
        nlvl+=1
    #pdb.set_trace()
    file1.write(str(nlvl)+" "+str(nwv)+" "+str(pressure_picaso[nlvl])+" "+str(wavenum[nnum])+" "+str(opd[i])+" "+str(g0[i])+" "+str(w0[i])+" \n")
    nwv+=1
    nnum+=1
file1.close()
pdb.set_trace()