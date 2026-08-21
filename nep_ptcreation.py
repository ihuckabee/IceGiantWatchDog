#import os
#os.environ["picaso_refdata"]="/Users/ihuckabee/Documents/picaso-master/reference"
#os.environ['PYSYN_CDBS']="/Users/ihuckabee/Documents/picaso-master/grp/redcat/trds"
import numpy as np
import os
import pdb
#refdata = os.getenv("picaso_refdata")
#from picaso import justdoit as jdi
#from picaso import justplotit as jpi
#pt file creation

file = "/Users/ihuckabee/Downloads/neptunemodel.txt"
f=open(file,"r")
lines=f.readlines()[4:]
result=[]
for x in lines:
    result.append(x.split())
pressure = []; temp_dalr = []; temp_salr = []; ch4 = []; h2s = []
molwt = []; cpr = []; brunt1 = []; brunt2 = []; eddy1 = []; eddy2 = []
count = 1
for x in result:
    if count%2 != 0:
        pressure.append(float(x[0])); temp_dalr.append(float(x[1])); temp_salr.append(float(x[2])); ch4.append(float(x[3])); h2s.append(float(x[4])); molwt.append(float(x[5]))
    else:
        cpr.append(float(x[0])); brunt1.append(float(x[1])); brunt2.append(float(x[2])); eddy1.append(float(x[3])); eddy2.append(float(x[4]))
    count+=1
f.close()
chem_other = 1 - np.array(ch4)+np.array(h2s)
h2 = 0.81*chem_other
he = 0.19*chem_other
pressure_abridged = []; temp_salr_abridged = []; ch4_abridged = []; h2s_abridged = []; h2_abridged = []; he_abridged = [] 
for i in range(len(pressure)):
    if i%15 == 0:
        pressure_abridged.append(pressure[i])
        temp_salr_abridged.append(temp_salr[i])
        ch4_abridged.append(ch4[i])
        h2s_abridged.append(h2s[i])
        h2_abridged.append(h2[i])
        he_abridged.append(he[i])

#file1 = open("nepforpicaso_salr_abridged.txt", "w")
#file1.write("pressure temperature CH4 H2S H2 He\n")
#pdb.set_trace()
pressure_abridged = pressure_abridged[::-1]
temp_salr_abridged = temp_salr_abridged[::-1]
ch4_abridged = ch4_abridged[::-1]
h2s_abridged = h2s_abridged[::-1]
h2_abridged = h2_abridged[::-1]
he_abridged = he_abridged[::-1]
#pdb.set_trace()
#for i in range(len(pressure_abridged)):
    #file1.write(str(pressure_abridged[i])+" "+str(temp_salr_abridged[i])+" "+str(ch4_abridged[i])+" "+str(h2s_abridged[i])+" "+str(h2_abridged[i])+" "+str(he_abridged[i])+" \n")
#file1.close()
#pdb.set_trace()
import matplotlib.pyplot as plt
plt.figure(figsize=(4,6))

plt.plot(temp_salr_abridged,pressure_abridged, linewidth=3)
plt.ylim(10**-2,10**2)
plt.xlim(0,400)
# Log scale for pressure axis
plt.yscale('log')
# Invert y-axis so larger pressures are at the bottom
plt.gca().invert_yaxis()

#plt.show()
plt.savefig('PTPlot_Neptune.pdf')


