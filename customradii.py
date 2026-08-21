import numpy as np
import matplotlib.pyplot as plt
import pdb
from scipy.stats import describe

'''
The paper i followed gave me the radius distrituion (mean and var) for the particles in each cloud
Technically I have four flavors of clouds/hazes
but clouds 1 and 3 have the same radius distribution (all gamma)
Need this info to feed into virga 
'''
np.random.seed(42)
shape_13, scale_13 = 0.05, 1.  
#mean = shape*scale, var = shape*scale^2
dist = np.array(1000*(np.random.gamma(shape_13,scale_13,10000)),dtype = int)#1000*np.random.gamma(shape_13,scale_13,10000)
dist_13 =np.array(dist, dtype = int) # 72.5K and 165K
#^^ temp based on pt profile and where the clouds are in pressure space
import pdb; pdb.set_trace()
shape_2, scale_2 = 0.68**2/0.25, 0.25/0.68 #for gamma
#shape_2, scale_2 = 0.68, 0.25
dist_2 = np.array(1000*(np.random.gamma(shape_2,scale_2,10000)),dtype = int)#70K
#i want 0.68 = shape*scale, 0.25 = shape*scale^2

shape_4, scale_4 = 2.8**2/0.5, 0.5/2.8
dist_4 = np.array(1000*(np.random.gamma(shape_4,scale_4,10000)),dtype = int)#53.8K



bins_4 = np.arange(1,np.mean(dist_4)+3*(np.std(dist_4)),1)
bins_low = np.array([0.,0.0001,0.001,0.01,0.1]);bins_high = np.arange(1,np.mean(dist_13)+3*np.std(dist_13),1)
bins = np.concatenate((bins_low,bins_4))
counts_1, res_bins, patches = plt.hist(dist_13,bins = bins); counts_1 = np.array(counts_1/10000.) #146K
counts_2, res_bins, patches = plt.hist(dist_2,bins = bins); counts_2 = np.array(counts_2/10000.) #84K
counts_3_1, res_bins, patches = plt.hist(dist_13,bins = bins); counts_3_1 = np.array(counts_3_1/10000.)#79.8K
counts_4, res_bins, patches = plt.hist(dist_4,bins = bins); counts_4 = np.array(counts_4/10000.) #54.7K
counts_3_2, res_bins, patches = plt.hist(dist_13,bins = bins); counts_3_2 = np.array(counts_3_2/10000.)#72.5K

#plt.hist(dist_2,density=True)
#plt.hist(dist_4,density=True)
#plt.show()
#pdb.set_trace()

bins = bins[1:]
file1 = open("nepcustomradii.txt", "w")
file1.write("size_nm 146 84 79.8 54.7 72.5 \n")
for i in range(len(counts_1)):
    file1.write(str(bins[i])+" "+str(counts_1[i])+" "+str(counts_2[i])+" "+str(counts_3_1[i])+" "+str(counts_4[i])+" "+str(counts_3_2[i])+" \n")
file1.close()


f = open("nepcustomradii.txt", "r")
lines = f.readlines()
temp_labels = lines[0].split()
n_146 = []; n_84 = []; n_79 = []; n_54 = []; n_73 = []
size = []

for i in lines[1:]:
    size.append(float(i.split()[0]))
    n_146.append(float(i.split()[1]))
    n_84.append(float(i.split()[2]))
    n_79.append(float(i.split()[3]))
    n_54.append(float(i.split()[4]))
    n_73.append(float(i.split()[5]))
plt.close()
#plt.plot(size,n_146, label = temp_labels[1])
plt.plot(size, n_84, label = "Aerosol 2")
plt.plot(size, n_79, label = "Aerosol 1 and 3")
plt.plot(size, n_54,label = "Aerosol 4")
#plt.plot(size, n_73, label = temp_labels[5])
plt.ylabel("distribution percent total")
plt.xlabel("particle radius (nm)")
plt.ylim(0,0.004)
plt.legend()
plt.show()
pdb.set_trace()

