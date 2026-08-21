import numpy as np
import matplotlib.pyplot as plt
import pickle
import pdb
sm = 4.514565e12 #m
rp = 2.47639e7 #m
fnames = ['./RefFiles/10x10_cldtest.pic', './RefFiles/10x10cloudfree.pic']
wno_3d = []; alb_3d = []
for i in range(len(fnames)):
        infile = open(fnames[i], 'rb')
        out3d = pickle.load(infile)
        #pdb.set_trace()
        wno,fpfs = out3d[0],out3d[1]
        alb = fpfs*(sm/rp)**2
        wno_3d.append(1e4/wno)
        alb_3d.append(alb)

j_band = [1.1,1.4] #micrometers 
band_index = np.where(((wno_3d[0]<j_band[1]) & (wno_3d[0]>j_band[0])))[0]
band_sum = np.sum(alb_3d[0][band_index[0]:band_index[-1]])
band_width = wno_3d[0][band_index[0]] - wno_3d[0][band_index[-1]]
del_lam = wno_3d[0][band_index[0]] - wno_3d[0][band_index[1]]
phot = band_sum*del_lam / band_width 
pdb.set_trace()