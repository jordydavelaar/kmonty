

import numpy as np
import matplotlib.pyplot as plt
from pylab import *
from scipy import signal
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.gridspec as gridspec

start=0
Nimg=1
num_freq=1
Np=1
Nt=90
data=np.zeros(shape=(2,200))
for j in range(start,start+Nimg,1):
	data_i = np.loadtxt("spec_sphere_%d.dat"%(j),unpack=True)
	print("file=",j)
	data[0]=data_i[0]
	for k in range(0, Np):
		for l in range(0,Nt):
			index= int( 180 + k*2 + 360*l)
			#print(index,k,l,Np)
			data[1]+=(data_i[index])/(1.*Np*Nt*Nimg)
data=np.transpose(data)
np.savetxt("spec_total_sphere_pwl.dat",data)

