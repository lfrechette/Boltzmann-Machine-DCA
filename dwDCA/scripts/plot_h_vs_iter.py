#Created by Layne Frechette, March 11th, 2021

import numpy as np
import pylab as plt
import sys

'''
family = sys.argv[1]
init = sys.argv[2]
maxiter = int(sys.argv[3])
eps = float(sys.argv[4])
mcsteps = int(sys.argv[5])
lam = float(sys.argv[6])
'''

maxiter=1000
ulim=1.0

for i in range(maxiter):
    data = np.loadtxt('data/fields/h1_%d.txt' % i ,skiprows=2)
    plt.figure()
    plt.imshow(data,vmin=-ulim,vmax=ulim,cmap='RdBu_r')
    plt.savefig('data/fields/h1_smallrange_%04d.png' % i)
    plt.close()

