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
ulim=1.7

for i in range(100):
    data = np.loadtxt('data/scores/scores_%d.txt' % i ,skiprows=2)
    plt.figure()
    plt.imshow(data,vmin=-ulim,vmax=ulim,cmap='RdBu_r')
    plt.savefig('data/scores/scores_%04d.png' % i)
    plt.close()

