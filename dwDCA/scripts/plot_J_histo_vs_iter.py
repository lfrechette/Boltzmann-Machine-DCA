#Created by Layne Frechette, March 11th, 2021

import numpy as np
import pylab as plt
import sys

maxiter=300
ulim=2.5
nbins=30

for i in range(maxiter):
    data = np.loadtxt('data/couplings/J1_%d.txt' % i ,skiprows=2)
    plt.figure()
    hist,bins=np.histogram(data,bins=nbins,density=True)
    center=(bins[:-1]+bins[1:])/2
    plt.plot(center,hist)
    plt.xlim([-ulim,ulim])
    plt.ylim([-0.01,5.0])
    plt.savefig('data/couplings/J1_hist_%04d.png' % i)
    plt.close()

