#Created by Layne Frechette, March 11th, 2021

import numpy as np
import pylab as plt
import sys

maxiter=1000
ulim=10.0
nbins=30

for i in range(maxiter):
    data = np.loadtxt('data/fields/h1_%d.txt' % i ,skiprows=2)
    plt.figure()
    hist,bins=np.histogram(data,bins=nbins,density=True)
    center=(bins[:-1]+bins[1:])/2
    plt.plot(center,hist)
    plt.xlim([-ulim,ulim])
    plt.ylim([-0.01,1.0])
    plt.savefig('data/fields/h1_hist_%04d.png' % i)
    plt.close()

