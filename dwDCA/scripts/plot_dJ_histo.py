#Created by Layne Frechette, March 11th, 2021

import numpy as np
import pylab as plt
import sys

maxiter=100
ulim=0.25
nbins=30

for i in range(maxiter):
    data = np.loadtxt('data/derivs/J/dJ1_%d.txt' % i ,skiprows=2)
    plt.figure()
    hist,bins=np.histogram(data,bins=nbins,density=True)
    center=(bins[:-1]+bins[1:])/2
    plt.plot(center,hist)
    plt.xlim([-ulim,ulim])
    plt.ylim([-0.01,120.0])
    plt.title('iter=%d' % i)
    plt.xlabel(r'$dJ$')
    plt.ylabel(r'$P(dJ)$')
    plt.savefig('data/derivs/J/dJ1_hist_%04d.png' % i, bbox_inches='tight')
    plt.close()

