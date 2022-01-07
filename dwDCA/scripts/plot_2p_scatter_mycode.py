#Created by Layne Frechette, March 23rd, 2021

import numpy as np
import pylab as plt
import sys
import matplotlib.cm as cm
from matplotlib import colors

mydir = sys.argv[1]

align = np.loadtxt(mydir + 'stat_align_2p.txt', skiprows=2)
q=2
N = int(align.shape[0]/q)
align=align.reshape((N,q*q))
#align = align[:,1:]
xvals = np.arange(N)+1

lower = 0.1
upper = 0.9
mymap = np.zeros(align.shape)
mymap[align<lower]=1
mymap[(align>=lower) & (align<upper)]=2
mymap[align>=upper]=3

cmap = colors.ListedColormap(['b','k','r'])

ulim=0.5
maxiter=99

for i in range(maxiter):
    mc = np.loadtxt(mydir + ('stat_MC_2p_%d.txt' % (i+1)), skiprows=2)
    mc=mc.reshape((N,q**2))
    #h = h[:,1:]#.flatten()
    plt.figure()
    for j in range(q**2):
        plt.scatter(xvals,align[:,j]-mc[:,j],s=10,c=mymap[:,j],cmap=cmap)
    plt.ylim([-ulim,ulim])
    plt.xlabel(r'pair index')
    plt.ylabel(r'$f_{ij}^{MSA}-f_{ij}^{MC}$')
    plt.title('iter %d' % (i+1))
    plt.savefig(mydir+'images/freq_2p_scatter_%04d.png' % i, bbox_inches='tight')
    plt.close()

