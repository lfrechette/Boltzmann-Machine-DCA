#Created by Layne Frechette, March 23rd, 2021

import numpy as np
import pylab as plt
import sys
import matplotlib.cm as cm
from matplotlib import colors

mydir = sys.argv[1]

align = np.loadtxt(mydir + 'stat_align_2p.txt')
N = align.shape[0]
q=21
align = align[:,2:]
xvals = np.arange(N)+1

lower = 0.1
upper = 0.9
mymap = np.zeros(align.shape)
mymap[align<lower]=1
mymap[(align>=lower) & (align<upper)]=2
mymap[align>=upper]=3

cmap = colors.ListedColormap(['b','k','r'])

ulim=3
maxiter=99

for i in range(maxiter):
    h = np.loadtxt(mydir + ('parameters_J_%d.txt' % (i+1)))
    h = h[:,2:]#.flatten()
    plt.figure()
    for j in range(q**2):
        plt.scatter(xvals,h[:,j],s=1,c=mymap[:,j],cmap=cmap)
    plt.ylim([-ulim,ulim])
    plt.xlabel(r'pair index')
    plt.ylabel(r'$J_{ij}(\alpha,\beta)$')
    plt.title(r'iter %d' % (i+1))
    plt.savefig(mydir+'images/freq_J_scatter_%04d.png' % i, bbox_inches='tight')
    plt.close()

