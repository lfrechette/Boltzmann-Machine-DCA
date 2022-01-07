#Created by Layne Frechette, March 23rd, 2021

import numpy as np
import pylab as plt
import sys
import matplotlib.cm as cm
from matplotlib import colors

mydir = sys.argv[1]

align = np.loadtxt(mydir + 'stat_align_1p.txt')
N = align.shape[0]
q=21
align = align[:,1:]
xvals = np.arange(N)+1

lower = 0.1
upper = 0.9
mymap = np.zeros(align.shape)
mymap[align<lower]=1
mymap[(align>=lower) & (align<upper)]=2
mymap[align>=upper]=3

cmap = colors.ListedColormap(['b','k','r'])

ulim=10
maxiter=999

for i in range(maxiter):
    h = np.loadtxt(mydir + ('parameters_h_%d.txt' % (i+1)))
    h = h[:,1:]#.flatten()
    plt.figure()
    for j in range(q):
        plt.scatter(xvals,h[:,j],s=1,c=mymap[:,j],cmap=cmap)
    plt.ylim([-ulim,ulim])
    plt.xlabel(r'$i$')
    plt.ylabel(r'$h_i(\alpha)$')
    plt.savefig(mydir+'images/freq_h_scatter_%04d.png' % i, bbox_inches='tight')
    plt.close()

