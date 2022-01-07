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

ulim=0.5
maxiter=999

for i in range(maxiter):
    MC = np.loadtxt(mydir + ('stat_MC_1p_%d.txt' % (i+1)))
    MC = MC[:,1:]#.flatten()
    plt.figure()
    for j in range(q):
        plt.scatter(xvals,align[:,j]-MC[:,j],s=1,c=mymap[:,j],cmap=cmap)
    plt.ylim([-ulim,ulim])
    plt.xlabel('i')
    plt.ylabel(r'$f_i^{\text{MSA}}-f_i^{\text{MC}}$')
    plt.savefig(mydir+'images/freq_1p_diff_scatter_%04d.png' % i, bbox_inches='tight')
    plt.close()

