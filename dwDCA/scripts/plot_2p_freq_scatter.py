#Created by Layne Frechette, March 23rd, 2021

import numpy as np
import pylab as plt
import sys
import matplotlib.cm as cm
from matplotlib import colors

mydir = sys.argv[1]

align = np.loadtxt(mydir + 'stat_align_2p.txt')
q=21
Npair = align.shape[0]
align = align[:,2:]
#align=align.reshape((Npair,q*q))
#align = align.flatten()
xvals = np.arange(Npair)+1

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
    MC = np.loadtxt(mydir + ('stat_MC_2p_%d.txt' % (i+1)))
    MC = MC[:,2:]#.flatten()
    #MC = MC.reshape((Npair,q**2))
    plt.figure()
    for j in range(q**2):
        plt.scatter(xvals,align[:,j]-MC[:,j],s=1,c=mymap[:,j],cmap=cmap)
    #plt.legend(loc='upper center', ncol=6,fontsize='xx-small',columnspacing=0.5,bbox_to_anchor=(0.5,1.5))
    plt.ylim([-ulim,ulim])
    plt.xlabel('pair index')
    plt.ylabel(r'$f_{ij}^{\text{MSA}}-f_{ij}^{\text{MC}}$')
    plt.title('iter %d' % (i+1))
    plt.savefig(mydir+'images/freq_2p_diff_scatter_%04d.png' % i, bbox_inches='tight')
    plt.close()

