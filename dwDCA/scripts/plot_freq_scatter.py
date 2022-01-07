#Created by Layne Frechette, March 23rd, 2021

import numpy as np
import pylab as plt
import sys
import matplotlib.cm as cm

mydir = sys.argv[1]

align = np.loadtxt(mydir + 'stat_align_1p.txt')
N = align.shape[0]
q=21
align = align[:,1:]
#align = align.flatten()
xvals = np.arange(N)+1

cmap = plt.get_cmap('jet')
colors = [cmap(i) for i in np.linspace(0,1,q)]

aa = ['-','A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']
ulim=0.5
maxiter=99

for i in range(maxiter):
    MC = np.loadtxt(mydir + ('stat_MC_1p_%d.txt' % (i+1)))
    MC = MC[:,1:]#.flatten()
    plt.figure()
    for j in range(q):
        plt.scatter(xvals,align[:,j]-MC[:,j],s=1,color=colors[j],label='%s'%aa[j])
    plt.legend(loc='upper center', ncol=6,fontsize='xx-small',columnspacing=0.5,bbox_to_anchor=(0.5,1.5))
    plt.ylim([-ulim,ulim])
    plt.xlabel('i')
    plt.ylabel(r'$f_i^{\text{MSA}}-f_i^{\text{MC}}$')
    plt.savefig('data/ranganathan/images/freq_diff_scatter_%04d.png' % i, bbox_inches='tight')
    plt.close()

