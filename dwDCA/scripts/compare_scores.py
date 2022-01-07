#Created by Layne Frechette, March 11th, 2021

import numpy as np
import pylab as plt
import sys
plt.rcParams.update({'font.size': 8})

family = sys.argv[1]

none_eps01_mc1000 = np.loadtxt('data/scores_%s_init=none_maxiter=1000_stepsize=0.100000_mcsteps=1000.txt' % (family))
none_eps01_mc10000 = np.loadtxt('data/scores_%s_init=none_maxiter=1000_stepsize=0.100000_mcsteps=10000.txt' % (family))
none_eps01_mc100000 = np.loadtxt('data/scores_%s_init=none_maxiter=1000_stepsize=0.100000_mcsteps=100000.txt' % (family))
none_eps001_mc1000 = np.loadtxt('data/scores_%s_init=none_maxiter=1000_stepsize=0.010000_mcsteps=1000.txt' % (family))
none_eps001_mc10000 = np.loadtxt('data/scores_%s_init=none_maxiter=1000_stepsize=0.010000_mcsteps=10000.txt' % (family))
none_eps001_mc100000 = np.loadtxt('data/scores_%s_init=none_maxiter=1000_stepsize=0.010000_mcsteps=100000.txt' % (family))

bm_eps01_mc1000 = np.loadtxt('data/scores_%s_init=bm_params.txt_maxiter=1000_stepsize=0.100000_mcsteps=1000.txt' % (family))
bm_eps01_mc10000 = np.loadtxt('data/scores_%s_init=bm_params.txt_maxiter=1000_stepsize=0.100000_mcsteps=10000.txt' % (family))
bm_eps01_mc100000 = np.loadtxt('data/scores_%s_init=bm_params.txt_maxiter=1000_stepsize=0.100000_mcsteps=100000.txt' % (family))
bm_eps001_mc1000 = np.loadtxt('data/scores_%s_init=bm_params.txt_maxiter=1000_stepsize=0.010000_mcsteps=1000.txt' % (family))
bm_eps001_mc10000 = np.loadtxt('data/scores_%s_init=bm_params.txt_maxiter=1000_stepsize=0.010000_mcsteps=10000.txt' % (family))
bm_eps001_mc100000 = np.loadtxt('data/scores_%s_init=bm_params.txt_maxiter=1000_stepsize=0.010000_mcsteps=100000.txt' % (family))

#None init
fig,ax=plt.subplots(2,3,figsize=(8,5))

eps = [0.1,0.01]
steps = [1000,10000,100000]
nonelist = [[none_eps01_mc1000,none_eps01_mc10000,none_eps01_mc100000],[none_eps001_mc1000,none_eps001_mc10000,none_eps001_mc100000]]
maxlist = [[np.max(none_eps01_mc1000),np.max(none_eps01_mc10000),np.max(none_eps01_mc100000)],[np.max(none_eps001_mc1000),np.max(none_eps001_mc10000),np.max(none_eps001_mc100000)]]

for i in range(2):
    for j in range(3):
        data = nonelist[i][j]
        data = data-np.min(data)
        data = data/np.max(data)
        im=ax[i,j].imshow(data)
        ax[i,j].set_title(r'$%d$ MC sweeps, $\epsilon=%.02f$, $J_{\text{max}}=%.02f$' % (steps[j], eps[i], maxlist[i][j]),fontsize=8)
        ax[i,j].tick_params(axis='both', which='major', labelsize=8)
        if(i==1):
            ax[i,j].set_xlabel(r'site $i$',fontsize=8)
        if(j==0):
            ax[i,j].set_ylabel(r'site $j$', fontsize=8)

cbar_ax = fig.add_axes([0.925,0.15,0.025,0.7])
fig.colorbar(im,cax=cbar_ax)
plt.savefig('scores_init=none.png', dpi=300, bbox_inches='tight')


#bm init
fig,ax=plt.subplots(2,3,figsize=(8,5))

eps = [0.1,0.01]
steps = [1000,10000,100000]
bmlist = [[bm_eps01_mc1000,bm_eps01_mc10000,bm_eps01_mc100000],[bm_eps001_mc1000,bm_eps001_mc10000,bm_eps001_mc100000]]
maxlist = [[np.max(bm_eps01_mc1000),np.max(bm_eps01_mc10000),np.max(bm_eps01_mc100000)],[np.max(bm_eps001_mc1000),np.max(bm_eps001_mc10000),np.max(bm_eps001_mc100000)]]

for i in range(2):
    for j in range(3):
        data = bmlist[i][j]
        data = data-np.min(data)
        data = data/np.max(data)
        im=ax[i,j].imshow(data)
        ax[i,j].set_title(r'$%d$ MC sweeps, $\epsilon=%.02f$, $J_{\text{max}}=%.02f$' % (steps[j], eps[i], maxlist[i][j]),fontsize=8)
        ax[i,j].tick_params(axis='both', which='major', labelsize=8)
        if(i==1):
            ax[i,j].set_xlabel(r'site $i$',fontsize=8)
        if(j==0):
            ax[i,j].set_ylabel(r'site $j$', fontsize=8)

cbar_ax = fig.add_axes([0.925,0.15,0.025,0.7])
fig.colorbar(im,cax=cbar_ax)
plt.savefig('scores_init=bm.png', dpi=300, bbox_inches='tight')

plt.show()
