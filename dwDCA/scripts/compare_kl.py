#Created by Layne Frechette, March 11th, 2021

import numpy as np
import pylab as plt
import sys
plt.rcParams.update({'font.size': 8})

family = sys.argv[1]

none_eps01_mc1000 = np.loadtxt('data/kl_%s_none_1000_0.100000_1000_0.010000.txt' % (family))
none_eps01_mc10000 = np.loadtxt('data/kl_%s_none_1000_0.100000_10000_0.010000.txt' % (family))
none_eps01_mc100000 = np.loadtxt('data/kl_%s_none_1000_0.100000_100000_0.010000.txt' % (family))
none_eps001_mc1000 = np.loadtxt('data/kl_%s_none_1000_0.010000_1000_0.010000.txt' % (family))
none_eps001_mc10000 = np.loadtxt('data/kl_%s_none_1000_0.010000_10000_0.010000.txt' % (family))
none_eps001_mc100000 = np.loadtxt('data/kl_%s_none_1000_0.010000_100000_0.010000.txt' % (family))

bm_eps01_mc1000 = np.loadtxt('data/kl_%s_bm_params.txt_1000_0.100000_1000_0.010000.txt' % (family))
bm_eps01_mc10000 = np.loadtxt('data/kl_%s_bm_params.txt_1000_0.100000_10000_0.010000.txt' % (family))
bm_eps01_mc100000 = np.loadtxt('data/kl_%s_bm_params.txt_1000_0.100000_100000_0.010000.txt' % (family))
bm_eps001_mc1000 = np.loadtxt('data/kl_%s_bm_params.txt_1000_0.010000_1000_0.010000.txt' % (family))
bm_eps001_mc10000 = np.loadtxt('data/kl_%s_bm_params.txt_1000_0.010000_10000_0.010000.txt' % (family))
bm_eps001_mc100000 = np.loadtxt('data/kl_%s_bm_params.txt_1000_0.010000_100000_0.010000.txt' % (family))

###
'''
plt.figure()
plt.plot(data1[:,0],data1[:,1],label=r'$10^3 MC sweeps')
plt.plot(data2[:,0],data2[:,1],label=r'$10^4 MC sweeps')
plt.plot(data3[:,0],data3[:,1],label=r'$10^5 MC sweeps')
plt.savefig('data/freq_diff_1p_%s_init=%s_maxiter=%d_stepsize=%f_compare_mcsteps.png' % (family,init,maxiter,eps), dpi=300, bbox_inches='tight')

plt.figure()
plt.plot(data1[:,0],data1[:,2],label=r'$10^3 MC sweeps')
plt.plot(data2[:,0],data2[:,2],label=r'$10^4 MC sweeps')
plt.plot(data3[:,0],data3[:,2],label=r'$10^5 MC sweeps')
plt.savefig('data/freq_diff_2p_%s_init=%s_maxiter=%d_stepsize=%f_compare_mcsteps.png' % (family,init,maxiter,eps), dpi=300, bbox_inches='tight')
'''

#None init
fig,ax=plt.subplots(2,3,figsize=(8.5,4))

eps = [0.1,0.01]
steps = [1000,10000,100000]
nonelist = [[none_eps01_mc1000,none_eps01_mc10000,none_eps01_mc100000],[none_eps001_mc1000,none_eps001_mc10000,none_eps001_mc100000]]

for i in range(2):
    for j in range(3):
        data = nonelist[i][j]
        ax[i,j].plot(data[:,0],data[:,1],label=r'$%d$ MC sweeps, $\epsilon=%.02f$' % (steps[j], eps[i]))
        ax[i,j].legend(fontsize=6)
        ax[i,j].tick_params(axis='both', which='major', labelsize=8)
        if(i==1):
            ax[i,j].set_xlabel('iterations',fontsize=8)
        if(j==0):
            ax[i,j].set_ylabel(r'$D_1$', fontsize=8)
plt.savefig('kl_h_init=none.png', dpi=300, bbox_inches='tight')

fig,ax=plt.subplots(2,3,figsize=(8.5,4))

eps = [0.1,0.01]
steps = [1000,10000,100000]
nonelist = [[none_eps01_mc1000,none_eps01_mc10000,none_eps01_mc100000],[none_eps001_mc1000,none_eps001_mc10000,none_eps001_mc100000]]

for i in range(2):
    for j in range(3):
        data = nonelist[i][j]
        ax[i,j].plot(data[:,0],data[:,2],label=r'$%d$ MC sweeps, $\epsilon=%.02f$' % (steps[j], eps[i]),color='red')
        ax[i,j].legend(fontsize=6)
        ax[i,j].tick_params(axis='both', which='major', labelsize=8)
        if(i==1):
            ax[i,j].set_xlabel('iterations',fontsize=8)
        if(j==0):
            ax[i,j].set_ylabel(r'$D_2$', fontsize=8)
plt.savefig('kl_J_init=none.png', dpi=300, bbox_inches='tight')

#None init
fig,ax=plt.subplots(2,3,figsize=(8.5,4))

eps = [0.1,0.01]
steps = [1000,10000,100000]
bmlist = [[bm_eps01_mc1000,bm_eps01_mc10000,bm_eps01_mc100000],[bm_eps001_mc1000,bm_eps001_mc10000,bm_eps001_mc100000]]

for i in range(2):
    for j in range(3):
        data = bmlist[i][j]
        ax[i,j].plot(data[:,0],data[:,1],label=r'$%d$ MC sweeps, $\epsilon=%.02f$' % (steps[j], eps[i]))
        ax[i,j].legend(fontsize=6)
        ax[i,j].tick_params(axis='both', which='major', labelsize=8)
        if(i==1):
            ax[i,j].set_xlabel('iterations',fontsize=8)
        if(j==0):
            ax[i,j].set_ylabel(r'$D_1$', fontsize=8)
plt.savefig('kl_h_init=bm.png', dpi=300, bbox_inches='tight')

fig,ax=plt.subplots(2,3,figsize=(8.5,4))

eps = [0.1,0.01]
steps = [1000,10000,100000]

for i in range(2):
    for j in range(3):
        data = bmlist[i][j]
        ax[i,j].plot(data[:,0],data[:,2],label=r'$%d$ MC sweeps, $\epsilon=%.02f$' % (steps[j], eps[i]),color='red')
        ax[i,j].legend(fontsize=6)
        ax[i,j].tick_params(axis='both', which='major', labelsize=8)
        if(i==1):
            ax[i,j].set_xlabel('iterations',fontsize=8)
        if(j==0):
            ax[i,j].set_ylabel(r'$D_2$', fontsize=8)
plt.savefig('kl_J_init=bm.png', dpi=300, bbox_inches='tight')

plt.show()
