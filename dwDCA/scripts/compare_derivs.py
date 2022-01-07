#Created by Layne Frechette, March 11th, 2021

import numpy as np
import pylab as plt
import sys
plt.rcParams.update({'font.size': 8})


none_eps01_mc1000 = np.loadtxt('data/derivs_ga_init=none_maxiter=1000_stepsize=0.100000_mcsteps=1000.txt')
none_eps01_mc10000 = np.loadtxt('data/derivs_ga_init=none_maxiter=1000_stepsize=0.100000_mcsteps=10000.txt')
none_eps01_mc100000 = np.loadtxt('data/derivs_ga_init=none_maxiter=1000_stepsize=0.100000_mcsteps=100000.txt')
none_eps001_mc1000 = np.loadtxt('data/derivs_ga_init=none_maxiter=1000_stepsize=0.010000_mcsteps=1000.txt')
none_eps001_mc10000 = np.loadtxt('data/derivs_ga_init=none_maxiter=1000_stepsize=0.010000_mcsteps=10000.txt')
none_eps001_mc100000 = np.loadtxt('data/derivs_ga_init=none_maxiter=1000_stepsize=0.010000_mcsteps=100000.txt')

bm_eps01_mc1000 = np.loadtxt('data/derivs_ga_init=bm_maxiter=1000_stepsize=0.100000_mcsteps=1000.txt')
bm_eps01_mc10000 = np.loadtxt('data/derivs_ga_init=bm_maxiter=1000_stepsize=0.100000_mcsteps=10000.txt')
bm_eps01_mc100000 = np.loadtxt('data/derivs_ga_init=bm_maxiter=1000_stepsize=0.100000_mcsteps=100000.txt')
bm_eps001_mc1000 = np.loadtxt('data/derivs_ga_init=bm_maxiter=1000_stepsize=0.010000_mcsteps=1000.txt')
bm_eps001_mc10000 = np.loadtxt('data/derivs_ga_init=bm_maxiter=1000_stepsize=0.010000_mcsteps=10000.txt')
bm_eps001_mc100000 = np.loadtxt('data/derivs_ga_init=bm_maxiter=1000_stepsize=0.010000_mcsteps=100000.txt')

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
            ax[i,j].set_ylabel(r'$\partial L/\partial h_{i}(\alpha)$', fontsize=8)
plt.savefig('deriv_h_init=none.png', dpi=300, bbox_inches='tight')

fig,ax=plt.subplots(2,3,figsize=(8.5,4))

eps = [0.1,0.01]
steps = [1000,10000,100000]
nonelist = [[none_eps01_mc1000,none_eps01_mc10000,none_eps01_mc100000],[none_eps001_mc1000,none_eps001_mc10000,none_eps001_mc100000]]

for i in range(2):
    for j in range(3):
        data = nonelist[i][j]
        ax[i,j].plot(data[:,0],data[:,3],label=r'$%d$ MC sweeps, $\epsilon=%.02f$' % (steps[j], eps[i]),color='red')
        ax[i,j].legend(fontsize=6)
        ax[i,j].tick_params(axis='both', which='major', labelsize=8)
        if(i==1):
            ax[i,j].set_xlabel('iterations',fontsize=8)
        if(j==0):
            ax[i,j].set_ylabel(r'$\partial L/\partial J_{ij}(\alpha,\beta)$', fontsize=8)
plt.savefig('derivs_J_init=none.png', dpi=300, bbox_inches='tight')

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
            ax[i,j].set_ylabel(r'$\partial L/\partial h_{i}(\alpha)$', fontsize=8)
plt.savefig('derivs_h_init=bm.png', dpi=300, bbox_inches='tight')

fig,ax=plt.subplots(2,3,figsize=(8.5,4))

eps = [0.1,0.01]
steps = [1000,10000,100000]

for i in range(2):
    for j in range(3):
        data = bmlist[i][j]
        ax[i,j].plot(data[:,0],data[:,3],label=r'$%d$ MC sweeps, $\epsilon=%.02f$' % (steps[j], eps[i]),color='red')
        ax[i,j].legend(fontsize=6)
        ax[i,j].tick_params(axis='both', which='major', labelsize=8)
        if(i==1):
            ax[i,j].set_xlabel('iterations',fontsize=8)
        if(j==0):
            ax[i,j].set_ylabel(r'$\partial L/\partial J_{ij}(\alpha,\beta)$', fontsize=8)
plt.savefig('derivs_J_init=bm.png', dpi=300, bbox_inches='tight')

plt.show()
