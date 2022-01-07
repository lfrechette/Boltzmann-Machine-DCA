#Created by Layne Frechette, March 11th, 2021

import numpy as np
import pylab as plt
import sys
plt.rcParams.update({'font.size': 8})

reg0 = np.loadtxt('data/freq_diff_rang_lambda=0.0.txt')
reg001= np.loadtxt('data/freq_diff_rang_lambda=0.01.txt')
reg1= np.loadtxt('data/freq_diff_rang_lambda=1.0.txt')


fig,ax=plt.subplots(1,3,figsize=(7,2.5))

regs=[0.0,0.01,1.0]
mylist = [reg0,reg001,reg1]

###
for j in range(3):
    data = mylist[j]
    ax[j].plot(data[:,0],data[:,1])
    ax[j].tick_params(axis='both', which='major', labelsize=8)
    ax[j].set_xlabel('iterations',fontsize=8)
    if(j==0):
        ax[j].set_ylabel(r'$\sum_{i,\alpha}(f_i^{\text{MSA}}(\alpha)-f_i^{\text{MC}}(\alpha))^2$', fontsize=8)
    ax[j].set_title(r'$\lambda=%.02f$' % regs[j])
    ax[j].set_xlim([0,100])
plt.savefig('freq_diff_h_rang_regs.png', dpi=300, bbox_inches='tight')


fig,ax=plt.subplots(1,3,figsize=(7,2.5))
for j in range(3):
    data = mylist[j]
    ax[j].plot(data[:,0],data[:,2],color='red')
    ax[j].tick_params(axis='both', which='major', labelsize=8)
    ax[j].set_xlabel('iterations',fontsize=8)
    if(j==0):
        ax[j].set_ylabel(r'$\sum_{i,j,\alpha,\beta}(f_{ij}^{\text{MSA}}(\alpha,\beta)-f_{ij}^{\text{MC}}(\alpha,\beta))^2$', fontsize=8)
    ax[j].set_title(r'$\lambda=%.02f$' % regs[j])
    ax[j].set_xlim([0,100])
plt.savefig('freq_diff_J_rang_regs.png', dpi=300, bbox_inches='tight')

'''
for j in range(3):
    data = mylist[j]
    ax[j].plot(data[:,0],data[:,1])
    ax[j].tick_params(axis='both', which='major', labelsize=8)
    ax[j].set_xlabel('iterations',fontsize=8)
    if(j==0):
        ax[j].set_ylabel(r'$\sum_{i,\alpha}(f_i^{\text{MSA}}(\alpha)-f_i^{\text{MC}}(\alpha))^2$', fontsize=8)
    ax[j].set_title(r'$\lambda=%.02f$' % regs[j])
    ax[j].set_xlim([0,100])
    ax[j].set_ylim([0,15])
plt.savefig('freq_diff_h_regs_zoom.png', dpi=300, bbox_inches='tight')


fig,ax=plt.subplots(1,3,figsize=(7,2.5))
for j in range(3):
    data = mylist[j]
    ax[j].plot(data[:,0],data[:,2],color='red')
    ax[j].tick_params(axis='both', which='major', labelsize=8)
    ax[j].set_xlabel('iterations',fontsize=8)
    if(j==0):
        ax[j].set_ylabel(r'$\sum_{i,j,\alpha,\beta}(f_{ij}^{\text{MSA}}(\alpha,\beta)-f_{ij}^{\text{MC}}(\alpha,\beta))^2$', fontsize=8)
    ax[j].set_title(r'$\lambda=%.02f$' % regs[j])
    ax[j].set_xlim([0,100])
    ax[j].set_ylim([0,350])
plt.savefig('freq_diff_J_regs_zoom.png', dpi=300, bbox_inches='tight')
'''
plt.show()
