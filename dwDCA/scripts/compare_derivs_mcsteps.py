#Created by Layne Frechette, March 11th, 2021

import numpy as np
import pylab as plt
import sys

family = sys.argv[1]
init = sys.argv[2]
maxiter = int(sys.argv[3])
eps = float(sys.argv[4])
#mcsteps = int(sys.argv[5])

data1 = np.loadtxt('data/derivs_%s_init=%s_maxiter=%d_stepsize=%f_mcsteps=1000.txt' % (family,init,maxiter,eps))
data2 = np.loadtxt('data/derivs_%s_init=%s_maxiter=%d_stepsize=%f_mcsteps=10000.txt' % (family,init,maxiter,eps))
data3 = np.loadtxt('data/derivs_%s_init=%s_maxiter=%d_stepsize=%f_mcsteps=100000.txt' % (family,init,maxiter,eps))

###
plt.figure()
plt.plot(data1[:,0],data1[:,1],label=r'$10^3 MC sweeps')
plt.plot(data2[:,0],data2[:,1],label=r'$10^4 MC sweeps')
plt.plot(data3[:,0],data3[:,1],label=r'$10^5 MC sweeps')
plt.savefig('data/dh_norm_%s_init=%s_maxiter=%d_stepsize=%f_compare_mcsteps.png' % (family,init,maxiter,eps), dpi=300, bbox_inches='tight')

plt.figure()
plt.plot(data1[:,0],data1[:,2],label=r'$10^3 MC sweeps')
plt.plot(data2[:,0],data2[:,2],label=r'$10^4 MC sweeps')
plt.plot(data3[:,0],data3[:,2],label=r'$10^5 MC sweeps')
plt.savefig('data/dh_max_%s_init=%s_maxiter=%d_stepsize=%f_compare_mcsteps.png' % (family,init,maxiter,eps), dpi=300, bbox_inches='tight')

plt.figure()
plt.plot(data1[:,0],data1[:,3],label=r'$10^3 MC sweeps')
plt.plot(data2[:,0],data2[:,3],label=r'$10^4 MC sweeps')
plt.plot(data3[:,0],data3[:,3],label=r'$10^5 MC sweeps')
plt.savefig('data/dJ_norm_%s_init=%s_maxiter=%d_stepsize=%f_compare_mcsteps.png' % (family,init,maxiter,eps), dpi=300, bbox_inches='tight')

plt.figure()
plt.plot(data1[:,0],data1[:,4],label=r'$10^3 MC sweeps')
plt.plot(data2[:,0],data2[:,4],label=r'$10^4 MC sweeps')
plt.plot(data3[:,0],data3[:,4],label=r'$10^5 MC sweeps')
plt.savefig('data/dJ_max_%s_init=%s_maxiter=%d_stepsize=%f_compare_mcsteps.png' % (family,init,maxiter,eps), dpi=300, bbox_inches='tight')

plt.show()
