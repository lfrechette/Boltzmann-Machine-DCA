#Created by Layne Frechette, March 11th, 2021

import numpy as np
import pylab as plt
import sys

family = sys.argv[1]
init = sys.argv[2]
maxiter = int(sys.argv[3])
eps = float(sys.argv[4])
mcsteps = int(sys.argv[5])
lam = float(sys.argv[6])

msa = np.loadtxt('data/%s_init=%s_maxiter=%d_stepsize=%f_mcsteps=%d_lambda=%f/stat_align_1p.txt' % (family, init, maxiter, eps, mcsteps, lam),skiprows=2)

ulim=1.0

for i in range(maxiter):
    data = np.loadtxt('data/%s_init=%s_maxiter=%d_stepsize=%f_mcsteps=%d_lambda=%f/stat_MC_1p_%d.txt' % (family, init, maxiter, eps, mcsteps, lam, i),skiprows=2)
    plt.figure()
    plt.imshow(data-msa,vmin=-ulim,vmax=ulim,cmap='RdBu_r')
    plt.savefig('data/images/freq_diff_map_%04d.png' % i)
    plt.close()

