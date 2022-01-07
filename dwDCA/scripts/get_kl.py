#Created by Layne Frechette, March 11th, 2021

import numpy as np
import pylab as plt
import sys
import os.path
from numpy import linalg as LA

msa = sys.argv[1]
input_name = sys.argv[2]
maxiter = int(sys.argv[3])
eps = float(sys.argv[4])
mcsteps = int(sys.argv[5])
lam = float(sys.argv[6])

scratch_dir = '/scratch/frechettelb/protein_evolution/dwDCA/data/%s_init=%s_maxiter=%d_stepsize=%f_mcsteps=%d_lambda=%f' % (msa, input_name, maxiter, eps, mcsteps, lam)
out_dir = 'data/%s_init=%s_maxiter=%d_stepsize=%f_mcsteps=%d_lambda=%f' % (msa, input_name, maxiter, eps, mcsteps, lam)

iters = np.arange(maxiter)
D1 = np.zeros(maxiter)
D2 = np.zeros(maxiter)

msa_1p = np.loadtxt(scratch_dir + '/stat_align_1p.txt', skiprows=2)
msa_2p = np.loadtxt(scratch_dir + '/stat_align_2p.txt', skiprows=2)

tiny=1e-5
print(msa_1p.shape)
L=msa_1p.shape[0]
q=msa_1p.shape[1]

for i in range(maxiter):
    print(i)
    if(os.path.isfile(scratch_dir + '/stat_MC_1p_%d.txt' % i)):
        mc_1p = np.loadtxt(scratch_dir + '/stat_MC_1p_%d.txt' % i, skiprows=2)
        mc_2p = np.loadtxt(scratch_dir + '/stat_MC_2p_%d.txt' % i, skiprows=2)
        #Compute KL divergence
        D1[i] = np.sum(np.multiply(msa_1p,np.log(np.divide(msa_1p+tiny,mc_1p+tiny))))/L
        D2[i] = np.sum(np.multiply(msa_2p,np.log(np.divide(msa_2p+tiny,mc_2p+tiny))))*(2/(L*(L-1)))
        print('%f %f' % (D1[i], D2[i]))

np.savetxt('data/kl_%s_%s_%d_%f_%d_%f.txt' % (msa, input_name, maxiter, eps, mcsteps, lam), np.c_[iters, D1, D2])

