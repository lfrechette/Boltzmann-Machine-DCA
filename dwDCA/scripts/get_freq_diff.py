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
p1_norm = np.zeros(maxiter)
p2_norm = np.zeros(maxiter)

msa_1p = np.loadtxt(scratch_dir + '/stat_align_1p.txt', skiprows=2)
msa_2p = np.loadtxt(scratch_dir + '/stat_align_2p.txt', skiprows=2)

for i in range(maxiter):
    print(i)
    if(os.path.isfile(scratch_dir + '/stat_MC_1p_%d.txt' % i)):
        mc_1p = np.loadtxt(scratch_dir + '/stat_MC_1p_%d.txt' % (i+1), skiprows=2)
        mc_2p = np.loadtxt(scratch_dir + '/stat_MC_2p_%d.txt' % (i+1), skiprows=2)
        diff_1p = msa_1p-mc_1p
        diff_2p = msa_2p-mc_2p
        p1_norm[i] = np.sum(diff_1p**2)
        p2_norm[i] = np.sum(diff_2p**2) 
        print('%f %f' % (p1_norm[i], p2_norm[i]))

np.savetxt('data/freq_diff_%s_%s_%d_%f_%d_%f.txt' % (msa, input_name, maxiter, eps, mcsteps, lam), np.c_[iters, p1_norm, p2_norm])

