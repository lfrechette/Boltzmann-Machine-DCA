#Created by Layne Frechette, March 30th, 2021

import numpy as np
import pylab as plt
import sys
import os.path
import glob
from numpy import linalg as LA

scratch_dir = '/scratch/frechettelb/protein_evolution/dwDCA/data/'

folder = sys.argv[1]

scratch_dir = scratch_dir + folder

maxiter = len(glob.glob(scratch_dir + '/stat_MC_1p*.txt'))-2
print(maxiter)

iters = np.arange(maxiter)
p1_norm = np.zeros(maxiter)
p2_norm = np.zeros(maxiter)
p1_max = np.zeros(maxiter)
p2_max = np.zeros(maxiter)

msa_1p = np.loadtxt(scratch_dir + '/stat_align_1p.txt', skiprows=2)
msa_2p = np.loadtxt(scratch_dir + '/stat_align_2p.txt', skiprows=2)

N=56
q=21
Npair=N*(N-1)/2

for i in range(maxiter):
    print(i)
    if(os.path.isfile(scratch_dir + '/stat_MC_1p_%d.txt' % i)):
        mc_1p = np.loadtxt(scratch_dir + '/stat_MC_1p_%d.txt' % (i+1), skiprows=2)
        mc_2p = np.loadtxt(scratch_dir + '/stat_MC_2p_%d.txt' % (i+1), skiprows=2)
        diff_1p = msa_1p-mc_1p
        diff_2p = msa_2p-mc_2p
        p1_norm[i] = np.sqrt(np.sum(diff_1p**2)/(N*q))
        p2_norm[i] = np.sqrt(np.sum(diff_2p**2)/(Npair*q**2))
        p1_max[i] = np.max(np.abs(diff_1p))
        p2_max[i] = np.max(np.abs(diff_2p))
        print('%f %f %f %f' % (p1_norm[i], p2_norm[i], p1_max[i], p2_max[i]))

out_dir = 'data/' + folder + '/freq_diff.txt'
np.savetxt(out_dir, np.c_[iters, p1_norm, p2_norm, p1_max, p2_max])

