#Created by Layne Frechette, March 12th, 2021

import numpy as np
import pylab as plt
import sys
import os.path
from numpy import linalg as LA

'''
msa = sys.argv[1]
input_name = sys.argv[2]
maxiter = int(sys.argv[3])
eps = float(sys.argv[4])
mcsteps = int(sys.argv[5])
'''

scratch_dir = sys.argv[1]

maxiter=1000

#scratch_dir = 'ranganathan'

iters = np.arange(maxiter)#*20+20
p1_norm = np.zeros(maxiter)
p2_norm = np.zeros(maxiter)
p1_norm_thresh = np.zeros((9,maxiter))
p2_norm_thresh = np.zeros((9,maxiter))

msa_1p = np.loadtxt(scratch_dir + '/stat_align_1p.txt', skiprows=2)
msa_2p = np.loadtxt(scratch_dir + '/stat_align_2p.txt', skiprows=2)

print(iters)

threshes=[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]
cnt=0
for i in range(maxiter):
    print(i)
    mc_1p = np.loadtxt(scratch_dir + '/stat_MC_1p_%d.txt' % (i+1), skiprows=2)
    mc_2p = np.loadtxt(scratch_dir + '/stat_MC_2p_%d.txt' % (i+1), skiprows=2)
    diff_1p = msa_1p-mc_1p
    diff_2p = msa_2p-mc_2p
    p1_norm[cnt] = np.sum(diff_1p**2)
    p2_norm[cnt] = np.sum(diff_2p**2) 
    for j in range(len(threshes)):
        thresh_indices_1 = np.where(msa_1p>threshes[j])
        thresh_indices_2 = np.where(msa_2p>threshes[j])
        p1_norm_thresh[j,cnt] = np.sum(diff_1p[thresh_indices_1]**2)
        p2_norm_thresh[j,cnt] = np.sum(diff_2p[thresh_indices_2]**2)
    print('%f %f' % (p1_norm[cnt], p2_norm[cnt]))
    cnt+=1

np.savetxt(scratch_dir + '/freq_diff.txt', np.c_[iters, p1_norm, p2_norm, p1_norm_thresh[0,:], p2_norm_thresh[0,:], p1_norm_thresh[1,:], p2_norm_thresh[1,:], p1_norm_thresh[2,:], p2_norm_thresh[2,:], p1_norm_thresh[3,:], p2_norm_thresh[3,:], p1_norm_thresh[4,:], p2_norm_thresh[4,:], p1_norm_thresh[5,:], p2_norm_thresh[5,:], p1_norm_thresh[6,:], p2_norm_thresh[6,:], p1_norm_thresh[7,:], p2_norm_thresh[7,:], p1_norm_thresh[8,:], p2_norm_thresh[8,:]])

