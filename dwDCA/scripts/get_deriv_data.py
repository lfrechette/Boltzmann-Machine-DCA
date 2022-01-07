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

scratch_dir = '/scratch/frechettelb/protein_evolution/dwDCA/data/%s_init=%s_maxiter=%d_stepsize=%f_mcsteps=%d' % (msa, input_name, maxiter, eps, mcsteps)
out_dir = 'data/%s_init=%s_maxiter=%d_stepsize=%f_mcsteps=%d' % (msa, input_name, maxiter, eps, mcsteps)

iters = np.arange(maxiter)
dh_norm = np.zeros(maxiter)
dh_max = np.zeros(maxiter)
dJ_norm = np.zeros(maxiter)
dJ_max = np.zeros(maxiter)

N = 0
q = 0
Npair = 0

for i in range(maxiter):
    print(i)
    if(os.path.isfile(scratch_dir + '/dh1_%d.txt' % i)):
        with open(scratch_dir + '/dh1_%d.txt' % i) as file:
            line = next(file).strip()
            line = next(file).strip()
            line = line.split( )
            N = int(line[0])
            q = int(line[1])
        Npair = int(N*(N-1)/2)
        dh = np.loadtxt(scratch_dir + '/dh1_%d.txt' % i, skiprows=2)
        dJ = np.loadtxt(scratch_dir + '/dJ1_%d.txt' % i, skiprows=2)
        for j in range(N):
            dh_norm[i] += LA.norm(dh[j,:])
        dh_norm[i] = dh_norm[i]/N
        for j in range(Npair):
            dJmat = dJ[j*q:(j+1)*q,:]
            dJ_norm[i] += LA.norm(dJmat)
        dJ_norm[i] = dJ_norm[i]/Npair
        dh_max[i] = np.max(np.abs(dh))
        dJ_max[i] = np.max(np.abs(dJ))
        print('%f %f %f %f' % (dh_norm[i], dh_max[i], dJ_norm[i], dJ_max[i]))

np.savetxt(out_dir+'/derivs.txt', np.c_[iters, dh_norm, dh_max,dJ_norm,dJ_max])

#plt.figure()
#plt.imshow(scores,vmin=0,vmax=1.0)
#plt.savefig('%s_scores.png' % family, dpi=300, bbox_inches='tight')
#plt.show()
