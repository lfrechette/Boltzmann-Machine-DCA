#Get contact scores for a given fitness landscape

import numpy as np
import sys
import math

msa = sys.argv[1]
init = sys.argv[2]
maxiter = int(sys.argv[3])
eps = float(sys.argv[4])
mcsteps = int(sys.argv[5])
lam = float(sys.argv[6])

scratch_dir = '/scratch/frechettelb/protein_evolution/dwDCA/data/%s_init=%s_maxiter=%d_stepsize=%f_mcsteps=%d_lambda=%f' % (msa, init, maxiter, eps, mcsteps, lam)

N=56
q=21

Npair = int(N*(N-1)/2)

for m in range(maxiter):

  print('Iteration %d' % m)
  J=np.zeros((Npair,q,q))
  scores=np.zeros((N,N))
  fn_vals=np.zeros((N,N))

  with open(scratch_dir + '/J1_%d.txt' % m) as f:
    lines = [line.rstrip() for line in f]

  l=0
  linecount=2
  for i in range(N-1):
    for j in range(i+1,N):
      for q1 in range(q):
        for q2 in range(q):
          J[l,q1,q2] = lines[linecount].split()[q2]
        linecount = linecount+1
      l = l+1

  #Compute scores...
  #First convert to zero-sum gauge
  l=0
  for i in range(N-1):
    for j in range(i+1,N):
      for q1 in range(q):
        for q2 in range(q):
          J[l,q1,q2] = J[l,q1,q2] - np.mean(J[l,:,q2]) - np.mean(J[l,q1,:]) + np.mean(J[l,:,:])
  #Then compute Froebenius norm of J matrices
  l=0
  for i in range(N-1):
    for j in range(i+1,N):
      fn_vals[i,j] = np.linalg.norm(J[l,1:,1:])
      fn_vals[j,i] = fn_vals[i,j]
      l = l+1

  #Compute the "corrected norm"
  for i in range(N-1):
    for j in range(i+1,N):
      scores[i,j] = fn_vals[i,j]-np.sum(fn_vals[i,:])*np.sum(fn_vals[:,j])/np.sum(fn_vals)
      scores[j,i] = scores[i,j]
    
  np.savetxt(scratch_dir + '/scores_%d.txt' % m, scores)
