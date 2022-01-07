#Get contact scores for a given fitness landscape

import numpy as np
import sys
import math

method = sys.argv[1]

N=56
q=21

Npair = int(N*(N-1)/2)

J=np.zeros((Npair,q,q))
scores=np.zeros((N,N))
fn_vals=np.zeros((N,N))

#First change to "zero sum" gauge

print(np.sum(J))

print('Loading parameter file...')
    
with open('%s_params.txt' % (method)) as f:
  lines = [line.rstrip() for line in f]

l=0
linecount=0
for i in range(N-1):
  for j in range(i+1,N):
    for q1 in range(q):
      for q2 in range(q):
        J[l,q1,q2] = lines[linecount].split()[5]
        linecount = linecount+1
    l = l+1


print('Loaded file. Computing scores...')

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
print('before switching to zero-sum gauge:')
print(np.max(fn_vals))

l=0
for i in range(N-1):
    for j in range(i+1,N):
        for q1 in range(q):
            for q2 in range(q):
                J[l,q1,q2] = J[l,q1,q2] - np.mean(J[l,:,q2]) - np.mean(J[l,q1,:]) + np.mean(J[l,:,:])

l=0
for i in range(N-1):
  for j in range(i+1,N):
    fn_vals[i,j] = np.linalg.norm(J[l,1:,1:])
    fn_vals[j,i] = fn_vals[i,j]
    l = l+1
print('after switching to zero-sum gauge:')
print(np.max(fn_vals))

#Compute the "corrected norm"
for i in range(N-1):
  for j in range(i+1,N):
    scores[i,j] = fn_vals[i,j]-np.sum(fn_vals[i,:])*np.sum(fn_vals[:,j])/np.sum(fn_vals)
    scores[j,i] = scores[i,j]
    
print('Scores Computed. Saving to file.')
np.savetxt('%s_scores.txt' % (method), scores)
