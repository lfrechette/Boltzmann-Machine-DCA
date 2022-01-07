#Get contact scores for a given fitness landscape

import numpy as np
import sys
import math

method = sys.argv[1]

N=56
q=21

Npair = int(N*(N-1)/2)

J=np.zeros((Npair,q,q))
h=np.zeros((N,q))

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
for i in range(N):
  for q1 in range(q):
    h[i,q1] = lines[linecount].split()[3]
    linecount = linecount+1

print(h)


print('Loaded file. Computing scores...')

np.savetxt('%s_h.txt' % (method), h)
