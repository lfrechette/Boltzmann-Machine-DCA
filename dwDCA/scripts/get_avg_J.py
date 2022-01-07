#Get average J over top N/2 GA predicted contacts for different methods
#Created by Layne Frechette, Feb. 18th, 2021

import numpy as np
import pylab as plt
import sys

N = 56
Npair = int(N*(N-1)/2)
J=np.zeros((Npair,21,21))
Javg=np.zeros((21,21))

method = sys.argv[1]


with open('%s_params.txt' % method) as f:
    lines = [line.rstrip() for line in f]

l=0
linecount=0
for i in range(N-1):
    for j in range(i+1,N):
        for q1 in range(21):
            for q2 in range(21):
                J[l,q1,q2] = lines[linecount].split()[5]
                linecount = linecount+1
        l = l+1

l=0
for i in range(N-1):
    for j in range(i+1,N):
        for q1 in range(21):
            for q2 in range(21):
                J[l,q1,q2] = J[l,q1,q2] - np.mean(J[l,:,q2]) - np.mean(J[l,q1,:]) + np.mean(J[l,:,:])

scores = np.loadtxt('%s_scores.txt' % method)




#Compute predicted contacts
cutoff = np.partition(scores.flatten(),-N)[-N]
pred=np.where(scores>cutoff,1,0)

l=0
for i in range(N-1):
    for j in range(i+1,N):
        if(pred[i,j]==1):
            Javg += J[l,:,:]
        l=l+1
Javg = Javg/(N)

np.savetxt('Javg_%s.txt' % method, Javg)

plt.figure()
plt.imshow(Javg[1:,1:])
plt.show()
