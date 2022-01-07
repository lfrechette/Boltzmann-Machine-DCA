#Created by Layne Frechette, March 12th, 2021

import numpy as np
import pylab as plt
import sys

data = np.loadtxt('data/freq_diff_rang_mimic.txt')

###
plt.figure()
plt.plot(data[:,0],data[:,1])
'''
plt.plot(data[:,0],data[:,3])
plt.plot(data[:,0],data[:,5])
plt.plot(data[:,0],data[:,7])
plt.plot(data[:,0],data[:,9])
plt.plot(data[:,0],data[:,11])
plt.plot(data[:,0],data[:,13])
plt.plot(data[:,0],data[:,15])
plt.plot(data[:,0],data[:,17])
plt.plot(data[:,0],data[:,19])
'''
plt.xlabel('iterations')
plt.ylabel(r'$\sum_{i,\alpha}(f_i^{\text{MSA}}(\alpha)-f_i^{\text{MC}}(\alpha))^2$')
plt.title('Ranganathan code')
plt.savefig('data/freq_diff_1p_rang_mimic.png', dpi=300, bbox_inches='tight')

plt.figure()
plt.plot(data[:,0],data[:,2],color='red')
'''
plt.plot(data[:,0],data[:,4])
plt.plot(data[:,0],data[:,6])
plt.plot(data[:,0],data[:,8])
plt.plot(data[:,0],data[:,10])
plt.plot(data[:,0],data[:,12])
plt.plot(data[:,0],data[:,14])
plt.plot(data[:,0],data[:,16])
plt.plot(data[:,0],data[:,18])
plt.plot(data[:,0],data[:,20])
'''
plt.xlabel('iterations')
plt.ylabel(r'$\sum_{i,j,\alpha,\beta}(f_{ij}^{\text{MSA}}(\alpha,\beta)-f_{ij}^{\text{MC}}(\alpha,\beta))^2$')
plt.title('Ranganathan code')
plt.savefig('data/freq_diff_2p_rang_mimic.png', dpi=300, bbox_inches='tight')

plt.show()
