#Created by Layne Frechette, March 11th, 2021

import numpy as np
import pylab as plt
import sys

family = sys.argv[1]
init = sys.argv[2]
maxiter = int(sys.argv[3])
eps = float(sys.argv[4])
mcsteps = int(sys.argv[5])

data = np.loadtxt('data/derivs_%s_init=%s_maxiter=%d_stepsize=%f_mcsteps=%d.txt' % (family,init,maxiter,eps,mcsteps))

###
plt.figure()
plt.plot(data[:,0],data[:,1])
plt.savefig('data/derivs_%s_init=%s_maxiter=%d_stepsize=%f_mcsteps=%d.png' % (family,init,maxiter,eps,mcsteps), dpi=300, bbox_inches='tight')

###
plt.figure()



plt.show()
