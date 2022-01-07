import numpy as np
import pylab as plt

none_eps01_mc1000 = np.loadtxt('data/scores_msa_ga.fas_init=none_maxiter=1000_stepsize=0.100000_mcsteps=1000.txt')
bm = np.loadtxt('bm_scores.txt')

max1 = np.amax(none_eps01_mc1000)
max2 = np.amax(bm)
ulim=max(max1,max2)
print(ulim)

###
fig, ax = plt.subplots(1,2,figsize=(7,3))

im=ax[0].imshow(none_eps01_mc1000, vmin=-ulim, vmax=ulim, cmap='RdBu_r')
im=ax[1].imshow(bm, vmin=-ulim, vmax=ulim, cmap='RdBu_r')
ax[0].set_title(r'$\epsilon=0.1$, 1000 MC steps, 1000 iters', fontsize=12)
ax[1].set_title(r'Ranganathan', fontsize=12)

cbar_ax = fig.add_axes([0.925,0.15,0.025,0.7])
fig.colorbar(im,cax=cbar_ax)

plt.savefig('compare_scores_new.png', dpi=300, bbox_inches='tight')

plt.show()
