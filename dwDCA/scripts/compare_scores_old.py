import numpy as np
import pylab as plt

ga = np.loadtxt('ga_scores.txt')
#bmstart = np.loadtxt('ga_bmstart_scores.txt')
#adapt = np.loadtxt('ga_adaptive_scores.txt')
ga_01 = np.loadtxt('ga_01_scores.txt')
ga_01_200 = np.loadtxt('ga_01_200_scores.txt')
bm = np.loadtxt('bm_scores.txt')

allscores = [ga,ga_01,ga_01_200,bm]
maxes = [np.max(ga),np.max(ga_01),np.max(ga_01_200),np.max(bm)]

for i in range(4):
  allscores[i] = allscores[i]-np.min(allscores[i])
  allscores[i] = allscores[i]/np.max(allscores[i])

###
fig, ax = plt.subplots(2,2,figsize=(7,7))

im=ax[0,0].imshow(allscores[0], vmin=0, vmax=1)
im=ax[0,1].imshow(allscores[1], vmin=0, vmax=1)
im=ax[1,0].imshow(allscores[2], vmin=0, vmax=1)
im=ax[1,1].imshow(allscores[3], vmin=0, vmax=1)
ax[0,0].set_title(r'$\epsilon=0.01$, 100 steps, $J_{max}=$%.02f' % maxes[0], fontsize=12)
ax[0,1].set_title(r'$\epsilon=0.1$, 100 steps, $J_{max}=$%.02f' % maxes[1], fontsize=12)
ax[1,0].set_title(r'$\epsilon=0.1$, 200 steps, $J_{max}=$%.02f' % maxes[2], fontsize=12)
ax[1,1].set_title(r'Ranganathan, $J_{max}=$%.02f' % maxes[3], fontsize=12)

cbar_ax = fig.add_axes([0.925,0.15,0.025,0.7])
fig.colorbar(im,cax=cbar_ax)

plt.savefig('compare_scores.png', dpi=300, bbox_inches='tight')

plt.show()
