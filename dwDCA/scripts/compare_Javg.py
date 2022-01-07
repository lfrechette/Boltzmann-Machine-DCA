import numpy as np
import pylab as plt

ga = np.loadtxt('Javg_ga.txt')
bmstart = np.loadtxt('Javg_ga_bmstart.txt')
adapt = np.loadtxt('Javg_ga_adaptive.txt')
bm = np.loadtxt('Javg_bm.txt')

allscores = [ga,bmstart,adapt,bm]

idx = [0,1,2,5,8,10,11,18,19,20,13,7,9,15,3,4,12,14,16,17,6]
for i in range(4):
  allscores[i] = allscores[i][idx,:]
  allscores[i] = allscores[i][:,idx]
  allscores[i] = allscores[i][1:,1:]
  allscores[i] = (allscores[i]+allscores[i].T)/2
  #allscores[i] = allscores[i]-np.min(allscores[i])
  #allscores[i] = allscores[i]/np.max(allscores[i])

fig, ax = plt.subplots(2,2,figsize=(7,7))

y_label_list = ['A','C','F','I','L','M','V','W','Y','P','H','K','R','D','E','N','Q','S','T','G']
ax[0,0].set_yticks([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20])
ax[0,0].set_yticklabels(y_label_list[::-1])
ax[0,0].set_xticks([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20])
ax[0,0].set_xticklabels(y_label_list)
ax[1,0].set_yticks([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20])
ax[1,0].set_yticklabels(y_label_list[::-1])
ax[1,0].set_xticks([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20])
ax[1,0].set_xticklabels(y_label_list)
ax[0,1].set_yticks([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20])
ax[0,1].set_yticklabels(y_label_list[::-1])
ax[0,1].set_xticks([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20])
ax[0,1].set_xticklabels(y_label_list)
ax[1,1].set_yticks([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20])
ax[1,1].set_yticklabels(y_label_list[::-1])
ax[1,1].set_xticks([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20])
ax[1,1].set_xticklabels(y_label_list)

im=ax[0,0].imshow(allscores[0], extent=[0.5,20.5,0.5,20.5], cmap='RdBu_r', vmin=-0.1, vmax=0.1)
im=ax[0,1].imshow(allscores[1], extent=[0.5,20.5,0.5,20.5], cmap='RdBu_r', vmin=-0.1, vmax=0.1)
im=ax[1,0].imshow(allscores[2], extent=[0.5,20.5,0.5,20.5], cmap='RdBu_r', vmin=-0.1, vmax=0.1)
im=ax[1,1].imshow(allscores[3], extent=[0.5,20.5,0.5,20.5], cmap='RdBu_r', vmin=-0.1, vmax=0.1)
ax[0,0].set_title(r'fixed step size')
ax[0,1].set_title(r'adaptive, Rang. start')
ax[1,0].set_title(r'adaptive step size')
ax[1,1].set_title(r'Ranganathan')

cbar_ax = fig.add_axes([0.925,0.15,0.025,0.7])
fig.colorbar(im,cax=cbar_ax)

plt.savefig('compare_Javg.png', dpi=300, bbox_inches='tight')

plt.show()
