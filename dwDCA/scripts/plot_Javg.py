#Plot <J> for different methods for GA
#Created by Layne Frechette on Feb. 18, 2021

import numpy as np
import pylab as plt
import matplotlib as mpl

mpl.rc('xtick', labelsize=7)
mpl.rc('ytick', labelsize=7)

N=56
q=21

Npair = int(N*(N-1)/2)

bm = np.loadtxt('Javg_ga.txt')

fig, ax = plt.subplots(1,1,figsize=(3,3))

myvmin=-0.001
myvmax=0.001

idx = [0,1,2,5,8,10,11,18,19,20,13,7,9,15,3,4,12,14,16,17,6]
bm = bm[idx,:]
bm = bm[:,idx]
bm = bm[1:,1:]
bm = (bm+bm.T)/2

im=ax.imshow(bm, extent=[0.5,20.5,0.5,20.5], cmap='RdBu_r')#, vmin=myvmin, vmax=myvmax)

ax.set_ylabel(r'amino acid $\alpha$')
ax.set_xlabel(r'amino acid $\beta$')

y_label_list = ['A','C','F','I','L','M','V','W','Y','P','H','K','R','D','E','N','Q','S','T','G']
ax.set_yticks([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20])
ax.set_yticklabels(y_label_list[::-1])
ax.set_xticks([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20])
ax.set_xticklabels(y_label_list)

cbar_ax = fig.add_axes([0.925, 0.15, 0.025, 0.7])
fig.colorbar(im, cax=cbar_ax)

plt.savefig('Javg.png', dpi=300, bbox_inches='tight')

plt.show()
