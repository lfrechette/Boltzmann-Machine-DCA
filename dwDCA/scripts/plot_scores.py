#Created by Layne Frechette, Feb. 25th, 2021

import numpy as np
import pylab as plt
import sys

family = sys.argv[1]

scores = np.loadtxt('%s_scores.txt' % family)

scores = scores-np.min(scores)
scores = scores/np.max(scores)

plt.figure()
plt.imshow(scores,vmin=0,vmax=1.0)
plt.savefig('%s_scores.png' % family, dpi=300, bbox_inches='tight')
plt.show()
