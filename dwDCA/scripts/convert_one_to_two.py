import numpy as np
import sys

fam1 = sys.argv[1]
fam2 = sys.argv[2]

with open('params/%s_params.txt' % fam1) as f:
    lines1 = f.readlines()

with open('params/%s_params.txt' % fam2) as f:
    lines2 = f.readlines()


with open('params/%s%s_2well_params.txt' % (fam1, fam2), 'w') as f:
    for line in lines1:
        f.write((line.replace('J', 'J1')).replace('h','h1'))
    for line in lines2:
        f.write((line.replace('J', 'J2')).replace('h','h2'))
