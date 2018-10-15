# -*- coding: utf-8 -*-
"""
Created on Fri Oct 12 14:39:05 2018

@author: gregz
"""

import numpy as np
import slurmfile
from astropy.io import fits

filenames = [line.rstrip('\n').split()
             for line in open('/work/03730/gregz/maverick/ifuslots.dat', 'r')]


F = [filename[1] for filename in filenames]

chunks = np.array_split(F, len(F) / 20 + 1)
for j, chunk in enumerate(chunks):
    n = len(chunk)
    name = 'rtwi_%s_%i' % ('virus', j+1)
    f = open(name+'.slurm', 'w')
    s = slurmfile.slurmstring % (n, '%j', name)
    f.write(s)
    f.close()
    f = open(name, 'w')
    s = []
    for ifuslot in chunk:
        cmd = ('python /work/03730/gregz/maverick/Panacea/panacea2.py '
               ' --ifuslot %03d -rt -td 20181003 -to 2 -te 1' % int(ifuslot))
        s.append(cmd)
    f.write('\n'.join(s))
    f.close()
    print('sbatch %s.slurm' % name)
