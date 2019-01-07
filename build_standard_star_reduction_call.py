# -*- coding: utf-8 -*-
"""
Created on Thu Oct  4 13:30:06 2018

@author: gregz
"""
import numpy as np
import os.path as op
import sys
import slurmfile


standard_names = ['HD_19445', 'SA95-42', 'GD50', 'G191B2B', 'FEIGE_25',
                  'HILTNER_600', 'G193-74', 'PG0823+546', 'HD_84937',
                  'GD108', 'FEIGE_34', 'HD93521', 'GD140', 'HZ_21',
                  'FEIGE_66', 'FEIGE_67', 'G60-54', 'HZ_44', 'GRW+70_5824',
                  'BD+26+2606', 'BD+33_2642', 'G138-31', 'WOLF_1346',
                  'BD_+17_4708', 'FEIGE_110', 'GD248', 'HZ_4',
                  'BD+40_4032', 'HILTNER_102']

object_table = [line.rstrip('\n').split() for line in open(sys.argv[1])]
filenames = []
objects = []
dates = []
standard_call = []
basecall = ('python /work/03730/gregz/maverick/Panacea/full_lrs2_reduction.py '
            '-d %s -s %s -o %s')
slot = sys.argv[2]
if slot == '056':
    names = "\"uv,orange\""
else:
    names = "\"red,farred\""
for _object in object_table:
    if slot not in _object[1]:
        continue
    filename = _object[0]
    objectname = _object[1].split('_%s' % slot)[0]
    exptime = float(_object[2])
    if objectname == '':
        continue
    for standard in standard_names:
        if objectname.lower() in standard.lower():
            obsid = op.basename(op.dirname(op.dirname(op.dirname(filename)))).split('lrs2')[1]
            date = op.basename(op.dirname(op.dirname(op.dirname(op.dirname(op.dirname(filename))))))
            standard_call.append(basecall %
                                 (date, names, objectname))

for f, basename in zip([standard_call], ['rstan']):
    chunks = np.array_split(f, len(f) / 20 + 1)
    for j, chunk in enumerate(chunks):
        n = len(chunk)
        name = basename+'_%s_%i' % (slot, j+1)
        f = open(name+'.slurm', 'w')
        s = slurmfile.slurmstring % (n, '%j', name)
        f.write(s)
        f.close()
        f = open(name, 'w')
        s = []
        for call in chunk:
            s.append(call)
        f.write('\n'.join(s))
        f.close()
        print('sbatch %s.slurm' % name)            