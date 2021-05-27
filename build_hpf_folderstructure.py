#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 18 09:11:02 2019

@author: gregz
"""

import os.path as op
import subprocess
import numpy as np
import glob
import re
import sys
from astropy.table import Table
from datetime import datetime, timedelta
from distutils.dir_util import mkpath
from input_utils import setup_logging

basedir = '/work2/03946/hetdex/maverick/HPF_Processed'
reducdir = '/work2/03946/hetdex/maverick/HPF'

log = setup_logging('hpf_data')

SlopeImages = op.join(basedir, 'SlopeImages')
WaveCal = op.join(basedir, 'WavelengthCalibrated_e2ds')
Manifests = op.join(basedir, 'Manifests')

# Get all dates in SlopeImages folder
dates = [op.basename(fn) 
         for fn in sorted(glob.glob(op.join(SlopeImages, '*')))]

if len(sys.argv) > 1:
    dates = [date for date in dates if int(date) > int(sys.argv[1])]

def create_table(filename):
    headers = [line.rstrip('\n').split() for line in open(filename)]
    headers0 = headers[0]
    init_string = [line.rstrip('\n') for line in open(filename)]
    init_string0 = init_string[0]
    indices = [0]
    biglist = []
    for h in headers0:
        indices.append(re.search(h, init_string0).end())
        biglist.append([])
    for k in np.arange(1, len(init_string)):
        string = init_string[k]
        for i in np.arange(len(indices)-1):
            biglist[i].append(string[indices[i]:indices[i+1]].replace(' ', ''))
    t = Table(biglist, names=headers0)
    return t

for date in dates:
    log.info('Working on %s.' % date)
    fns = sorted(glob.glob(op.join(SlopeImages, date, '*')))
    try:
        T1 = create_table(op.join(Manifests, 'hpf_%s.list' % date))
    except:
        datec_ = datetime(int(date[:4]), int(date[4:6]), int(date[6:]))
        daten_ = datec_ + timedelta(days=1)
        daten = '%04d%02d%02d' % (daten_.year, daten_.month, daten_.day)
        try:
            T2 = create_table(op.join(Manifests, 'hpf_%s.list' % daten))
        except:
            continue
        T1 = T2
    datec_ = datetime(int(date[:4]), int(date[4:6]), int(date[6:]))
    daten_ = datec_ + timedelta(days=1)
    daten = '%04d%02d%02d' % (daten_.year, daten_.month, daten_.day)
    try:
        T2 = create_table(op.join(Manifests, 'hpf_%s.list' % daten))
    except:
        T2 = T1
    for fn in fns:
        obs = op.basename(fn)
        T = T1
        sel = T1['UT-Date'] == date
        sel2 = T1['ObsNum'] == obs
        if (sel * sel2).sum() == 0:
            sel = T2['UT-Date'] == date
            sel2 = T2['ObsNum'] == obs
            T = T2
        ind = np.where(sel*sel2)[0]
        if len(ind) == 0:
            log.warning('Could not find %s and %s in manifest' % (date, obs))
            continue
        ind = ind[0]
        name = T['Frame'][ind]
        oname = name.replace('.fits', '.optimal.fits')
        if T['ObsType'][ind] == 'Cal':
            folder = op.join(reducdir, 'CALS')
        if T['ObsType'][ind] == 'Sci':
            if T['QProg'][ind] != '':
                folder = op.join(reducdir, T['QProg'][ind])
            else:
                folder = op.join(reducdir, 'ORPHANS')
        mkpath(folder)
        loc = op.join(SlopeImages, date, obs, name)
        oloc = op.join(WaveCal, date, obs, oname)
        sname = op.join(folder, name[:-5] + '_' + obs + '.fits')
        soname = op.join(folder, name[:-5] + '_' + obs + '.optimal.fits')
        if op.exists(loc):
            cmd = 'cp %s %s' % (loc, sname)
            V = subprocess.call(cmd, shell=True)
        else:
            log.warning('Could not copy %s because does not exist' % loc)
        if op.exists(oloc):
            cmd = 'cp %s %s' % (oloc, soname)
            V = subprocess.call(cmd, shell=True)
        else:
            log.warning('Could not copy %s because does not exist' % oloc)
        
            
        