#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  7 09:31:56 2019

@author: gregz
"""

import glob
import datetime
import os.path as op
import numpy as np
import argparse as ap
from astropy.io import fits
from input_utils import setup_logging


parser = ap.ArgumentParser(add_help=True)

parser.add_argument("directory",
                    help='''base directory for reductions''', type=str)

parser.add_argument("outname",
                    help='''Name of output file''', type=str)

parser.add_argument("--object",
                    help='''Name of Object''', type=str,
                    default=None)

parser.add_argument("--caldirectory", 
                    default='/work/03946/hetdex/maverick/LRS2/CALS',
                    help='''cal directory for reductions''', type=str)

parser.add_argument("--standirectory", 
                    default='/work/03946/hetdex/maverick/LRS2/STANDARDS',
                    help='''Standards directory for reductions''', type=str)

args = parser.parse_args(args=None)

call = 'python Panacea/lrs2_experiment.py  %s -d %s -c %s'

args.log = setup_logging('lrs2_experiment')

filenames = sorted(glob.glob(op.join(args.directory, 'm*uv.fits')))

#da = bname.split('_')[1]
obj, ra, dec, ifuslot = ([], [], [], [])

channels = ['uv', 'orange']
make_calls = []
for filename in filenames:
    f = fits.open(filename)
    try:
        n = f[0].header['OBJECT']
    except:
        continue
    try:
        r = f[0].header['QRA']
    except:
        r = '12:00:00'
    try:
        d = f[0].header['QDEC']
    except:
        d = '+00:00:00'
    try:
        ifuslot.append(n.split('_')[-2])
    except:
        continue
    if args.object is not None:
        try: 
            st = n.split(n[-6:])[0]
        except:
            continue
        if args.object.lower() not in st.lower():
            continue
    date = filename.split('_')[1]
    standards = sorted(glob.glob(op.join(args.standirectory, 'm*%s*uv.fits' % date)))
    calls = []
    for chan in channels:
        calls.append(call % (op.basename(filename.replace('uv', chan)), args.directory, args.caldirectory))
    for stan in standards:
        for chan in channels:
            calls.append(call % (op.basename(stan.replace('uv', chan)),
                                 args.standirectory, args.caldirectory))
    make_calls.append('; '.join(calls))

N = int(np.ceil(len(make_calls) / 20.))
chunks = np.array_split(make_calls, N)
for j, chunk in enumerate(chunks):
    with open('%s_%i' % (args.outname, j+1), 'w') as out_file:
        for prog in chunk:
            out_file.write(prog + '\n')