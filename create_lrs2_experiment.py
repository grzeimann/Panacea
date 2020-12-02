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

call = 'python3 /work/03730/gregz/maverick/Panacea/lrs2_experiment.py  %s -d %s -c %s'

com_call = 'python3 /work/03730/gregz/maverick/Panacea/combine_lrs2_experiment.py  %s %s'

args.log = setup_logging('lrs2_experiment')

filenames = sorted(glob.glob(op.join(args.directory, 'm*exp01*uv.fits')))

#da = bname.split('_')[1]
obj, ra, dec, ifuslot = ([], [], [], [])

def get_standards(date):
    standards = []
    while len(standards) == 0:
        standards = sorted(glob.glob(op.join(args.standirectory, 'm*%s*uv.fits' % date)))
        datet = datetime.datetime(int(date[:4]), int(date[4:6]), int(date[6:]))
        datet = datetime.timedelta(days=-1) + datet
        date = '%04d%02d%02d' % (datet.year, datet.month, datet.day)
    return standards

channels = ['uv', 'orange', 'red', 'farred']
make_calls = []
for filename in filenames:
    allfilenames = sorted(glob.glob(filename.replace('exp01', 'exp*')))
    calls = []
    flag = True
    for fn in allfilenames:
        f = fits.open(fn)
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
            flag = False
            continue
        if args.object is not None:
            try: 
                st = n.split(n[-6:])[0]
            except:
                continue
            if args.object.lower() not in st.lower():
                flag = False
                continue
        for chan in channels:
            calls.append(call % (op.basename(fn.replace('uv', chan)), args.directory, args.caldirectory))
    if flag:
        date = filename.split('_')[1]
        date = '20191118'
        standards = get_standards(date)
        name1 = '_'.join(op.basename(filename).split('_')[:3])  
        name2 = '_'.join(op.basename(standards[0]).split('_')[:3])
        name1 = name1.replace('multi', 'spectrum')
        name2 = name2.replace('multi', 'spectrum')
        for stan in standards:
            for chan in channels:
                calls.append(call % (op.basename(stan.replace('uv', chan)),
                                     args.standirectory, args.caldirectory))
        calls.append(com_call % (name1, name2))
        make_calls.append('; '.join(calls))

N = int(np.ceil(len(make_calls) / 20.))
chunks = np.array_split(make_calls, N)
for j, chunk in enumerate(chunks):
    with open('%s_%i' % (args.outname, j+1), 'w') as out_file:
        for prog in chunk:
            out_file.write(prog + '\n')