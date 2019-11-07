#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  7 09:31:56 2019

@author: gregz
"""

import glob

import os.path as op
import numpy as np
import argparse as ap
from astropy.io import fits
from input_utils import setup_logging


parser = ap.ArgumentParser(add_help=True)

parser.add_argument("directory",
                    help='''base directory for reductions''', type=str)

parser.add_argument("caldirectory",
                    help='''cal directory for reductions''', type=str)

parser.add_argument("outname",
                    help='''Name of output file''', type=str)
args = parser.parse_args(args=None)

args.log = setup_logging('advance_cube_creation')

filenames = sorted(glob.glob(op.join(args.directory, 'm*uv.fits')))

#da = bname.split('_')[1]
obj, ra, dec, ifuslot = ([], [], [], [])
keep_files = []
for filename in filenames:
    f = fits.open(filename)    
    n, r, d = (f[0].header['OBJECT'], f[0].header['QRA'], f[0].header['QDEC'])
    
    st = n.split(n[-6:])[0]
    try:
        ifuslot.append(n.split('_')[-2])
    except:
        continue
    ra.append(r)
    dec.append(d)
    obj.append(st)
    keep_files.append(filename)
uobj = np.unique(obj)
calls = []
for o in uobj:
    inds = [i for i, ob in enumerate(obj) if o == ob]
    blue, red, sky = ([], [], [])
    raspl = ra[inds[0]].split(':')
    decspl = dec[inds[0]].split(':')
    rah = raspl[0] + 'h' + raspl[1] + 'm' + raspl[2] + 's'
    dech = decspl[0] + 'd' + decspl[1] + 'm' + decspl[2] + 's'
    if '-' in decspl[0]:
        dech = ' ' + dech
    for ind in inds:
        filename = keep_files[ind]
        bname = op.basename(filename)
        bsky = bname.replace('uv', 'red')
        rname = bsky
        rsky = bname
        if ifuslot[ind] == '056':
            blue.append(bname)
            sky.append(bsky)
        if ifuslot[ind] == '066':
            red.append(rname)
            sky.append(rsky)
    blue = ','.join(blue)
    red = ','.join(red)
    sky = ','.join(sky)
    call = ('python %s Panacea/advanced_cube_creation.py "' + blue + '" "' +
            red + '" "' + sky + '" "' + rah + '" "' + dech + '" ' + 
            "-d %s -c %s") % (o, args.directory, args.caldirectory)
    calls.append(call)
    
with open(args.outname, 'w') as out_file:
    for call in calls:
        out_file.write(call + '\n')
        
        


            
    
