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

parser.add_argument("caldirectory",
                    help='''cal directory for reductions''', type=str)

parser.add_argument("outname",
                    help='''Name of output file''', type=str)

parser.add_argument("--atfile",
                    help='''Name of at file for re-reduction''', type=str,
                    default=None)

parser.add_argument("--object",
                    help='''Name of Object''', type=str,
                    default=None)

parser.add_argument("-sd", "--sep_date",
                    help='''Separate Dates''',
                    action="count", default=0)
args = parser.parse_args(args=None)

atcall = 'echo "source ~hetdex/.bashrc_greg; runlrs2wranglergeneral %s %s" | at %s'

args.log = setup_logging('advance_cube_creation')

filenames = sorted(glob.glob(op.join(args.directory, 'm*uv.fits')))

#da = bname.split('_')[1]
obj, ra, dec, ifuslot = ([], [], [], [])
keep_files = []
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
    st = n.split(n[-6:])[0]
    try:
        ifuslot.append(n.split('_')[-2])
    except:
        continue
    if args.object is not None:
        if args.object.lower() not in st.lower():
            continue
    ra.append(r)
    dec.append(d)
    obj.append(st)
    keep_files.append(filename)
uobj = np.unique(obj)
calls, atcalls = ([], [])
now = datetime.datetime.now()
for o in uobj:
    inds = [i for i, ob in enumerate(obj) if o == ob]
    blue, red, sky = ([], [], [])
    raspl = ra[inds[0]].split(':')
    decspl = dec[inds[0]].split(':')
    rah = raspl[0] + 'h' + raspl[1] + 'm' + raspl[2] + 's'
    dech = decspl[0] + 'd' + decspl[1] + 'm' + decspl[2] + 's'
    if '-' in decspl[0]:
        dech = ' ' + dech
    dates = []
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
        date = bname.split('_')[1]
        dates.append(date)
    rdates = [name.split('_')[1] for name in red]
    bdates = [name.split('_')[1] for name in blue]
    sdates = [name.split('_')[1] for name in sky]

    udates = np.unique(np.hstack([rdates, bdates, sdates]))
    if args.sep_date:
        for j, udate in enumerate(udates):
            di = [i for i, d in enumerate(dates) if d == udate]
            B = [blue[i] for i, d in enumerate(bdates) if d == udate]
            R = [red[i] for i, d in enumerate(rdates) if d == udate]
            S = [sky[i] for i, d in enumerate(sdates) if d == udate]
            B = ','.join(B)
            R = ','.join(R)
            S = ','.join(S)
            call = ('python Panacea/advanced_cube_creation.py %s "' + B + '" "' +
                    R + '" "' + S + '" "' + rah + '" "' + dech + '" ' + 
                    "-d %s -c %s -dw 0.7") % (o+'_%s' % udate , args.directory, args.caldirectory)
            calls.append(call)
            now = now + datetime.timedelta(seconds=240)
            tim = now.strftime('%H:%M %B %d')
            atcalls.append(atcall % (udate, o, tim))
    else:
        B = ','.join(blue)
        R = ','.join(red)
        S = ','.join(sky)
        call = ('python Panacea/advanced_cube_creation.py %s "' + B + '" "' +
                R + '" "' + S + '" "' + rah + '" "' + dech + '" ' + 
                "-d %s -c %s -dw 0.7") % (o, args.directory, args.caldirectory)
        calls.append(call)
        for udate in udates:
            now = now + datetime.timedelta(seconds=120)
            tim = now.strftime('%H:%M %B %d')
            atcalls.append(atcall % (udate, o, tim))
with open(args.outname, 'w') as out_file:
    for call in calls:
        out_file.write(call + '\n')
if args.atfile is not None:
    with open(args.atfile, 'w') as out_file:
        for call in atcalls:
            out_file.write(call + '\n')        


            
    
