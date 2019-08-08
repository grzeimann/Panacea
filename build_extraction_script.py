#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  8 10:48:33 2019

@author: gregz
"""

import argparse as ap
import numpy as np

from astropy.io import fits

parser = ap.ArgumentParser(add_help=True)

parser.add_argument("filename",
                    help='''Name of fits file with Bin Table''',
                    type=str)

parser.add_argument("chunks",
                    help='''Break file into this number of chunks''',
                    type=int)

args = parser.parse_args(args=None)


fitsfile = fits.open(args.filename)
bintable = fitsfile[1].data

inds = np.argsort(bintable['ra'])

s = []
ind_chunks = np.array_split(inds, args.chunks)
for i, ind in enumerate(ind_chunks):
    List = [fits.PrimaryHDU(), fits.BinTableHDU(bintable[ind])]
    F = fits.HDUList(List)
    F.writeto('data_chunk_%03d.fits' % (i + 1), overwrite=True)
    s.append('python /work/03730/gregz/maverick/Panacea/get_extractions.py data_chunk_%03d.fits '
             'gaia_hetdex_spectra_%03d.fits' % (i + 1, i + 1))
f = open('extraction_script', 'w')
f.write('\n'.join(s) + '\n')