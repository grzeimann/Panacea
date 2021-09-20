#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 20 10:46:25 2021

@author: gregz
"""

from astropy.io import fits
import numpy as np
import glob
import os.path as op
from input_utils import setup_logging

log = setup_logging('adr')
basedir = '/work/03946/hetdex/maverick/LRS2/STANDARDS'

filenames = glob.glob(op.join(basedir, 'multi_2021*uv.fits'))

N = 11
M = len(filenames)
k=0
XC, YC, WC = (np.zeros((M, 2*N)), np.zeros((M, 2*N)), np.zeros((2*N,)))

for filename in filenames:
    log.info('Working on %s' % filename)
    j = 0
    for name in ['uv', 'orange']:
        f = fits.open(filename)
        wave = f[6].data[0]
        x = f[5].data[:, 0]
        y = f[5].data[:, 1]
        skysub = f[2].data 
        chunks = [np.nanmedian(xo, axis=0) for xo in np.array_split(skysub, N,
                                                                    axis=0)]
        wc = [np.nanmedian(xo) for xo in np.array_split(wave, N)]
        xc = np.zeros((N,))
        yc = np.zeros((N,))
        wc = np.array(wc)
        for i, chunk in enumerate(chunks):
            ind = np.nanargmax(chunk)
            d = np.sqrt((x-x[ind])**2 + (y-y[ind])**2)
            sel = (d < 2.) * np.isfinite(chunk)
            xc[i] = np.sum(x[sel]*chunk[sel]) / np.sum(chunk[sel])
            yc[i] = np.sum(y[sel]*chunk[sel]) / np.sum(chunk[sel])
        XC[k, j:j+N] = xc
        YC[k, j:j+N] = yc
        WC[j:j+N] = wc
        j += N
    
    xc = np.interp(5500., WC, XC[k])
    yc = np.interp(5500., WC, YC[k])
    XC[k] = XC[k] - xc
    YC[k] = YC[k] - yc
    log.info('Centroid: %0.2f %0.2f' % (xc, yc))
    k += 1

fits.HDUList([fits.PrimaryHDU(XC), fits.ImageHDU(YC), fits.ImageHDU(WC)]).writeto('test.fits', overwrite=True)
    
        