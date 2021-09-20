#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 20 10:46:25 2021

@author: gregz
"""

from astropy.io import fits
from astropy.modeling.models import Gaussian2D
from astropy.modeling.fitting import LevMarLSQFitter
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

fitter = LevMarLSQFitter()

G = Gaussian2D()

for filename in filenames:
    log.info('Working on %s' % filename)
    j = 0
    for name in ['uv', 'orange']:
        f = fits.open(filename.replace('uv', name))
        wave = f[6].data[0]
        x = f[5].data[:, 0]
        y = f[5].data[:, 1]
        skysub = f[2].data 
        chunks = [np.nanmedian(xo, axis=1) for xo in np.array_split(skysub, N,
                                                                    axis=1)]
        wc = [np.nanmedian(xo) for xo in np.array_split(wave, N)]
        xc = np.zeros((N,))
        yc = np.zeros((N,))
        wc = np.array(wc)
        for i, chunk in enumerate(chunks):
            ind = np.nanargmax(chunk)
            d = np.sqrt((x-x[ind])**2 + (y-y[ind])**2)
            xc[i] = x[ind]
            yc[i] = y[ind]
            G.amplitude.value = np.nanmax(chunk)
            G.x_mean.value = xc[i]
            G.y_mean.value = yc[i]
            sel = (d < 3.0) * np.isfinite(chunk)
            fit = fitter(G, x[sel], y[sel], chunk[sel])
            xc[i] = fit.x_mean.value * 1.
            yc[i] = fit.y_mean.value * 1.
        XC[k, j:j+N] = xc
        YC[k, j:j+N] = yc
        WC[j:j+N] = wc
        j += N
    
    XC[k, WC<4650.] += 0.20
    YC[k, WC<4650.] -= 0.24
    sel = np.isfinite(XC[k]) * np.isfinite(YC[k])
    p0 = np.polyfit(WC[sel], XC[k][sel], 3)
    p1 = np.polyfit(WC[sel], YC[k][sel], 3)
    xc = np.polyval(p0, 5500)
    yc = np.polyval(p1, 5500)
    XC[k] = XC[k] - xc
    YC[k] = YC[k] - yc
    if np.sqrt(xc**2 + yc**2)>2.5:
        XC[k] = np.nan
        YC[k] = np.nan
    log.info('Centroid: %0.2f %0.2f' % (xc, yc))
    k += 1

fits.HDUList([fits.PrimaryHDU(XC), fits.ImageHDU(YC), fits.ImageHDU(WC)]).writeto('test.fits', overwrite=True)
    
        