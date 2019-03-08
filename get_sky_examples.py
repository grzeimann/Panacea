#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  8 15:37:24 2019

@author: gregz
"""
import glob
from astropy.io import fits
import numpy as np

path = '/work/03946/hetdex/maverick/LRS2/[UT,PSU,M]*/multi_201902*uv.fits'
files = glob.glob(path)

names = ['uv', 'orange', 'red', 'farred']
s, w = ({}, {})

for name in names:
    s[name] = []
fn = files[0]
for name in names:
    F = fits.open(fn.replace('uv', name))
    w[name] = F[6].data[0]
for fn in files:
    F = fits.open(fn)
    if F[0].header['EXPTIME'] > 120.:
        for name in names:
            F = fits.open(fn.replace('uv', name))
            sky = np.median(F[1].data, axis=0)
            s[name].append(sky)

for name in names:
    fits.PrimaryHDU(np.vstack([w[name], s[name]])).writeto('sky_%s.fits' %
                   name, overwrite=True)