#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  8 15:37:24 2019

@author: gregz
"""
import glob
from astropy.io import fits
import numpy as np
import os.path as op

path = '/work/03946/hetdex/maverick/LRS2/*/multi_201902*uv.fits'
files = glob.glob(path)

names = ['uv', 'orange', 'red', 'farred']
s, w = ({}, {})

def get_fits(fn, name):
    F = fits.open(fn.replace('uv', name))
    return F[1].data * 1.

for name in names:
    s[name] = []
fn = files[0]
for name in names:
    F = fits.open(fn.replace('uv', name))
    w[name] = F[6].data[0]
for fn in files:
    F = fits.open(fn)
    flag = True
    for name in names:
        flag *= op.exists(fn.replace('uv', name))
    if F[0].header['EXPTIME'] > 120. and flag:
        for name in names:
            data = get_fits(fn, name)
            sky = np.median(data, axis=0)
            s[name].append(sky)

for name in names:
    fits.PrimaryHDU(np.vstack([w[name], s[name]])).writeto('sky_%s.fits' %
                   name, overwrite=True)