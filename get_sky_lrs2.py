#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun  4 14:48:58 2022

@author: gregz
"""

import numpy as np
from astropy.io import fits
import glob
import os.path as op


filenames = np.array(sorted(glob.glob('multi*uv.fits')))
dates = [op.basename(f).split('_')[1] for f in filenames]
dates = np.array(dates)

wave = np.arange(3650, 10500, 0.7)
skies = []
for date in np.unique(dates):
    inds = np.where(date == dates)[0]
    mnames = filenames[inds]
    for filename in mnames:
        f = fits.open(filename)
        name = f[0].header['OBJECT']
        millum = f[0].header['MILLUM']
        throughp = f[0].header['THROUGHP']
        if millum == 51e4:
            continue
        if (throughp < 0.1) + (throughp == 1.):
            continue
        slot = name.split('_')[-2]
        if slot == '056':
            f = fits.open(filename.replace('uv', 'red'))
            g = fits.open(filename.replace('uv', 'farred'))
            skyred = f[1].data[125]
            skyfarred = g[1].data[125]
            wavered = f[6].data[0]
            wavefarred = g[6].data[0]
            skyre = np.interp(wave, wavered, skyred, 
                              left=np.nan, right=np.nan)
            skyfr = np.interp(wave, wavefarred, skyfarred, 
                              left=np.nan, right=np.nan)
        if slot == '066':
            f = fits.open(filename.replace('uv', 'uv'))
            g = fits.open(filename.replace('uv', 'orange'))
            skyuv = f[1].data[125]
            skyorange = g[1].data[125]
            waveuv = f[6].data[0]
            waveorange = g[6].data[0]
            skyuv = np.interp(wave, waveuv, skyuv, 
                              left=np.nan, right=np.nan)
            skyor = np.interp(wave, waveorange, skyorange, 
                              left=np.nan, right=np.nan)
    wsel = (wave>6530) * (wave < 6880)
    norm = np.nanmedian(skyre[wsel] / skyor[wsel])
    print('Norm for %s %s: %0.2f' % (filename, name, norm))
    skyuv *= norm
    skyor *= norm
    sky = np.nanmean([skyuv, skyor, skyre, skyfr], axis=0)
    skies.append(sky)
skies = np.array(skies)
fits.PrimaryHDU(skies).writeto('skyfile.fits', overwrite=True)