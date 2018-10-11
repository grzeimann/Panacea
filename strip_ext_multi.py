# -*- coding: utf-8 -*-
"""
Created on Wed Sep 26 09:15:18 2018

@author: gregz
"""

import numpy as np
from astropy.io import fits

filenames = [line.rstrip('\n').split() for line in open('/work/03730/gregz/maverick/test_2.dat', 'r')]

ext = 'spectrum'

fitslist = []
for filename, amp in zip(filenames, ['LL', 'LU', 'RU', 'RL']):
    F = fits.open(filename[0])
    if amp == 'LL':
        f = fits.PrimaryHDU(F[ext].data)
        f.header['EXTNAME']=amp
    else:
        f = fits.ImageHDU(F[ext].data)
        f.header['EXTNAME']=amp
    fitslist.append(f)

fits.HDUList(fitslist).writeto('sky_august2018.fits', overwrite=True)
