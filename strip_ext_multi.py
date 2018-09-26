# -*- coding: utf-8 -*-
"""
Created on Wed Sep 26 09:15:18 2018

@author: gregz
"""

import numpy as np
from astropy.io import fits

filenames = np.loadtxt('/work/03730/gregz/maverick/multi_orange_list.dat')
ext = 'skysub'

fitslist = []
for i, filename in enumerate(filenames):
    F = fits.open(filename)
    if i == 0:
        f = fits.PrimaryHDU(F[ext].data)
    else:
        f = fits.ImageHDU(F[ext].data)
    fitslist.append(f)

fits.HDUList(fitslist).writeto('all_multi.fits', overwrite=True)