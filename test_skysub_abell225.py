# -*- coding: utf-8 -*-
"""
Created on Tue Jan  8 06:44:36 2019

@author: gregz
"""

import numpy as np
import warnings
from scipy.interpolate import LSQBivariateSpline
from scipy.signal import medfilt
from astropy.io import fits
from astropy.table import Table
def get_selection(array1, array2):
    m1 = medfilt(array1, 5)
    m2 = medfilt(array2, 5)
    y1 = np.abs(array1 - m1)
    y2 = np.abs(array2 - m2)
    mad1 = np.nanmedian(y1)
    mad2 = np.nanmedian(y2)
    return (y1 < (5 * mad1)) * (y2 < (5 * mad2))

def solve_system(sci_list, sky_list, x, y, xoff, yoff, sci_image):
    norm1 = np.zeros((sci_list.shape[1],))
    norm2 = np.zeros((sci_list.shape[1],))
    newsci = sci_list * 0.
    newsky = sky_list * 0.
    C = np.zeros((len(x), 2))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        I = LSQBivariateSpline(x, y, sci_image, np.linspace(-6.0, 6.0, 27),
                               np.linspace(-3.5, 3.5, 15))
    for j in np.arange(sci_list.shape[1]):
        if sci_image.ndim == 1:
            xnew = x - xoff[j]
            ynew = y - yoff[j]
            C[:, 0] = I(xnew, ynew, grid=False)
        else:
            C[:, 0] = sci_image[:, j] * 1.
        sel = (np.isfinite(sci_list[:, j]) *
               np.isfinite(sky_list[:, j]))
        C[:, 1] = sky_list[:, j] * 1.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            sol = np.linalg.lstsq(C[sel], sci_list[sel, j])[0]
        norm1[j] = sol[0]
        norm2[j] = sol[1]
        newsci[:, j] = C[:, 0] * sol[0]
        newsky[:, j] = C[:, 1] * sol[1]
    return newsky, newsci, norm1, norm2

F = fits.open('FeRpses20171226T043454.3_056_sci_R.fits')
fdata = F[0].data
G = fits.open('FeRpses20171226T044832.4_056_sci_R.fits')
gdata = G[0].data
wave = np.zeros((F[0].data.shape[1],))
wave[0] = F[0].header['CRVAL1']
for i in np.arange(1, F[0].data.shape[1]):
    wave[i] = wave[i-1] * (1 + 30. / 299792.458)
wave_0 = np.mean(wave[1955:2412])
darfile = '/Users/gregz/cure/panacea/lrs2_config/dar_BR.dat'
T = Table.read(darfile, format='ascii.fixed_width_two_line')
xoff = (np.interp(wave, T['wave'], T['x_0']) -
        np.interp(wave_0, T['wave'], T['x_0']))
yoff = (np.interp(wave, T['wave'], T['y_0']) -
        np.interp(wave_0, T['wave'], T['y_0']))
X = np.loadtxt('/Users/gregz/cure/LRS2_reduction/lrs2_config/mapping_files/LRS2_B_OR_mapping.txt', skiprows=11, usecols=(1, 2))
sci_image = np.median((fdata - gdata)[:, 1955:2412], axis=1)
sky, dummy1, dummy2, dummy3 = solve_system(fdata, gdata, X[:, 0], X[:, 1],
                                           xoff, yoff, sci_image)
f1 = fits.PrimaryHDU(fdata)
f2 = fits.ImageHDU(gdata)
f3 = fits.ImageHDU(sky)
f4 = fits.ImageHDU(fdata - sky)
fits.HDUList([f1, f2, f3, f4]).writeto('test.fits', overwrite=True)
