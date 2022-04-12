#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 15 09:49:22 2022

@author: gregz
"""

import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import os.path as op
from astropy.convolution import convolve, Gaussian2DKernel, Gaussian1DKernel
from astropy.table import Table
from scipy.interpolate import interp1d, griddata
from astropy.modeling.models import Gaussian2D
from astropy.modeling.fitting import LevMarLSQFitter

from sklearn.decomposition import PCA
import seaborn as sns



_lines = [3726.1, 3729.1, 3889., 4101.76, 4340.5,  4363.2, 4471., 4861.3, 
          4958.9, 5006.8, 5875.7, 6300.3, 6312.1, 6548., 6562.8, 6583.4, 
          6678., 6716.5, 6730.8]

# Plotting style to mended depending on necessity
sns.set_context('poster')
sns.set_style('ticks')

plt.rcParams["font.family"] = "Times New Roman"

def get_continuum(y, sel, bins=25):
    yz = y * 1.
    yz[sel] = np.nan
    x = np.array(np.arange(len(y)), dtype=float)
    xc = np.array([np.nanmean(xi) for xi in np.array_split(x, bins)])
    yc = np.array([np.nanmedian(xi) for xi in np.array_split(y, bins)])
    sel = np.isfinite(yc)
    I = interp1d(xc[sel], yc[sel], kind='linear', bounds_error=False,
                 fill_value='extrapolate')
    return I(x)

def pca_fit(H, data, sel):
    sel = sel * np.isfinite(data)
    sol = np.linalg.lstsq(H.T[sel], data[sel])[0]
    res = np.dot(H.T, sol)
    return res

def get_sky_pixels(sky, init_sel):
    for j in np.arange(2):
        cont = get_continuum(sky, init_sel, bins=15)
        mad = np.nanmedian(np.abs(sky-cont - np.nanmedian(sky-cont)))
        mask = sky-cont > 5. * mad
        for i in np.arange(1, 6):
            mask[i:] += mask[:-i]
            mask[:-i] += mask[i:]
        init_sel += mask
    return mask

def make_new_cube(cube, xgrid, ygrid, xc, yc, size=7):
    S = np.zeros((xgrid.shape[0]*xgrid.shape[1], 2))
    S[:, 0] = xgrid.ravel()-xc
    S[:, 1] = ygrid.ravel()-yc
    xg, yg = np.meshgrid(np.linspace(-size, size, size*2*4+1),
                         np.linspace(-size, size, size*2*4+1))
    newcube = np.ones((cube.shape[0], size*2*4+1, size*2*4+1))*np.nan
    for i in np.arange(cube.shape[0]):
        sel = np.isfinite(cube[i].ravel())
        if sel.sum()>20:
            newcube[i] = griddata(S[sel], cube[i].ravel()[sel], (xg, yg),
                                  method='linear')
    return newcube, xg, yg

def get_data_and_grid(filename, xc, yc, sky=False):
    k = fits.open(filename)
    x = np.arange(k[0].header['NAXIS1'])*k[0].header['CDELT1'] + k[0].header['CRVAL1']
    y = np.arange(k[0].header['NAXIS2'])*k[0].header['CDELT2'] + k[0].header['CRVAL2']
    w = np.arange(k[0].header['NAXIS3'])*k[0].header['CDELT3'] + k[0].header['CRVAL3']
    xgrid, ygrid = np.meshgrid(x, y)
    if sky:
        k2 = fits.open(filename.replace('_cube.fits', '_sky_cube.fits'))
        data = k[0].data + k2[0].data
    else:
        data = k[0].data
    k2 = fits.open(filename.replace('_cube.fits', '_error_cube.fits'))
    datae = k2[0].data
    return w, xgrid, ygrid, data, datae, k[0].header

def write_cube(wave, xgrid, ygrid, Dcube, outname, he):
    '''
    Write data cube to fits file
    
    Parameters
    ----------
    wave : 1d numpy array
        Wavelength for data cube
    xgrid : 2d numpy array
        x-coordinates for data cube
    ygrid : 2d numpy array
        y-coordinates for data cube
    Dcube : 3d numpy array
        Data cube, corrected for ADR
    outname : str
        Name of the outputted fits file
    he : object
        hdu header object to carry original header information
    '''
    hdu = fits.PrimaryHDU(np.array(Dcube, dtype='float32'))
    hdu.header['CRVAL1'] = xgrid[0, 0]
    hdu.header['CRVAL2'] = ygrid[0, 0]
    hdu.header['CRVAL3'] = wave[0]
    hdu.header['CRPIX1'] = 1
    hdu.header['CRPIX2'] = 1
    hdu.header['CRPIX3'] = 1
    hdu.header['CTYPE1'] = 'pixel'
    hdu.header['CTYPE2'] = 'pixel'
    hdu.header['CTYPE3'] = 'pixel'
    hdu.header['CDELT1'] = xgrid[0, 1] - xgrid[0, 0]
    hdu.header['CDELT2'] = ygrid[1, 0] - ygrid[0, 0]
    hdu.header['CDELT3'] = wave[1] - wave[0]
    for key in he.keys():
        if key in hdu.header:
            continue
        if ('CCDSEC' in key) or ('DATASEC' in key):
            continue
        if ('BSCALE' in key) or ('BZERO' in key):
            continue
        hdu.header[key] = he[key]
    hdu.writeto(outname, overwrite=True)

base = '/Users/gregz/cure/Remedy/for_ashley/'
bfile = op.join(base, 'SN2022erw_20220314_LRS2B_cube.fits')
rfile = op.join(base, '2022gnp_20220411_LRS2R_cube.fits')
objname = '2022gnp_20220411'
#posb = (-2.00, -2.5)
posr = (0.25, -0.75)
#posr = (-1.6, -0.8)
#posbh = (-3.75, -0.25)
posrh = (-3.75, -0.25)
#posrh = (-3.35, 1.45)
radius = 2.0
redshift = 0.01086
allwave = np.arange(3650., 10500.7, 0.7)
allwave = np.arange(6450., 10500.7, 0.7)

spec, SKY, BACK, error = ([], [], [], [])
for filename, pos, posh in zip([rfile], [posr], [posrh]):
    xc, yc = pos
    xc2, yc2 = posh
    wave, xgrid, ygrid, data, datae, header = get_data_and_grid(filename, xc, yc,
                                                                sky=True)
    orign = np.zeros(wave.shape, dtype=bool)
    for line in _lines:
        orign[np.abs(wave-line*(1+redshift))<10.] = True
    uvmask = np.abs(wave-3736.0) < 1.6
    data[uvmask] = np.nan
    datae[uvmask] = np.nan
    d = np.sqrt((xgrid-xc)**2 + (ygrid-yc)**2)
    d2 = np.sqrt((xgrid-xc2)**2 + (ygrid-yc2)**2)

    skysel = (d > 4.) * (xgrid > -1.) * (np.isfinite(data).sum(axis=0) > 3000.)
    detwave = (1+redshift) * 6563.
    wsel = np.abs(wave-detwave) < 5.
    image = np.nanmean(data[wsel], axis=0)
    skysel = (np.isfinite(data).sum(axis=0) > 3800.)* (image < np.nanpercentile(image, 25))
    plt.figure()
    plt.imshow(image, origin='lower')
    plt.imshow(skysel, origin='lower', alpha=0.5)
    plt.show()
    sky = np.nanmedian(data[:, skysel], axis=1)
    
    sky[wave<4650.] = sky[wave<4650]*0.92
    skydata = data * 1.
    skydata[:] = sky[:, np.newaxis, np.newaxis]
    skydata[np.isnan(data)] = np.nan
    data = data - sky[:, np.newaxis, np.newaxis]
    # New smoothed/interpolated spectra
    for i in np.arange(data.shape[1]):
        for j in np.arange(data.shape[2]):
            sel = np.isnan(data[:, i, j])
            if (~sel).sum() > 100.:
                data[:, i, j] = np.interp(wave, wave[~sel], data[~sel, i, j], left=np.nan, 
                                 right=np.nan)
                datae[:, i, j] = np.interp(wave, wave[~sel], datae[~sel, i, j], left=np.nan, 
                                 right=np.nan)
                datae[sel, i, j] = datae[sel, i, j]*1.5
    
    cont_cube = data * np.nan
    for i in np.arange(data.shape[1]):
        for j in np.arange(data.shape[2]):
            if np.isfinite(data[:,i,j]).sum() > 1000.:
                fsel = np.isnan(data[:, i, j])
                cont_cube[:, i, j] = get_continuum(data[:, i, j], orign, bins=100)
                cont_cube[fsel, i, j] = np.nan
    
    cont_sub = data - cont_cube
    cont_sub[-1] = 0.0
    objsel = np.isfinite(cont_sub).sum(axis=0) == len(cont_sub)
    objsel = objsel  * (d > radius)
    # Fit PCA Model
    pcawave = np.ones(wave.shape, dtype=bool)
    pca = PCA(n_components=5).fit(cont_sub[pcawave][:, objsel].T)
    Hk = pca.components_
    
    em_model = data * 0.
    # Fit residuals for PCA eigenvalues and subtract model
    for i in np.arange(data.shape[1]):
        for j in np.arange(data.shape[2]):
            yp = cont_sub[:, i, j]
            res = pca_fit(Hk, yp[pcawave], 
                          np.ones((pcawave.sum(),), dtype=bool))
            ycopy = np.nan * yp
            ycopy[pcawave] = res
            ycopy[np.isnan(yp)] = np.nan
            em_model[:, i, j] = ycopy
    
    background = cont_cube * 0.
    pixsize = xgrid[0, 1] - xgrid[0, 0]
    for i in np.arange(cont_cube.shape[0]):
        image = cont_cube[i] * 1.
        nanmask = np.isnan(image)
        mask = d <= radius
        image[mask] = np.nan
        smooth = convolve(image, Gaussian2DKernel(radius*1.0/2.35/pixsize), 
                          boundary='fill', fill_value=0.0)
        smooth[nanmask] = np.nan
        background[i] = smooth
        
    outname = op.join(base, '%s_sky.fits' % op.basename(filename).split('_cube')[0])
    write_cube(wave, xgrid, ygrid, skydata, outname, header)
    outname = op.join(base, '%s_pca.fits' % op.basename(filename).split('_cube')[0])
    write_cube(wave, xgrid, ygrid, em_model, outname, header)
    outname = op.join(base, '%s_smooth.fits' % op.basename(filename).split('_cube')[0])
    write_cube(wave, xgrid, ygrid, background, outname, header)
    outname = op.join(base, '%s_sub.fits' % op.basename(filename).split('_cube')[0])
    sub = data-background-em_model
    back = background + em_model
    write_cube(wave, xgrid, ygrid, sub, outname, header)
    fitter = LevMarLSQFitter()
    image = np.nanmedian(sub[100:300], axis=0) * 1e18
    sel = (d < 3.) * np.isfinite(image)
    G = Gaussian2D(amplitude=np.nanmax(image), x_mean=xc, y_mean=yc)
    fit = fitter(G, xgrid[sel], ygrid[sel], image[sel])
    for n in fit.param_names:
        print('%s: %0.2f' % (n, getattr(fit, n).value))
    W = fit(xgrid, ygrid)
    W = W / W.sum()
    ispec = np.zeros((sub.shape[0],))
    iback = np.zeros((sub.shape[0],))
    isky = np.zeros((sub.shape[0],))
    ierror = np.zeros((sub.shape[0],))
    for i in np.arange(sub.shape[0]):
        dsel = (d<=radius) * np.isfinite(sub[i])
        d2sel = (d2<=radius) * np.isfinite(sub[i])

        cor = W[dsel].sum()
        ispec[i] = (np.nansum(sub[i, dsel] * W[dsel]) /
                       np.sum(W[dsel]**2)) / cor
        isky[i] = np.nansum(skydata[i, dsel])
        iback[i] = np.nansum(back[i, d2sel])
        ierror[i] = (np.sqrt(np.nansum(datae[i, dsel]**2 *
                                      W[dsel]**2)) /
                     np.sum(W[dsel]**2)) / cor
    spec.append(np.interp(allwave, wave, ispec, left=np.nan, right=np.nan))
    error.append(np.interp(allwave, wave, ierror, left=np.nan, right=np.nan))
    SKY.append(np.interp(allwave, wave, isky, left=np.nan, right=np.nan))
    BACK.append(np.interp(allwave, wave, iback, left=np.nan, right=np.nan))

# overlap = np.isfinite(spec[0]) * np.isfinite(spec[1])
# norm = np.nanmedian(spec[0][overlap] / spec[1][overlap])
# normback = np.nanmedian(BACK[0][overlap] / BACK[1][overlap])

# normsky = np.nanmedian(SKY[0][overlap] / SKY[1][overlap])

# G = Gaussian1DKernel(2.5)
# spec[1] = convolve(spec[1], G)  * norm
# error[1] = error[1] * norm
# SKY[1] = convolve(SKY[1], G) * normsky
# BACK[1] = convolve(BACK[1], G) * normback

# avg = np.nanmean(spec, axis=0)
# sky = np.nanmean(SKY, axis=0)
# avgerror = np.nanmean(error, axis=0)
# back = np.nanmean(BACK, axis=0)
avg = spec[0]
sky = SKY[0]
back = BACK[0]
avgerror = error[0]

T = Table([allwave, avg, avgerror, back, sky], names=['wavelength', 'f_lam', 'e_lam', 'b_lam', 's_lam'])
T.write(op.join(base, '%s_spec.dat' % objname), format='ascii.fixed_width_two_line')

wran = [[3650, 10500], [3650, 3900], [3900, 5200], [6600, 6950],
        [8000, 10000]]
wran = [[6450, 10500], [6600, 6950], [7300, 7800], [7900, 8800],
        [8800, 9500]]
fig, ax = plt.subplots(5, 1, figsize=(20, 10))
ax[0].set_position([0.1, 0.53, 0.86, 0.42])
ax[1].set_position([0.1, 0.08, 0.19, 0.42])
ax[2].set_position([0.323, 0.08, 0.19, 0.42])
ax[3].set_position([0.546, 0.08, 0.19, 0.42])
ax[4].set_position([0.77, 0.08, 0.19, 0.42])
for i, wr in enumerate(wran):
    wsel = (allwave>wr[0]) * (allwave<wr[1])
    ax[i].plot(allwave[wsel], sky[wsel]*1e17, color='olive', lw=0.5, alpha=0.5, label='Sky Model')
    ax[i].plot(allwave[wsel], back[wsel]*1e17, color='salmon', lw=0.5, label='Background Galaxy + HII region + sky residuals')
    ax[i].plot(allwave[wsel], avg[wsel]*1e17, color='k', lw=0.5, label='Target')
    f_ax10 = ax[i]
    f_ax10.tick_params(axis='both', which='both', direction='in')
    f_ax10.tick_params(axis='y', which='both', left=True, right=True)
    f_ax10.tick_params(axis='x', which='both', bottom=True, top=True)
    f_ax10.tick_params(axis='both', which='major', length=8, width=2)
    f_ax10.tick_params(axis='both', which='minor', length=5, width=1)
    f_ax10.minorticks_on()
    ax[i].set_xlim([wr[0], wr[1]])
    ax[i].set_ylim([-5, 75])
ax[0].legend()
ax[0].set_ylabel(r'F$_{\lambda}$ (10$^{-17}$erg s$^{-1}$ cm$^{-1}$ ${\AA}^{-1}$)', labelpad=10, fontsize=22)
plt.savefig(op.join(base, '%s_plot.png' % objname), dpi=150)
