#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  3 08:27:47 2020

@author: gregz
"""

import argparse as ap
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
import os.path as op
import seaborn as sns
import warnings
from astropy.convolution import convolve, Gaussian2DKernel
from astropy.io import fits
from astropy.modeling.models import Gaussian1D, Moffat2D
from astropy.modeling.fitting import LevMarLSQFitter
from astropy.visualization.stretch import AsinhStretch, LinearStretch
from astropy.visualization import ImageNormalize
from astropy.table import Table
from math_utils import biweight
from matplotlib import rc
from input_utils import setup_logging
from scipy.interpolate import griddata, interp1d
warnings.simplefilter("ignore")

parser = ap.ArgumentParser(add_help=True)


parser.add_argument("filename",
                    help='''Science Star File Name''',
                    type=str)

parser.add_argument("standard_name",
                    help='''Standard Star Name''',
                    type=str)

parser.add_argument("-sd", "--standard_directory",
                    help='''Standard Star Directory''',
                    type=str, default='/work/03730/gregz/maverick/standards')

parser.add_argument("--xc",
                    help='''Science x-centroid''',
                    type=float, default=0.0)

parser.add_argument("--yc",
                    help='''Science y-centroid''',
                    type=float, default=0.0)

parser.add_argument("-mr", "--mask_radius",
                    help='''Mask Radius''',
                    type=float, default=3)

parser.add_argument("-er", "--extraction_radius",
                    help='''Extraction Radius''',
                    type=float, default=3)

parser.add_argument("-sr", "--smoothing_radius",
                    help='''Smothing Radius for background''',
                    type=float, default=2.5)

parser.add_argument("-wi", "--wave_inspect",
                    help='''Inspection Wavelength''',
                    type=float, default=6810.)

parser.add_argument("-iw", "--inspect_width",
                    help='''Inspection Wavelength''',
                    type=float, default=10.)

args = parser.parse_args(args=None)
args.log = setup_logging('extract_standard_from_cube')


sns.set_context('talk')
sns.set_style('ticks')

standard_file = op.join(args.standard_directory, args.standard_name.lower()+'_cal.fits')
A = fits.open(standard_file)

def background_model(cube, maskpix, radius=1.5):
    G = Gaussian2DKernel(radius)
    cube_temp = cube * 1.
    cube_temp[:, maskpix] = np.nan
    for i in np.arange(cube.shape[0]):
        cube_temp[i] = convolve(cube_temp[i], G, boundary='extend')
        #cube_temp[i] = biweight(cube_temp[i])
    cube_temp[np.isnan(cube)] = np.nan
    return cube_temp

def identify_low_background(image, pthresh=7.5):
    thresh = np.nanpercentile(image, pthresh)
    return image < thresh

def get_background_correction(cube, backpix):
    backcor = biweight(cube[:, backpix], axis=1)
    return backcor

def extract_spectrum(cube, pixsel, error=False):
    if error:
        return np.sqrt(np.nansum(cube[:, pixsel]**2, axis=1))
    else:
        return np.nansum(cube[:, pixsel], axis=1)
    
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

def model_telluric_absorption(wave, ratio):
    x, y = (wave, ratio)
    waves = [3670., 4050., 4200., 4420., 4500., 4560., 4750., 5050., 5600.,
             6000., 6650., 6700., 7100., 7450., 7700., 8000., 8450, 8600.,
             8800., 9900., 10100.]
    norms = [1., 1., 1., 1., 1., 1.]
    windows = [[6510., 6650], [6850., 6970.], [7100., 7400.], [7550., 7720.],
               [8000., 8450.], [8850., 9900.]]
    v = []
    for wave in waves:
        v.append(biweight(y[np.abs(x-wave)<10.]))
    I = interp1d(waves, v, kind='linear', fill_value='extrapolate', bounds_error=False)
    continuum = I(x)
    z = np.array(y / continuum)
    windows = [[6510., 6650], [6850., 6970.], [7100., 7400.], [7550., 7720.], [8000., 8450.], [8850., 9900.]]
    mask = np.ones(continuum.shape, dtype=bool)
    for window, norm in zip(windows, norms):
        sel = (x>=window[0]) * (x<=window[1])
        mask[sel] = False
        z[sel] = z[sel] * norm
    z[mask] = 1.
    return z

def get_data_and_grid(filename, xc, yc):
    k = fits.open(filename)
    x = np.arange(k[0].header['NAXIS1'])*k[0].header['CDELT1'] + k[0].header['CRVAL1']
    y = np.arange(k[0].header['NAXIS2'])*k[0].header['CDELT2'] + k[0].header['CRVAL2']
    w = np.arange(k[0].header['NAXIS3'])*k[0].header['CDELT3'] + k[0].header['CRVAL3']
    xgrid, ygrid = np.meshgrid(x, y)
    data = k[0].data
    #data, xgrid, ygrid = make_new_cube(data, xgrid, ygrid, xc, yc)
    return w, xgrid, ygrid, data, k[0].header

def get_aperture_correction(wave, cube, xgrid, ygrid, mask_sel, radius,
                            n=10):
    M1 = Moffat2D(x_0=0., y_0=0.0, alpha=3.5, gamma=1.8)
    for M in [M1]:
        M.alpha.fixed = True
        M.x_0.fixed = True
        M.y_0.fixed = True
        M.gamma.bounds = (0., np.inf)
    fitter = LevMarLSQFitter()
    waves = [np.mean(wc) for wc in np.array_split(wave, n)]
    d = np.sqrt(xgrid**2 + ygrid**2)
    apcor = []
    gamma = []
    for chunk in np.array_split(cube, n, axis=0):
        chunk[chunk==0.] = np.nan
        subimage = biweight(chunk, axis=0)
        mask = mask_sel * np.isfinite(subimage)
        fit = fitter(M1, xgrid[mask], ygrid[mask], subimage[mask])
        gamma.append(fit.gamma.value * 1.)
        model = fit(xgrid, ygrid)
        apcor.append(model[d < radius].sum() / model.sum())
    apcor, waves, gamma = [np.array(x) for x in [apcor, waves, gamma]]
    sel = np.abs(gamma - 1.8)>0.01
    I = interp1d(waves[sel], apcor[sel], kind='linear', fill_value='extrapolate',
                 bounds_error=False)
    return I(wave)

def make_plots(filename, wave, cube, backcube, mask_sel, xgrid, ygrid,
               backpix=None):
    fig, ax = plt.subplots(1, 3, figsize=(12, 5))
    fig.subplots_adjust(right=0.8)
    wsel = np.abs(wave - args.wave_inspect) < args.inspect_width/2.
    G = Gaussian1D(mean=args.wave_inspect, stddev=2.)
    g = G(wave)[wsel] / G(wave)[wsel].sum()
    image = np.nansum(cube[wsel] * g[:, np.newaxis, np.newaxis], axis=0)
    backimage = np.nansum(backcube[wsel] * g[:, np.newaxis, np.newaxis], axis=0)
    subimage = np.nansum((cube[wsel] - backcube[wsel]) *
                         g[:, np.newaxis, np.newaxis], axis=0)
    vmin = np.nanpercentile(image, 2)
    vmax = np.nanpercentile(image, 98)
    norm = ImageNormalize(stretch=LinearStretch(), 
                          vmin=vmin,
                          vmax=vmax)
    cnt = 0
    for a, im, title in zip(ax, [image, backimage, subimage], 
                        ['Data', 'Background Model', 'Background Subtracted']):
        plt.sca(a)
        f_ax10 = plt.gca()
        f_ax10.tick_params(axis='both', which='both', direction='in')
        f_ax10.tick_params(axis='y', which='both', left=True, right=True)
        f_ax10.tick_params(axis='x', which='both', bottom=True, top=True)
        f_ax10.tick_params(axis='both', which='major', length=8, width=2)
        f_ax10.tick_params(axis='both', which='minor', length=5, width=1)
        dx = xgrid[0, 1] - xgrid[0, 0]
        dy = ygrid[1, 0] - ygrid[0, 0]
        aim = a.imshow(im, origin='lower', norm=norm,
                 cmap=plt.get_cmap('Spectral'),
                 extent=[xgrid.min(), xgrid.max()+dx, ygrid.min(), ygrid.max()+dy])
        if cnt == 2:
            cbar_ax = fig.add_axes([0.85, 0.26, 0.05, 0.48])
            cbar = fig.colorbar(aim, cax=cbar_ax)
            cbar.ax.get_yaxis().labelpad = 30
            cbar.ax.set_ylabel(r'F$_{\lambda}$ (ergs / s / cm$^2$ / $\AA$)', rotation=270, fontsize=20)
        cnt += 1
        a.contour(xgrid, ygrid, np.array(mask_sel, dtype=float), levels=[1])

        plt.minorticks_on()
        a.set_title(title)
        a.set_xlim([xgrid.min(), xgrid.max()+dx])
        a.set_ylim([ygrid.min(), ygrid.max()+dy])
    fig.savefig('%s_diagnostic_plot.png' %
                op.basename(filename).split('_cube.fits')[0], dpi=300)

def run_reduction(filename, xc, yc):
    w, xgrid, ygrid, cube, header = get_data_and_grid(filename, xc, yc)
    w, xgrid, ygrid, ecube, header = get_data_and_grid(filename.replace(
                                                       '_cube.fits',
                                                       '_error_cube.fits'),
                                                       xc, yc)
    xc, yc = (0., 0.)
    mask_sel = (np.sqrt((xgrid - xc)**2 + (ygrid - yc)**2) <
                 args.mask_radius)
    extract_sel = (np.sqrt((xgrid - xc)**2 + (ygrid - yc)**2) <
                   args.extraction_radius)
    
    args.log.info('Getting Background Model for %s' % op.basename(filename))
    backcube = background_model(cube, mask_sel, radius=args.smoothing_radius)
    make_plots(filename, w, cube, backcube, mask_sel, xgrid, ygrid)
    try:
        apcor = get_aperture_correction(w, cube-backcube, xgrid, ygrid, mask_sel, 
                                        args.extraction_radius)
    except:
        args.log.warning('Aperture Correction Failed for %s' % op.basename(filename))
        apcor = np.ones(w.shape)
    args.log.info('Getting Spectra for %s' % op.basename(filename))
    spectrum = extract_spectrum(cube-backcube, extract_sel)
    back_spectrum = extract_spectrum(backcube, extract_sel)
    error_spectrum = extract_spectrum(ecube, extract_sel, error=True)
    outname = '%s_new.fits' % op.basename(filename).split('_cube')[0]
    #write_cube(w, xgrid, ygrid, cube-backcube, outname, header)
    outname = '%s_error_new.fits' % op.basename(filename.replace(
                                                       '_cube.fits',
                                                       '_error_cube.fits')).split('_error_cube')[0]
    #write_cube(w, xgrid, ygrid, ecube, outname, header)
    wsel = np.abs(w-3736.) <= 4.
    spectrum[wsel] = np.nan
    error_spectrum[wsel] = np.nan
    back_spectrum[wsel] = np.nan
    return w, spectrum, error_spectrum, back_spectrum, apcor

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

wave, spec, error, back, apcor = run_reduction(args.filename, args.xc,
                                        args.yc)
sel = np.abs(A[1].data['WAVELENGTH']-8200)<800
p = np.polyfit(A[1].data['WAVELENGTH'][sel], 
               np.log10(A[1].data['FLUX'][sel]), 2)
bsel = np.abs(wave-7400)<25
new = np.interp(wave, A[1].data['WAVELENGTH'], A[1].data['FLUX'],
                right=np.nan, left=np.nan)
model = 10**np.polyval(p, wave)
norm = biweight(new[bsel] / model[bsel])
standard_flam = np.zeros(wave.shape)
standard_flam[wave<7400] = new[wave<7400]
standard_flam[wave>=7400] = model[wave>=7400]    
telluric_correction = model_telluric_absorption(wave, spec / apcor /
                                                standard_flam)
telluric_correction[np.isnan(telluric_correction)] = 1.0
T = Table([wave, spec/apcor, error/apcor, back, standard_flam, apcor, telluric_correction],
          names=['wavelength', 'f_lam', 'e_lam', 'b_lam', 's_lam', 'apcor', 'telcor'])
T.write('%s_spec.dat' % (op.basename(args.filename).split('_cube.fits')[0]), 
        format='ascii.fixed_width_two_line', overwrite=True)
