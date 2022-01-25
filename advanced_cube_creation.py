# -*- coding: utf-8 -*-
"""
Created on Fri Oct  5 14:53:10 2018

@author: gregz
"""

import matplotlib
matplotlib.use('agg')
import argparse as ap
import matplotlib.pyplot as plt
import numpy as np
import os.path as op
import seaborn as sns
import sys
import warnings

from astrometry import Astrometry
from astropy import units as U
from astropy.coordinates import SkyCoord
from astropy.convolution import Gaussian1DKernel, convolve, Gaussian2DKernel
from astropy.convolution import interpolate_replace_nans
from astropy.io import fits
from astropy.modeling.models import Gaussian2D, Polynomial1D, Polynomial2D
from astropy.modeling.fitting import LevMarLSQFitter, FittingWithOutlierRemoval
from astropy.stats import biweight_midvariance, sigma_clipped_stats, mad_std
from astropy.stats import sigma_clip
from astropy.table import Table
from input_utils import setup_logging
from matplotlib.ticker import MultipleLocator
from scipy.interpolate import LSQBivariateSpline, griddata
from scipy.signal import medfilt, savgol_filter
from scipy.ndimage import percentile_filter
from sklearn.decomposition import PCA
from math_utils import biweight
from fiber_utils_remedy import find_peaks


warnings.filterwarnings("ignore")



def get_script_path():
    return op.dirname(op.realpath(sys.argv[0]))

DIRNAME = get_script_path()

parser = ap.ArgumentParser(add_help=True)

parser.add_argument("-d", "--directory",
                    help='''base directory for reductions''',
                    type=str, default="")

parser.add_argument("-c", "--caldirectory",
                    help='''cal directory for reductions''',
                    type=str, default="")

parser.add_argument("galaxyname",  help='''Name of Galaxy''', type=str)

parser.add_argument("bluesciobs",
                    help='''e.g., multi_20170126_0000011_exp02_orange.fits''',
                    type=str)

parser.add_argument("redsciobs",
                    help='''e.g., multi_20170126_0000011_exp02_farred.fits''',
                    type=str)

parser.add_argument("skyobs",
                    help='''e.g., multi_20170126_0000011_exp01_farred.fits''',
                    type=str, default=None)

parser.add_argument("ra",  help='''RA of Galaxy''', type=str)

parser.add_argument("dec",  help='''Dec of Galaxy''', type=str)

parser.add_argument("-dw", "--delta_wavelength",
                    help='''Delta Wavelength in linear units for output''',
                    default=None, type=float)

parser.add_argument("-bw", "--blue_wavelength",
                    help='''blue wavelength''',
                    default=None, type=float)

parser.add_argument("-rw", "--red_wavelength",
                    help='''Red Wavelength for object''',
                    default=None, type=float)

parser.add_argument("-no", "--normalization",
                    help='''Absolute Normalization''',
                    default=None, type=float)

parser.add_argument("-ew", "--emission_width",
                    help='''Emission width''',
                    default=5., type=float)

parser.add_argument("-dss", "--dont_subtract_sky",
                    help='''Don't Subtract Sky''',
                    action="count", default=0)

parser.add_argument("-ss", "--simple_sky",
                    help='''Simple Sky''',
                    action="count", default=0)

parser.add_argument("-mac", "--make_arc_cube",
                    help='''Make Arc Cube''',
                    action="count", default=0)

parser.add_argument("-uda", "--use_default_adr",
                    help='''Use Default ADR (only works for side)''',
                    action="count", default=0)

parser.add_argument("-cws", "--correct_wavelength_to_sky",
                    help='''Correct wavelength to sky''',
                    action="count", default=0)

args = parser.parse_args(args=None)
args.log = setup_logging('advance_cube_creation')



def fit_continuum(wv, sky, skip=3, fil_len=95, func=np.array):
    skym_s = 1. * sky
    sky_sm = savgol_filter(skym_s, fil_len, 1)
    allind = np.arange(len(wv), dtype=int)
    for i in np.arange(5):
        mad = np.median(np.abs(sky - sky_sm))
        outlier = func(sky - sky_sm) > 1.5 * mad
        sel = np.where(outlier)[0]
        for j in np.arange(1, skip+1):
            sel = np.union1d(sel, sel + 1)
            sel = np.union1d(sel, sel - 1)
        sel = np.sort(np.unique(sel))
        sel = sel[skip:-skip]
        good = np.setdiff1d(allind, sel)
        skym_s = 1.*sky
        skym_s[sel] = np.interp(wv[sel], wv[good], sky_sm[good])
        sky_sm = savgol_filter(skym_s, fil_len, 1)
    return sky_sm


def make_skyline_model(wave, skylines, norm, dw_pix=None, kernel_size=2.1):
        kernel_size = kernel_size * dw_pix
        skymodel = np.zeros(wave.shape)
        for line in skylines:
            G = (norm * line[1] / np.sqrt((np.pi * 2. * kernel_size**2)) *
                 np.exp(-1. * (wave-line[0])**2 / (2. * kernel_size**2)))
            skymodel += G
        return skymodel


def convert_vac_to_air(skyline):
    s2 = (1e4 / skyline[:, 0])**2
    n = (1 + 0.0000834254 + 0.02406147 / (130 - s2) + 0.00015998 /
         (38.9 - s2))
    skyline[:, 0] = skyline[:, 0] / n
    return skyline


def get_skyline_file(skyline_file):
    V = np.loadtxt(skyline_file)
    V[:, 0] = V[:, 0] * 1e4
    skyline = convert_vac_to_air(V)
    return skyline


def create_header_objection(wave, image, func=fits.ImageHDU):
    hdu = func(np.array(image, dtype='float32'))
    hdu.header['CRVAL1'] = wave[0]
    hdu.header['CRVAL2'] = 1
    hdu.header['CRPIX1'] = 1
    hdu.header['CRPIX2'] = 1
    hdu.header['CTYPE1'] = 'pixel'
    hdu.header['CTYPE2'] = 'pixel'
    hdu.header['CDELT2'] = 1
    hdu.header['CDELT1'] = wave[1] - wave[0]
    return hdu


def get_selection(array1, array2):
    m1 = medfilt(array1, 5)
    m2 = medfilt(array2, 5)
    y1 = np.abs(array1 - m1)
    y2 = np.abs(array2 - m2)
    mad1 = np.nanmedian(y1)
    mad2 = np.nanmedian(y2)
    return (y1 < (5 * mad1)) * (y2 < (5 * mad2))

def solve_system(sci_list, sky_list, x, y, xoff, yoff, sci_image):
    norm1 = np.zeros((sci_list[1].shape[1],))
    norm2 = np.zeros((sci_list[1].shape[1],))
    newsci = sci_list[1] * 0.
    newsky = sky_list[1] * 0.
    C = np.zeros((len(x), 2))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        I = LSQBivariateSpline(x, y, sci_image, np.linspace(-6.0, 6.0, 27),
                               np.linspace(-3.5, 3.5, 15))
    for j in np.arange(sci_list[1].shape[1]):
        if sci_image.ndim == 1:
            xnew = x - xoff[j]
            ynew = y - yoff[j]
            C[:, 0] = I(xnew, ynew, grid=False)
        else:
            C[:, 0] = sci_image[:, j]
        sel = get_selection(sci_list[1][:, j], sky_list[1][:, j])
        sel = (np.isfinite(sci_list[1][:, j]) *
               np.isfinite(sky_list[1][:, j]))
        C[:, 1] = sky_list[1][:, j]
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            sol = np.linalg.lstsq(C[sel], sci_list[1][sel, j])[0]
        norm1[j] = sol[0]
        norm2[j] = sol[1]
        newsci[:, j] = C[:, 0] * sol[0]
        newsky[:, j] = C[:, 1] * sol[1]
    return newsky, newsci, norm1, norm2


def execute_sigma_clip(y, sigma=3):
    try:
        mask = sigma_clip(y, masked=True, maxiters=None, stdfunc=mad_std,
                          sigma=sigma)
    except:
        mask = sigma_clip(y, iters=None, stdfunc=mad_std, sigma=sigma)
    return mask


def correct_amplifier_offsets(y, xp, yp, order=1, kernel=12.):
    i_d = np.sqrt(xp**2 + yp**2)
    i_d_sel = i_d < 8.
    ind = np.where(i_d_sel)[0][np.argmax(y[i_d_sel])]
    xc = xp[ind]
    yc = yp[ind]
    d = np.sqrt((xp-xc)**2 + (yp-yc)**2)
    k = y * 1.
    k[y==0.] = np.nan
    def split_fit(var, ind=140):
        model = k * 0.
        model[:ind] = convolve(k[:ind], Gaussian1DKernel(kernel), boundary='wrap')
        model[ind:] = convolve(k[ind:], Gaussian1DKernel(kernel), boundary='wrap')
        return model
    for i in np.arange(5.):
        model = split_fit(k)
        mstd = mad_std(k - model, ignore_nan=True)
        k[(k-model)>2.5*mstd] = np.nan
    model = split_fit(k)
    loc = 2.5
    k[y==0.] = np.nan
    k[d < loc+1.0] = np.nan
    model = split_fit(k)
    bad = np.isnan(k)
    good = ~bad
    fitter = FittingWithOutlierRemoval(LevMarLSQFitter(), sigma_clip,
                                       stdfunc=mad_std)
    if good.sum() < 5:
        args.log.warning("Cannot Make Correction")
        return np.ones(y.shape), k
    fit, mask = fitter(Polynomial2D(1), xp[good], yp[good], (y/model)[good])
    smodel = fit(xp, yp)
    cor = model * smodel 
    bl, bml = biweight(k[:140]-cor[:140], calc_std=True)
    bl, bmr = biweight(k[140:]-cor[140:], calc_std=True)
    if (bml>0.03) or (bmr>0.03):
        args.log.warning("Cannot Make Correction")
        return np.ones(y.shape), k
    return cor/biweight(cor), k

def estimate_sky(data):
    y = np.mean(data[:, 400:-400], axis=1)
    mask = sigma_clip(y, masked=True, maxiters=None)
    skyfibers = ~mask.mask
    sky = np.median(data[skyfibers], axis=0)
    return sky

def get_pca_sky_residuals(data, ncomponents=5):
    pca = PCA(n_components=ncomponents)
    H = pca.fit_transform(data)
    A = np.dot(H, pca.components_)
    return pca, A

def get_pca_fit_residuals(data, pca):
    coeff = pca.transform(data)
    model = np.dot(coeff, pca.components_)
    return model

def identify_sky_pixels(sky):
    cont = percentile_filter(sky, 5, size=100)
    try:
        mask = sigma_clip(sky - cont, masked=True, maxiters=None,
                          stdfunc=mad_std, sigma=5)
    except:
        mask = sigma_clip(sky - cont, iters=None, stdfunc=mad_std,
                          sigma=5) 
    return mask.mask, cont    

    

def correct_wavelength_to_sky(spectra, skylines, wave, thresh=4.5):
    nfibers = spectra.shape[0]
    sky_wave = skylines['wavelength']
    X = np.arange(len(sky_wave))
    CW = np.zeros((nfibers, len(sky_wave))) * np.nan
    for i in np.arange(nfibers):
        skypixels, back = identify_sky_pixels(spectra[i])
        pl, ps, pv, pw = find_peaks(spectra[i]-back, wave)
        if len(pl):
            V = pw[:, np.newaxis] - sky_wave
            d = np.abs(V)
            D = np.min(d, axis=0)
            L = np.argmin(d, axis=0)
            sel = D < thresh
            y = L[sel]
            x = X[sel]
            CW[i, sel] = V[y, x]
    return CW   

def make_skyline_wave_offset_vs_fiber_plot(wavecorrection_list, utc_list, wave,
                                           galname):
    ml = MultipleLocator(5)
    sns.set_context('talk')
    sns.set_style("ticks", {"xtick.major.size": 8, "ytick.major.size": 8,
                            "xtick.minor.size": 2, "xtick.minor.size": 2})
    plt.figure(figsize=(7, 5))
    plt.gca().set_position([0.2, 0.2, 0.7, 0.7])
    dw = wave[1] - wave[0]
    for wc, utc in zip(wavecorrection_list, utc_list):
        y = np.nanmean(wc, axis=1)
        plt.scatter(np.arange(280), y, alpha=0.1, s=20, edgecolor='none')
        sel = np.isfinite(y)
        P = np.polyval(np.polyfit(np.arange(280)[sel], y[sel], 1),
                       np.arange(280))
        plt.plot(np.arange(280), P, lw=2, label='UTC ' + utc[:-2])
    plt.axes().xaxis.set_minor_locator(ml)
    plt.ylim([dw*-4., dw*4.])
    plt.xlim([-2, 280])
    plt.xlabel(r'Fiber Number')
    plt.ylabel(r'$\lambda$ Offset from Sky Lines ($\AA$)')
    plt.legend()
    plt.savefig('Wave_offset_from_sky_vs_fiber_%s.png' % galname, dpi=300)

def make_skyline_wave_offset_vs_wave_plot(wavecorrection_list, utc_list, wave,
                                           galname, SkyLines):
    ml = MultipleLocator(100)
    sns.set_context('talk')
    sns.set_style("ticks", {"xtick.major.size": 8, "ytick.major.size": 8,
                            "xtick.minor.size": 2, "xtick.minor.size": 2})
    dw = wave[1]-wave[0]
    plt.figure(figsize=(7, 5))
    plt.gca().set_position([0.2, 0.2, 0.7, 0.7])
    for wc, utc in zip(wavecorrection_list, utc_list):
        y = np.nanmean(wc, axis=0)
        plt.scatter(SkyLines['wavelength'], y, alpha=0.4, s=40, edgecolor='none')
        sel = np.isfinite(y)
        P = np.polyval(np.polyfit(SkyLines['wavelength'][sel], y[sel], 2),
                       wave)
        plt.plot(wave, P, lw=2, label='UTC ' + utc[:-2])
    plt.axes().xaxis.set_minor_locator(ml)
    plt.ylim([dw*-4., dw*4.])
    plt.xlim([wave[0], wave[-1]])
    plt.xlabel(r'Wavelength ($\AA$)')
    plt.ylabel(r'$\lambda$ Offset from Sky Lines ($\AA$)')
    plt.legend()
    plt.savefig('Wave_offset_from_sky_vs_wave_%s.png' % galname, dpi=300)

def find_centroid(pos, y):
    grid_x, grid_y = np.meshgrid(np.linspace(-7., 7., (14*5+1)),
                                 np.linspace(-3.5, 3.5, 7*5+1))
    G = Gaussian2DKernel(4.)
    image = griddata(pos[y>0., :2], y[y>0.], (grid_x, grid_y), method='linear')
    i_d = np.sqrt(pos[:, 0]**2 + pos[:, 1]**2)
    i_d_sel = i_d < 5.0
    ind = np.where(i_d_sel)[0][np.nanargmax(y[i_d_sel])]
    xc, yc = (pos[ind, 0], pos[ind, 1])
    d = np.sqrt((grid_x - xc)**2 + (grid_y - yc)**2)
    sel = (d < 1.5) * np.isfinite(image)
    newimage = image * 1.
    newimage[np.isnan(newimage)] = 0.0
    newimage[sel] = np.nan
    back = interpolate_replace_nans(newimage, G)
    back[back==0.] = np.nan
    smooth_back = convolve(back, G, boundary='extend')
    image = image - smooth_back
    xc = np.sum(image[sel] * grid_x[sel]) / np.sum(image[sel])
    yc = np.sum(image[sel] * grid_y[sel]) / np.sum(image[sel])
    a = np.nanmax(image[sel])
    image = image / a
    G = Gaussian2D(x_mean=xc, y_mean=yc, amplitude=1.)
    fitter = FittingWithOutlierRemoval(LevMarLSQFitter(), sigma_clip,
                                       stdfunc=mad_std)
    fit, mask = fitter(G, grid_x[sel], grid_y[sel], image[sel])
    bl, std = biweight(image[~sel], calc_std=True)
    fitquality = False
    if fit.amplitude > 5 * std:
        fitquality = True
    return fit.x_mean.value, fit.y_mean.value, fitquality

def get_adr_curve(pos, data, ordery=1, orderx=0, bins=7):
    x = np.arange(data.shape[1])
    xc = [np.mean(xi) / 1000. for xi in np.array_split(x, bins)]
    yc = [biweight(di, axis=1) for di in np.array_split(data, bins, axis=1)]
    init_y = biweight(data[:, 200:-200], axis=1)
    xP, yP, isgood = find_centroid(pos, init_y)
    xk, yk = ([], [])
    flag = np.zeros((len(xc),), dtype=bool)
    for i, yi in enumerate(yc):
        xp, yp, isgood = find_centroid(pos, yi)
        D = np.sqrt((xP - xp)**2 + (yp - yP)**2)
        if (D < 1.0) and isgood:
            flag[i] = True
        xk.append(xp)
        yk.append(yp)
    fitter = FittingWithOutlierRemoval(LevMarLSQFitter(), sigma_clip,
                                       stdfunc=mad_std)
    fitter = LevMarLSQFitter()
    if flag.sum()>= ordery+1:
        fitx = fitter(Polynomial1D(orderx), np.array(xc)[flag], np.array(xk)[flag])
        fity = fitter(Polynomial1D(ordery), np.array(xc)[flag], np.array(yk)[flag])
        return fitx(x/1000.), fity(x/1000.)
    else:
        args.log.warning('Problem fitting ADC curve.  Not enough points.')
        return np.ones(x.shape) * xP, np.ones(x.shape)* yP
                
def make_cube(xloc, yloc, data, error, Dx, Dy, good, scale, ran, skysub=True,
              radius=0.7):
    '''
    Make data cube for a given ifuslot
    
    Parameters
    ----------
    xloc : 1d numpy array
        fiber positions in the x-direction [ifu frame]
    yloc : 1d numpy array
        fiber positions in the y-direction [ifu frame]
    data : 2d numpy array
        fiber spectra [N fibers x M wavelength]
    Dx : 1d numpy array
        Atmospheric refraction in x-direction as a function of wavelength
    Dy : 1d numpy array
        Atmospheric refraction in y-direction as a function of wavelength
    ftf : 1d numpy array
        Average fiber to fiber normalization (1 value for each fiber)
    scale : float
        Spatial pixel scale in arcseconds
    seeing_fac : float
    
    Returns
    -------
    Dcube : 3d numpy array
        Data cube, corrected for ADR
    xgrid : 2d numpy array
        x-coordinates for data cube
    ygrid : 2d numpy array
        y-coordinates for data cube
    '''
    a, b = data.shape
    N1 = int((ran[1] - ran[0]) / scale)+1
    N2 = int((ran[3] - ran[2]) / scale)+1
    xgrid, ygrid = np.meshgrid(np.linspace(ran[0], ran[1], N1),
                               np.linspace(ran[2], ran[3], N2))
    Dcube = np.zeros((b,)+xgrid.shape)
    Ecube = np.zeros((b,)+xgrid.shape)

    S = np.zeros((data.shape[0], 2))
    D = np.sqrt((xloc - xloc[:, np.newaxis])**2 +
                (yloc - yloc[:, np.newaxis])**2)
    area = scale**2 / (3. / 4. * np.sqrt(3.) * 0.59**2)

    W = np.zeros(D.shape, dtype=bool)
    W[D < radius] = True
    Gp = Gaussian1DKernel(5.)
    c = data * 0.
    for i in np.arange(data.shape[0]):
        c[i] = interpolate_replace_nans(data[i], Gp)
    data = c * 1.
    sky = np.zeros((data.shape[1],))
    for k in np.arange(b):
        S[:, 0] = xloc - Dx[k]
        S[:, 1] = yloc - Dy[k]
        d = np.sqrt(S[:, 0]**2 + S[:, 1]**2)
        sel = (data[:, k] / np.nansum(data[:, k] * W, axis=1)) <= 0.4
        sel *= np.isfinite(data[:, k]) * good * (error[:, k] > 0.)
        if np.sum(sel) > 15:
            if skysub:
                back = biweight(data[sel*(d>4.), k])
            else:
                back = 0.
            sky[k] = back
            Dcube[k, :, :] = (griddata(S[sel], data[sel, k]-back,
                                       (xgrid, ygrid), method='linear') * area)
            Ecube[k, :, :] = (griddata(S[sel], error[sel, k],
                                       (xgrid, ygrid), method='linear') * area)
    return Dcube, Ecube, xgrid, ygrid, sky * area

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

def get_ADR_RAdec(xoff, yoff, astrometry_object):
        '''
        Use an astrometry object to convert delta_x and delta_y into
        delta_ra and delta_dec
        
        Parameters
        ----------
        astrometry_object : Astrometry Class
            Contains the astrometry to convert x and y to RA and Dec
        '''
        tRA, tDec = astrometry_object.tp.wcs_pix2world(xoff, yoff, 1)
        ADRra = ((tRA - astrometry_object.ra0) * 3600. *
                      np.cos(np.deg2rad(astrometry_object.dec0)))
        ADRdec = (tDec - astrometry_object.dec0) * 3600.
        return ADRra, ADRdec
    
def correct_skyline_subtraction(y, sel, pca):
    if sel.sum()>5:
        back = biweight(y[sel])
        res = y - back
        coeff = np.dot(res[sel], pca.components_.T[sel])
        model = np.dot(coeff, pca.components_)
        return model
    else:
        return y * 0.

def make_cor_plot(cor, k, y, name):
    ml = MultipleLocator(5)
    sns.set_context('talk')
    sns.set_style("ticks", {"xtick.major.size": 8, "ytick.major.size": 8,
                            "xtick.minor.size": 2, "xtick.minor.size": 2})
    plt.figure(figsize=(7, 5))
    plt.gca().set_position([0.2, 0.2, 0.7, 0.7])
    norm = biweight(k)
    plt.plot(y/norm)
    plt.plot(k/norm, lw=1)
    plt.plot(cor/norm, 'r-', lw=1)
    
    plt.axes().xaxis.set_minor_locator(ml)
    plt.ylim([0.8, 1.2])
    plt.xlim([0, 280])
    plt.xlabel(r'Fiber Number')
    plt.ylabel(r'Average Value with Respect to Sky')
    plt.legend()
    plt.savefig('cor_%s.png' % name, dpi=300)

def get_arc_pca(spec, pos, good, components=15):
    spec[~good] = 0.
    sky = biweight(spec[good], axis=0)
    mask, cont = identify_sky_pixels(sky)
    ratio = biweight(spec[:, mask] / sky[mask], axis=1)
    cor, k = correct_amplifier_offsets(ratio, pos[:, 0], pos[:, 1])
    nsky = biweight(spec / ratio[:, np.newaxis], axis=0)
    X = (spec / ratio[:, np.newaxis] - nsky)
    X[:, ~mask] = 0.
    X[~good] = 0.
    X = X.swapaxes(0, 1)
    pca, A = get_pca_sky_residuals(X, ncomponents=components)
    return pca

def get_residual_map(data, pca, good):
    res = data * 0.
    for i in np.arange(data.shape[1]):
        sel = good * np.isfinite(data[:, i])
        coeff = np.dot(data[sel, i], pca.components_.T[sel])
        model = np.dot(coeff, pca.components_)
        res[:, i] = model
    return res

def get_cube(SciFits_List, CalFits_List, Pos, scale, ran, skies, waves, cnt,
             cors, def_wave, ems, chns, sky_subtract=True, cal=False,
             scale_sky=False):
    F = []
    info = []
    if cors is None:
        cors = [None] * len(skies)
    for _scifits, _calfits, P, skY, cor, wave, em, chn in zip(SciFits_List,
                                                     CalFits_List, Pos, skies,
                                                     cors, waves, ems, chns):
        args.log.info('Working on reduction for %s' % _scifits.filename())
        if not cal:
            SciSpectra = _scifits[0].data
            SciError = _scifits[3].data
        else:
            SciSpectra = _calfits['arcspec'].data
            SciError = 0. * SciSpectra
            sel = SciSpectra > 0.
            SciError[sel]= np.sqrt(SciSpectra[sel]/np.sqrt(2) + 3**2*2.)
            SciError[~sel] = np.sqrt(3**2*2.)
        if cor is not None:
            SciSpectra /= cor[:, np.newaxis]
            SciError /= cor[:, np.newaxis]
        good = (SciSpectra == 0.).sum(axis=1) < 200
        pos = _scifits[5].data
        pca = get_arc_pca(_calfits['arcspec'].data, pos, good,
                          components=75)
        if chn == 'orange':
            print(chn)
            SciSpectra[:140] = SciSpectra[:140] / 1.025
            SciSpectra[140:] = SciSpectra[140:] / 0.975
        if cor is None:
            sel = (SciSpectra == 0.).sum(axis=1) < 200
            y = biweight(SciSpectra[:, 200:-200], axis=1)
            cor, k = correct_amplifier_offsets(y, pos[:, 0], pos[:, 1])
            mask = execute_sigma_clip(y / cor)
            selm = mask.mask * sel
            d = np.sqrt((pos[:, 0, np.newaxis,] - pos[:, 0])**2 +
                        (pos[:, 1, np.newaxis,] - pos[:, 1])**2)
            for j in np.where(selm)[0]:
                selm = selm + (d[j] < 3.)
            sel = sel * ~selm
            Sky = biweight(SciSpectra[sel] / cor[sel, np.newaxis],
                           axis=0)
            y = biweight(SciSpectra[:, 200:-200] / Sky[200:-200], axis=1)
            cor, k = correct_amplifier_offsets(y, pos[:, 0], pos[:, 1])
            make_cor_plot(cor, k, y, op.basename(_scifits.filename()))
            SciSpectra /= cor[:, np.newaxis]
            SciError /= cor[:, np.newaxis]
        if sky_subtract:
            if skY is None:
                quick_sky = biweight(SciSpectra, axis=0)
                mask, cont = identify_sky_pixels(quick_sky)
                std_sky = mad_std((quick_sky-cont)[~mask])
                loc, values = find_peaks((quick_sky-cont), thresh=15*std_sky)
                loc = np.array(np.round(loc), dtype=int)
                loc = loc[(loc>10) * (loc<len(quick_sky)-10)]
                # Remove Continuum (gaussian filter)
                skysub_rect = SciSpectra - quick_sky[np.newaxis, :]
                Dummy = skysub_rect
                for i in np.arange(-6, 7):
                    Dummy[:, loc+i] = np.nan
                Smooth = Dummy * np.nan
                for i in np.arange(Dummy.shape[0]):
                    Smooth[i] = convolve(Dummy[i], Gaussian1DKernel(2.0), boundary='extend')
                    while np.isnan(Smooth[i]).sum():
                        Smooth[i] = interpolate_replace_nans(Smooth[i], Gaussian1DKernel(4.0))
                res = get_residual_map(skysub_rect-Smooth, pca, good)
                skyval = quick_sky+res
                if args.simple_sky:
                    d = np.sqrt(pos[:, 0]**2 + pos[:, 1]**2)
                    skyval = np.ones((280, 1)) * biweight(SciSpectra[d>3.], axis=0)[np.newaxis, :]
                SciSpectra = SciSpectra - skyval
            else:
                args.log.info('Using other sky for subtraction')
                sel = (SciSpectra == 0.).sum(axis=1) < 200
                if em is not None:
                    wsel = np.abs(wave-em)<args.emission_width
                    if wsel.sum() > 0:
                        y = biweight(SciSpectra[:, wsel], axis=1)
                y = biweight(SciSpectra[:, 200:-200], axis=1)
                cor, k = correct_amplifier_offsets(y, pos[:, 0], pos[:, 1])
                mask = execute_sigma_clip(y / cor)
                selm = mask.mask * sel
                d = np.sqrt((pos[:, 0, np.newaxis,] - pos[:, 0])**2 +
                            (pos[:, 1, np.newaxis,] - pos[:, 1])**2)
                for j in np.where(selm)[0]:
                    selm = selm + (d[j] < 3.)
                sel = sel * ~selm
                args.log.info('Number of fibers for sky: %i' % sel.sum())
                ratio = biweight(SciSpectra[sel] / skY[sel], axis=0)
                skyval = skY * ratio
                SciSpectra = SciSpectra - skyval
        else:
            skyval = 0. * SciSpectra
        SciSpectra[~good] = 0.
        zcube, ecube, xgrid, ygrid, scisky = make_cube(P[0], P[1],
                                               SciSpectra, SciError,
                                               P[2], P[3], good,
                                               scale, ran, skysub=False)
        scube = zcube * 0.
        scube, secube, xgrid, ygrid, dummy = make_cube(P[0], P[1],
                                               skyval, 1.*np.isfinite(skyval),
                                               P[2], P[3], good,
                                               scale, ran, skysub=False)
        scube = scube + scisky[:, np.newaxis, np.newaxis]
        skysub_cube = zcube
        newcube = np.zeros((len(def_wave), zcube.shape[1], zcube.shape[2]))
        newerrcube = np.zeros((len(def_wave), zcube.shape[1], zcube.shape[2]))
        skycube = np.zeros((len(def_wave), zcube.shape[1], zcube.shape[2]))

        for j in np.arange(zcube.shape[1]):
            for k in np.arange(zcube.shape[2]):
                newcube[:, j, k] = np.interp(def_wave, wave, skysub_cube[:, j, k],
                                             left=np.nan, right=np.nan)
                newerrcube[:, j, k] = np.interp(def_wave, wave, ecube[:, j, k],
                                             left=np.nan, right=np.nan)
                skycube[:, j, k] = np.interp(def_wave, wave, scube[:, j, k],
                                             left=np.nan, right=np.nan)
        info.append([newcube, newerrcube, skycube, xgrid, ygrid])
        
        F.append([])
        cnt += 1
    return F, info

def get_norm(cube, xgrid, ygrid, wave, dist=3.):
    xl = 250
    xh = 400
    image = biweight(cube[xl:xh], axis=0)
    image = image / np.nansum(image)
    fits.PrimaryHDU(image).writeto('debug.fits', overwrite=True)
    G = Gaussian2D(x_mean=0.0, y_mean=0.0, amplitude=0.005)
    G.x_mean.bounds = (-0.25, 0.25)
    G.y_mean.bounds = (-0.25, 0.25)
    G.x_stddev.bounds = (0.5, 5.)
    G.y_stddev.bounds = (0.5, 5.)
    G.amplitude.bounds = (0., 1.)
    fitter = FittingWithOutlierRemoval(LevMarLSQFitter(), sigma_clip,
                                       stdfunc=mad_std)
    X, Y = (xgrid.ravel(), ygrid.ravel())
    distsel = np.sqrt(X**2 + Y**2) < dist
    totsel = distsel * np.isfinite(image.ravel()) * (image.ravel() > 1e-6)
    rejected, fit = fitter(G, X[totsel], Y[totsel], image.ravel()[totsel])
    area = fit.amplitude / np.sqrt(2. * np.pi * fit.x_stddev * fit.y_stddev)
    print(fit)
    return area


def main():
    bluesciobs = [x.replace(' ', '') for x in args.bluesciobs.split(',')]
    redsciobs = [x.replace(' ', '') for x in args.redsciobs.split(',')]

    skyobs = [x.replace(' ', '') for x in args.skyobs.split(',')]
    side_dict = {'uv': 'LRS2B', 'orange': 'LRS2B', 'red': 'LRS2R', 'farred': 'LRS2R'}
    otherchannel_dict = {'uv': 'orange', 'orange': 'uv', 'red': 'farred',
                         'farred': 'red'}
    em_dict = {'uv': args.blue_wavelength, 'orange': args.blue_wavelength,
               'red': args.red_wavelength, 'farred': args.red_wavelength}
    sciobs = []
    B, R = (False, False)
    for OBS in [bluesciobs, redsciobs]:
        for obs in OBS:
            if obs != '':    
                channel = obs.split('_')[-1][:-5]
                otherchannel = otherchannel_dict[channel]
                side = side_dict[channel]
                if side == 'LRS2B':
                    B = True
                if side == 'LRS2R':
                    R = True
                sciobs.append(obs)
                sciobs.append(obs.replace(channel, otherchannel))
    skyobs2 = []
    for OBS in [skyobs]:
        for obs in OBS:
            if obs != '':    
                channel = obs.split('_')[-1][:-5]
                otherchannel = otherchannel_dict[channel]
                skyobs2.append(obs)
                skyobs2.append(obs.replace(channel, otherchannel))
    skyobs = skyobs2
    scale = 0.25
    ran = [-3.6, 3.6, -6.4, 6.4]
    
# =============================================================================
# Fitting for Sky
# =============================================================================
    sky, cor, chan = ([], [], [])
    SkyFits_List = []

    for _skyobs in skyobs:
        channel = _skyobs.split('_')[-1][:-5]
        chan.append(channel)
        SkyFits_List.append(fits.open(op.join(args.directory, _skyobs)))
        if channel == 'orange':
            SkyFits_List[-1][0].data[:140] = SkyFits_List[-1][0].data[:140] / 1.025
            SkyFits_List[-1][0].data[:140] = SkyFits_List[-1][0].data[:140] / 0.975
        args.log.info('Sky observation: %s loaded' % (_skyobs))
        SkySpectra = SkyFits_List[-1][0].data
        P = SkyFits_List[-1][5].data
        sel = (SkySpectra == 0.).sum(axis=1) < 200
        y = biweight(SkySpectra[:, 200:-200], axis=1)
        correction, k = correct_amplifier_offsets(y, P[:, 0], P[:, 1])
        mask = execute_sigma_clip(y / correction)
        selm = mask.mask * sel
        sel = sel * ~selm
        y = biweight(SkySpectra[:, 200:-200] /
                     biweight(SkySpectra[sel, 200:-200], axis=0), axis=1)
        correction, k = correct_amplifier_offsets(y, P[:, 0], P[:, 1])
        cor.append(correction)
        make_cor_plot(correction, k, y, op.basename(SkyFits_List[-1].filename()))
        X = SkySpectra / correction[:, np.newaxis]
        X[SkyFits_List[-1][3].data==0.] = np.nan
        sky.append(X)
        cor.append(correction)
    sky_dict = {'uv': None, 'orange': None, 'red': None, 'farred': None}
    cor_dict = {'uv': None, 'orange': None, 'red': None, 'farred': None}

    for uchan in np.unique(chan):
        sky_dict[uchan] = np.array(np.nanmean([sk for sk, ch in zip(sky, chan)
                                            if ch == uchan], axis=0))
        cor_dict[uchan] = np.array(np.mean([sk for sk, ch in zip(cor, chan)
                                            if ch == uchan], axis=0))
# =============================================================================
# Reading cooridinates for Astrometry
# =============================================================================
    try:
        S = SkyCoord(args.ra, args.dec)
    except:
        args.log.error('Coordinates need to be in format XXhXXmXX.Xs and ' 
                       '+/-XXdXXmXX.Xs')
        sys.exit(1)
        
# =============================================================================
# Gathering science and calibration information        
# =============================================================================
    ran_list = []
    SciFits_List = []
    CalFits_List = []
    T = []
    Pos = []
    info = []
    skies = []
    cors = []
    waves = []
    ems = []
    chns = []
    for j, _sciobs in enumerate(sciobs):
        channel = _sciobs.split('_')[-1][:-5]
        side = side_dict[channel]
        sky = sky_dict[channel]
        cor = cor_dict[channel]
        SciFits_List.append(fits.open(op.join(args.directory, _sciobs)))
        args.log.info('Science observation: %s loaded' % (_sciobs))
        date = _sciobs.split('_')[1]
        calname = 'cal_%s_%s.fits' % (date, channel)
        CalFits_List.append(fits.open(op.join(args.caldirectory, calname)))
        args.log.info('Cal observation: %s loaded' % calname)
        channel_dict = {'uv': 'BL', 'orange': 'BR', 'red': 'RL', 'farred': 'RR'}
        specinit = channel_dict[channel]
        darfile = op.join(DIRNAME, 'lrs2_config/dar_%s.dat' % specinit)
        T.append(Table.read(darfile, format='ascii.fixed_width_two_line'))
        iwave = SciFits_List[0][6].data[0]
        wave = SciFits_List[-1][6].data[0]
        waves.append(wave)
        pos = SciFits_List[0][5].data * 1.
        pos[:, 0] = SciFits_List[0][5].data[:, 0] * 1.
        pos[:, 1] = SciFits_List[0][5].data[:, 1] * 1.
        wave_0 = np.mean(iwave)
        xint =  np.interp(wave_0, T[0]['wave'], T[0]['x_0'])
        yint =  np.interp(wave_0, T[0]['wave'], T[0]['y_0'])
        xoff = np.interp(wave, T[-1]['wave'], T[-1]['x_0']) - xint
        yoff = np.interp(wave, T[-1]['wave'], T[-1]['y_0']) - yint
        if side == 'LRS2B':
            order = 3
            bins=11
        else:
            order = 1
            bins = 5
        if not args.use_default_adr:
            xoff, yoff = get_adr_curve(pos, SciFits_List[-1][0].data, 
                                       ordery=order, bins=bins)
        args.log.info('%s: %0.2f, %0.2f' % (_sciobs, np.mean(xoff), np.mean(yoff)))
        xc, yc = (0., 0.)
        A = Astrometry(S.ra.deg, S.dec.deg, SciFits_List[-1][0].header['PARANGLE'],
                       xc, yc, x_scale=1., y_scale=1., kind='lrs2')
        raoff, decoff = get_ADR_RAdec(xoff+xc, yoff+yc, A)
        ra, dec = A.tp.wcs_pix2world(pos[:, 0], pos[:, 1], 1)
        delta_ra = np.cos(np.deg2rad(S.dec.deg)) * (ra - S.ra.deg) * 3600.
        delta_dec = (dec - S.dec.deg) * 3600.
        ran1 = [np.min(delta_ra), np.max(delta_ra), np.min(delta_dec),
                np.max(delta_dec)]
        ran_list.append(ran1)
        Pos.append([delta_ra, delta_dec, raoff, decoff])
        skies.append(sky)
        cors.append(cor)
        ems.append(em_dict[channel])
        chns.append(channel)
    ran_array = np.array(ran_list)
    rmax = np.max(ran_array, axis=0)
    rmin = np.min(ran_array, axis=0)
    ran = [np.floor(rmin[0]/scale)*scale, np.ceil(rmax[1]/scale)*scale,
           np.floor(rmin[2]/scale)*scale, np.ceil(rmax[3]/scale)*scale]   
    args.log.info('Cube limits - x: [%0.2f, %0.2f], y: [%0.2f, %0.2f]' %
                  (ran[0], ran[1], ran[2], ran[3]))

#    if args.correct_wavelength_to_sky:
#        skylinefile = op.join(DIRNAME, 'lrs2_config/skylines_%s.dat' % channel)
#        SkyLines = Table.read(skylinefile, format='ascii.fixed_width_two_line')
#        wavecorrection_list, utc_list = ([], [])
#        for _scifits in SciFits_List:
#            args.log.info('Correcting wavelength for %s' % _scifits.filename())
#            SciSpectra = _scifits[0].data
#            utc_list.append(_scifits[0].header['UT'])
#            CW = correct_wavelength_to_sky(SciSpectra, SkyLines, wave)
#            wavecorrection_list.append(CW)
#        for _skyfits in SkyFits_List:
#            args.log.info('Correcting wavelength for %s' % _skyfits.filename())
#            SkySpectra = _skyfits[0].data
#            utc_list.append(_skyfits[0].header['UT'])
#            CW = correct_wavelength_to_sky(SkySpectra, SkyLines, wave)
#            wavecorrection_list.append(CW)
#        make_skyline_wave_offset_vs_fiber_plot(wavecorrection_list, utc_list,
#                                               wave, args.galaxyname)
#        make_skyline_wave_offset_vs_wave_plot(wavecorrection_list, utc_list,
#                                               wave, args.galaxyname,
#                                               SkyLines)
    
    
    

    if args.delta_wavelength is None:
        dw = np.min(np.diff(np.array(waves), axis=1))
    else:
        dw = args.delta_wavelength
    lw = np.min(np.array(waves))
    hw = np.max(np.array(waves))
    def_wave = np.arange(lw, hw+dw, dw)
    args.log.info('Wavelength Range: %0.1f - %0.1f A, %0.2f A' % (lw, hw, dw))
    ############################### Science ###################################
    cnt = 0
    chan_list = [_sciobs.split('_')[-1][:-5] for _sciobs in sciobs]
    side_list = [side_dict[_sciobs.split('_')[-1][:-5]] for _sciobs in sciobs]
    F, info = get_cube(SciFits_List, CalFits_List, Pos, scale, ran, skies, 
                       waves, cnt, cors, def_wave, ems, chns, sky_subtract=True)
    xgrid = info[0][3]
    ygrid = info[0][4]
    if B and R:
        l = []
        for i, chan, side in zip(info, chan_list, side_list):
            if chan == 'farred':
                norm = 0.93
            else:
                norm = 1.0
            if side == 'LRS2R':
                l.append(i[0] * norm)
        bcube = np.nanmean(np.array([i[0] for i, side in zip(info, side_list)
                                     if side=='LRS2B']), axis=0)
        rcube = np.nanmean(np.array(l), axis=0)
        becube = np.nanmean(np.array([i[1] for i, side in zip(info, side_list)
                                     if side=='LRS2B']), axis=0)
        recube = np.nanmean(np.array([i[1] for i, side in zip(info, side_list)
                                     if side=='LRS2R']), axis=0)
        bscube = np.nanmean(np.array([i[2] for i, side in zip(info, side_list)
                                     if side=='LRS2B']), axis=0)
        rscube = np.nanmean(np.array([i[2] for i, side in zip(info, side_list)
                                     if side=='LRS2R']), axis=0)
        gwave = np.array(np.isfinite(bcube) * np.isfinite(rcube), dtype=float)
        boverlap = np.nansum(bcube * gwave, axis=(1, 2))
        roverlap = np.nansum(rcube * gwave, axis=(1, 2))
        norm = biweight(boverlap / roverlap)
        args.log.info('Normalziation of Blue Side to Red Side: %0.2f' % norm)
        zcube = np.nanmean([bcube, rcube*norm], axis=0)
        ecube = np.sqrt(np.nanmean([becube**2, recube**2*norm**2], axis=0))
        scube = np.nanmean([bscube, rscube*norm], axis=0)
    else:
        l = []
        for i, chan in zip(info, chan_list):
            if chan == 'farred':
                norm = 0.93
            else:
                norm = 1.0
            l.append(i[0] * norm)
        zcube = np.nanmean(np.array(l), axis=0)
        ecube = np.sqrt(np.nanmean(np.array([i[1] for i in info])**2, axis=0))
        scube = np.nanmean(np.array([i[2] for i in info]), axis=0)
    Header = SciFits_List[0][0].header
    if B and R:
        name = 'LRS2B+R'
    if B and (not R):
        name = 'LRS2B'
    if R and (not B):
        name = 'LRS2R'
    outname = '%s_%s_cube.fits' % (args.galaxyname,  name)
    eoutname = '%s_%s_error_cube.fits' % (args.galaxyname,  name)
    soutname = '%s_%s_sky_cube.fits' % (args.galaxyname,  name)
    if args.make_arc_cube:
        F, info = get_cube(SciFits_List, CalFits_List, Pos, scale, ran, skies, 
                       waves, cnt, cors, def_wave, ems, sky_subtract=False,
                       cal=True)
        acube = np.nanmean(np.array([i[0] for i in info]), axis=0)
        aoutname = '%s_%s_lamp_cube.fits' % (args.galaxyname,  name)
        write_cube(def_wave, xgrid, ygrid, acube, aoutname, Header)
    if args.normalization is not None:
        zcube[:] *= args.normalization
        ecube[:] *= args.normalization
        scube[:] *= args.normalization
    write_cube(def_wave, xgrid, ygrid, zcube, outname, Header)
    write_cube(def_wave, xgrid, ygrid, ecube, eoutname, Header)
    write_cube(def_wave, xgrid, ygrid, scube, soutname, Header)
    



main()
