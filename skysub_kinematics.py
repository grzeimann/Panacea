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
from astropy.convolution import Gaussian1DKernel, convolve
from astropy.convolution import interpolate_replace_nans
from astropy.io import fits
from astropy.modeling.models import Gaussian2D
from astropy.modeling.fitting import LevMarLSQFitter, FittingWithOutlierRemoval
from astropy.stats import biweight_midvariance, sigma_clipped_stats, mad_std
from astropy.stats import sigma_clip
from astropy.table import Table
from input_utils import setup_logging
from matplotlib.ticker import MultipleLocator
from scipy.interpolate import LSQBivariateSpline, griddata
from scipy.signal import medfilt, savgol_filter
from sklearn.decomposition import PCA
from math_utils import biweight


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

parser.add_argument("sciobs",
                    help='''e.g., multi_20170126_0000011_exp02_farred.fits''',
                    type=str)

parser.add_argument("skyobs",
                    help='''e.g., multi_20170126_0000011_exp01_farred.fits''',
                    type=str, default=None)

parser.add_argument("ra",  help='''RA of Galaxy''', type=str)

parser.add_argument("dec",  help='''Dec of Galaxy''', type=str)

parser.add_argument("-dss", "--dont_subtract_sky",
                    help='''Don't Subtract Sky''',
                    action="count", default=0)

parser.add_argument("-ds", "--dont_correct_wavelength_to_sky",
                    help='''Don't correct wavelength to sky if used''',
                    action="count", default=0)

args = parser.parse_args(args=None)
args.log = setup_logging('skysub_jonelle')


def make_frame(xloc, yloc, data, wave, dw, Dx, Dy, wstart=5700.,
               wend=5800., scale=0.4, seeing_fac=1.3):
    seeing = seeing_fac * scale
    a, b = data.shape
    x = np.arange(xloc.min()-scale,
                  xloc.max()+1*scale, scale)
    y = np.arange(yloc.min()-scale,
                  yloc.max()+1*scale, scale)
    xgrid, ygrid = np.meshgrid(x, y)
    zgrid = np.zeros((b,)+xgrid.shape)
    area = 3. / 4. * np.sqrt(3.) * 0.59**2
    for k in np.arange(b):
        sel = np.isfinite(data[:, k])
        D = np.sqrt((xloc[:, np.newaxis, np.newaxis] - Dx[k] - xgrid)**2 +
                    (yloc[:, np.newaxis, np.newaxis] - Dy[k] - ygrid)**2)
        W = np.exp(-0.5 / (seeing/2.35)**2 * D**2)
        N = W.sum(axis=0)
        zgrid[k, :, :] = ((data[sel, k][:, np.newaxis, np.newaxis] *
                           W[sel]).sum(axis=0) / N / scale**2 / area)
    wi = np.searchsorted(wave, wstart, side='left')
    we = np.searchsorted(wave, wend, side='right')

    zimage = biweight_location(zgrid[wi:we+1], axis=(0,))
    return zgrid, zimage, xgrid, ygrid


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

def find_peaks(y, wave, thresh=8.):
    def get_peaks(flat, XN):
        YM = np.arange(flat.shape[0])
        inds = np.zeros((3, len(XN)))
        inds[0] = XN - 1.
        inds[1] = XN + 0.
        inds[2] = XN + 1.
        inds = np.array(inds, dtype=int)
        Peaks = (YM[inds[1]] - (flat[inds[2]] - flat[inds[0]]) /
                 (2. * (flat[inds[2]] - 2. * flat[inds[1]] + flat[inds[0]])))
        return Peaks
    diff_array = y[1:] - y[:-1]
    loc = np.where((diff_array[:-1] > 0.) * (diff_array[1:] < 0.))[0]
    peaks = y[loc+1]
    std = np.sqrt(biweight_midvariance(y))
    loc = loc[peaks > (thresh * std)]+1
    peak_loc = get_peaks(y, loc)
    peaks = y[np.round(peak_loc).astype(int)]
    peak_wave = np.interp(peak_loc, np.arange(len(wave)), wave)
    return peak_loc, peaks/std, peaks, peak_wave

def execute_sigma_clip(y):
    try:
        mask = sigma_clip(y, masked=True, maxiters=None, stdfunc=mad_std)
    except:
        mask = sigma_clip(y, iters=None, stdfunc=mad_std)
    return mask


def correct_amplifier_offsets(y, order=1):
    x = np.arange(len(y))
    maskL = execute_sigma_clip(y[:140])
    maskR = execute_sigma_clip(y[140:])
    modelL = np.polyval(np.polyfit(x[:140][~maskL.mask],
                                   y[:140][~maskL.mask], order), x[:140])
    modelR = np.polyval(np.polyfit(x[140:][~maskR.mask],
                                   y[140:][~maskR.mask], order), x[140:])
    avg = np.mean(np.hstack([modelL, modelR]))
    return np.hstack([modelL / avg, modelR / avg])

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
    G = Gaussian1DKernel(20.0)
    cont = convolve(sky, G)
    try:
        mask = sigma_clip(sky - cont, masked=True, maxiters=None,
                          stdfunc=mad_std)
    except:
        mask = sigma_clip(sky - cont, iters=None, stdfunc=mad_std) 
    for i in np.arange(5):
        nsky = sky * 1.
        mask.mask[1:] += mask.mask[:-1]
        mask.mask[:-1] += mask.mask[1:]
        nsky[mask.mask] = np.nan
        cont = convolve(nsky, G, boundary='extend')
        while np.isnan(cont).sum():
            cont = interpolate_replace_nans(cont, G)
        try:
            mask = sigma_clip(sky - cont, masked=True, maxiters=None,
                              stdfunc=mad_std)
        except:
            mask = sigma_clip(sky - cont, iters=None, stdfunc=mad_std) 
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
    image = griddata(pos[y>0., :2], y[y>0.], (grid_x, grid_y), method='cubic')
    xc, yc = (pos[np.nanargmax(y), 0], pos[np.nanargmax(y), 1])
    d = np.sqrt((grid_x - xc)**2 + (grid_y - yc)**2)
    sel = (d < 2.) * np.isfinite(image)
    xc = np.sum(image[sel] * grid_x[sel]) / np.sum(image[sel])
    yc = np.sum(image[sel] * grid_y[sel]) / np.sum(image[sel])
    return xc, yc

def make_cube(xloc, yloc, data, error, Dx, Dy, good, scale, ran,
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
    for k in np.arange(b):
        S[:, 0] = xloc - Dx[k]
        S[:, 1] = yloc - Dy[k]
        sel = (data[:, k] / np.nansum(data[:, k] * W, axis=1)) <= 0.4
        sel *= np.isfinite(data[:, k]) * good * (error[:, k] > 0.)
        if np.sum(sel) > 15:
            Dcube[k, :, :] = (griddata(S[sel], data[sel, k],
                                       (xgrid, ygrid), method='cubic') * area)
            Ecube[k, :, :] = (griddata(S[sel], error[sel, k],
                                       (xgrid, ygrid), method='cubic') * area)
    return Dcube, Ecube, xgrid, ygrid

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

def get_cube(SciFits_List, Pos, scale, ran, sky, wave, cnt, cor=None,
             sky_subtract=True, cal=False):
    F = []
    info = []
    for _scifits, P in zip(SciFits_List, Pos):
        args.log.info('Working on reduction for %s' % _scifits.filename())
        if not cal:
            SciSpectra = _scifits[0].data
            SciError = _scifits[3].data
        else:
            SciSpectra = _scifits['arcspec'].data
            SciError = 0. * SciSpectra
            sel = SciSpectra > 0.
            SciError[sel]= np.sqrt(SciSpectra[sel]/np.sqrt(2) + 3**2*2.)
            SciError[~sel] = np.sqrt(3**2*2.)
        if cor is not None:
            SciSpectra /= cor[:, np.newaxis]
            SciError /= cor[:, np.newaxis]
        good = (SciSpectra == 0.).sum(axis=1) < 200
        zcube, ecube, xgrid, ygrid = make_cube(P[0], P[1],
                                               SciSpectra, SciError,
                                               P[2], P[3], good,
                                               scale, ran)
        d = np.sqrt(P[0]**2 + P[1]**2)
        skysel = (d > np.max(d) - 1.5)
        pixsel = np.zeros(xgrid.shape, dtype=bool)
        for fib in np.where(skysel)[0]:    
            D = np.sqrt((xgrid-P[0][fib])**2 + (ygrid-P[1][fib])**2)
            pixsel += D < 0.6
        if sky_subtract:
            scisky = np.nanmedian(zcube[:, pixsel], axis=1)
            if sky is not None:
                ratio = biweight(zcube[:, pixsel] / sky[:, np.newaxis], axis=1)
                scisky = sky * ratio
        else:
            scisky = np.zeros((SciSpectra.shape[1],))
        skysub_cube = zcube - scisky[:, np.newaxis, np.newaxis]
        info.append([skysub_cube, ecube, xgrid, ygrid])
        
        # Subtract sky in fiber space rather than on the cube
        if sky_subtract:
            skytemp = np.nanmedian(SciSpectra[skysel], axis=0)
            if sky is not None:
                ratio = biweight(SciSpectra[skysel] / sky, axis=0)
                skytemp = sky * ratio
        else:
            skytemp = np.zeros((SciSpectra.shape[1],))
        if cnt == 0:
            func = fits.PrimaryHDU
        else:
            func = fits.ImageHDU
        f1 = create_header_objection(wave, SciSpectra - skytemp,
                                     func=func)
        F.append(f1)
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
    sciobs = [x.replace(' ', '') for x in args.sciobs.split(',')]
    skyobs = [x.replace(' ', '') for x in args.skyobs.split(',')]
    SciFits_List = []
    for _sciobs in sciobs:
        SciFits_List.append(fits.open(op.join(args.directory, _sciobs)))
        args.log.info('Science observation: %s loaded' % (_sciobs))
    SkyFits_List = []
    for _skyobs in skyobs:
        skyflag = _skyobs.lower() != 'none'
        if skyflag:
            SkyFits_List.append(fits.open(op.join(args.directory, _skyobs)))
            args.log.info('Sky observation: %s loaded' % (_skyobs))
    CalFits_List = []
    for _sciobs in sciobs:
        date = _sciobs.split('_')[1]
        channel = _sciobs.split('_')[-1][:-5]
        calname = 'cal_%s_%s.fits' % (date, channel)
        CalFits_List.append(fits.open(op.join(args.caldirectory, calname)))
        args.log.info('Cal observation: %s loaded' % calname)
    
    # Basic Info Dump
    wave = SciFits_List[0][6].data[0]
    channel = sciobs[0].split('_')[-1][:-5]
    channel_dict = {'uv': 'BL', 'orange': 'BR', 'red': 'RL', 'farred': 'RR'}
    specinit = channel_dict[channel]
    darfile = op.join(DIRNAME, 'lrs2_config/dar_%s.dat' % specinit)
    T = Table.read(darfile, format='ascii.fixed_width_two_line')
    
    if not args.dont_correct_wavelength_to_sky:
        skylinefile = op.join(DIRNAME, 'lrs2_config/skylines_%s.dat' % channel)
        SkyLines = Table.read(skylinefile, format='ascii.fixed_width_two_line')
        wavecorrection_list, utc_list = ([], [])
        for _scifits in SciFits_List:
            args.log.info('Correcting wavelength for %s' % _scifits.filename())
            SciSpectra = _scifits[0].data
            utc_list.append(_scifits[0].header['UT'])
            CW = correct_wavelength_to_sky(SciSpectra, SkyLines, wave)
            wavecorrection_list.append(CW)
        for _skyfits in SkyFits_List:
            args.log.info('Correcting wavelength for %s' % _skyfits.filename())
            SkySpectra = _skyfits[0].data
            utc_list.append(_skyfits[0].header['UT'])
            CW = correct_wavelength_to_sky(SkySpectra, SkyLines, wave)
            wavecorrection_list.append(CW)
        make_skyline_wave_offset_vs_fiber_plot(wavecorrection_list, utc_list,
                                               wave, args.galaxyname)
        make_skyline_wave_offset_vs_wave_plot(wavecorrection_list, utc_list,
                                               wave, args.galaxyname,
                                               SkyLines)
    
    # Get Fiber Positions and ADR Correction
    pos = SciFits_List[0][5].data * 1.
    pos[:, 0] = SciFits_List[0][5].data[:, 0] * 1.
    pos[:, 1] = SciFits_List[0][5].data[:, 1] * 1.
    wave_0 = np.mean(wave)
    xoff = (np.interp(wave, T['wave'], T['x_0']) -
            np.interp(wave_0, T['wave'], T['x_0']))
    yoff = (np.interp(wave, T['wave'], T['y_0']) -
            np.interp(wave_0, T['wave'], T['y_0']))
    
    sky = []
    cor = []
    for _skyfits in SkyFits_List:
        SkySpectra = _skyfits[0].data
        sel = (SkySpectra == 0.).sum(axis=1) < 200
        y = np.nanmedian(SkySpectra[:, 410:440], axis=1)
        mask = execute_sigma_clip(y)
        sel = sel * ~mask.mask
        correction = correct_amplifier_offsets(y)
        cor.append(correction)
        sky.append(np.median(SkySpectra[sel], axis=0))
    if len(sky) > 0:
        sky = np.sum(sky, axis=0)
        cor = np.mean(cor, axis=0)
    else:
        sky = None
        cor = None
    info = []
    scale = 0.25
    ran = [-3.6, 3.6, -6.4, 6.4]
    ran_list = []
    try:
        S = SkyCoord(args.ra, args.dec)
    except:
        args.log.error('Coordinates need to be in format XXhXXmXX.Xs and ' 
                       '+/-XXdXXmXX.Xs')
        sys.exit(1)
    Pos = []
    Xc, Yc = ([], [])
    for _scifits in SciFits_List:
        SciSpectra = _scifits[0].data
        y = np.median(SciSpectra[:, 410:440], axis=1)
        xc, yc = find_centroid(pos, y)
        Xc.append(xc)
        Yc.append(yc)
        args.log.info('%s has object centroid at: %0.2f, %0.2f' %
                      (_scifits.filename(), xc, yc))
        A = Astrometry(S.ra.deg, S.dec.deg, _scifits[0].header['PARANGLE'],
                       xc, yc, x_scale=1., y_scale=1., kind='lrs2')
        raoff, decoff = get_ADR_RAdec(xoff+xc, yoff+yc, A)
        ra, dec = A.tp.wcs_pix2world(pos[:, 0], pos[:, 1], 1)
        delta_ra = np.cos(np.deg2rad(S.dec.deg)) * (ra - S.ra.deg) * 3600.
        delta_dec = (dec - S.dec.deg) * 3600.
        ran1 = [np.min(delta_ra), np.max(delta_ra), np.min(delta_dec),
                np.max(delta_dec)]
        ran_list.append(ran1)
        Pos.append([delta_ra, delta_dec, raoff, decoff])
    ran_array = np.array(ran_list)
    rmax = np.max(ran_array, axis=0)
    rmin = np.min(ran_array, axis=0)
    ran = [np.floor(rmin[0]/scale)*scale, np.ceil(rmax[1]/scale)*scale,
           np.floor(rmin[2]/scale)*scale, np.ceil(rmax[3]/scale)*scale]   
    args.log.info('Cube limits - x: [%0.2f, %0.2f], y: [%0.2f, %0.2f]' %
                  (ran[0], ran[1], ran[2], ran[3]))
    
    ############################### Science ###################################
    cnt = 0
    F, info = get_cube(SciFits_List, Pos, scale, ran, sky, wave, cnt,
                       sky_subtract=True)
    xgrid = info[0][2]
    ygrid = info[0][3]
    norms = np.array([get_norm(i[0], xgrid, ygrid, wave) for i in info])
    norms = norms / np.mean(norms)
    for j, norm in enumerate(norms):
        args.log.info('Normalization for frame %i: %0.2f' % (j, norm))
    zcube = np.nanmean(np.array([i[0]/norm for i, norm in zip(info, norms)]), axis=0)
    ecube = np.sqrt(np.nanmean(np.array([i[1] for i in info])**2, axis=0))
    Header = _scifits[0].header
    outname = '%s_%s_cube.fits' % (args.galaxyname,  channel)
    eoutname = '%s_%s_error_cube.fits' % (args.galaxyname,  channel)
    write_cube(wave, xgrid, ygrid, zcube, outname, Header)
    write_cube(wave, xgrid, ygrid, ecube, eoutname, Header)
    
    ################################# Sky #####################################
    cnt = 1
    F1, info = get_cube(SkyFits_List, Pos, scale, ran, sky, wave, cnt, cor=cor,
                        sky_subtract=False)
    if len(info):
        xgrid = info[0][2]
        ygrid = info[0][3]
        zcube = np.nanmean(np.array([i[0] for i in info]), axis=0)
        ecube = np.sqrt(np.nanmean(np.array([i[1] for i in info])**2, axis=0))
        Header = _scifits[0].header
        outname = '%s_%s_sky_cube.fits' % (args.galaxyname,  channel)
        eoutname = '%s_%s_sky_error_cube.fits' % (args.galaxyname,  channel)
        write_cube(wave, xgrid, ygrid, zcube, outname, Header)
        write_cube(wave, xgrid, ygrid, ecube, eoutname, Header)
    outname = '%s_%s_multi.fits' % (args.galaxyname,  channel)
    fits.HDUList(F+F1).writeto(outname, overwrite=True)
    
    ################################# Cal #####################################
    cnt = 1
    F2, info = get_cube(CalFits_List, Pos, scale, ran, sky, wave, cnt,
                        sky_subtract=False, cal=True)
    if len(info):
        xgrid = info[0][2]
        ygrid = info[0][3]
        zcube = np.nanmean(np.array([i[0] for i in info]), axis=0)
        ecube = np.sqrt(np.nanmean(np.array([i[1] for i in info])**2, axis=0))
        Header = _scifits[0].header
        outname = '%s_%s_lamp_cube.fits' % (args.galaxyname,  channel)
        eoutname = '%s_%s_lamp_error_cube.fits' % (args.galaxyname,  channel)
        write_cube(wave, xgrid, ygrid, zcube, outname, Header)
        write_cube(wave, xgrid, ygrid, ecube, eoutname, Header)
    outname = '%s_%s_multi.fits' % (args.galaxyname,  channel)
    
    fits.HDUList(F+F1+F2).writeto(outname, overwrite=True)

main()
