# -*- coding: utf-8 -*-
"""
Created on Fri Oct  5 14:53:10 2018

@author: gregz
"""

import argparse as ap
import numpy as np
import os.path as op
import sys
import warnings

from astropy.convolution import Gaussian1DKernel, convolve
from astropy.convolution import interpolate_replace_nans
from astropy.modeling.models import Polynomial2D
from astropy.modeling.fitting import LevMarLSQFitter
from astropy.io import fits
from astropy.stats import biweight_midvariance, sigma_clipped_stats
from astropy.stats import sigma_clip
from astropy.table import Table
from input_utils import setup_logging
from scipy.interpolate import LSQBivariateSpline, interp1d
from scipy.signal import medfilt, savgol_filter
from skimage.feature import register_translation
from sklearn.decomposition import PCA
from utils import biweight_location


warnings.filterwarnings("ignore")


def get_script_path():
    return op.dirname(op.realpath(sys.argv[0]))

DIRNAME = get_script_path()

parser = ap.ArgumentParser(add_help=True)

parser.add_argument("-d", "--directory",
                    help='''base directory for reductions''',
                    type=str, default="")
parser.add_argument("sciobs",
                    help='''e.g., multi_20170126_0000011_exp02_farred.fits''',
                    type=str)
parser.add_argument("skyobs",
                    help='''e.g., multi_20170126_0000011_exp01_farred.fits''',
                    type=str, default=None)
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



def write_cube(wave,  xgrid, ygrid, zgrid, outname):
    hdu = fits.PrimaryHDU(np.array(zgrid, dtype='float32'))
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
    hdu.writeto(outname, overwrite=True)


def create_image_header(wave, xgrid, ygrid, zgrid, func=fits.ImageHDU):
    hdu = func(np.array(zgrid, dtype='float32'))
    hdu.header['CRVAL1'] = xgrid[0, 0]
    hdu.header['CRVAL2'] = ygrid[0, 0]
    hdu.header['CRPIX1'] = 1
    hdu.header['CRPIX2'] = 1
    hdu.header['CTYPE1'] = 'pixel'
    hdu.header['CTYPE2'] = 'pixel'
    hdu.header['CDELT1'] = xgrid[0, 1] - xgrid[0, 0]
    hdu.header['CDELT2'] = ygrid[1, 0] - ygrid[0, 0]
    return hdu


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


def correct_amplifier_offsets(data):
    y = np.mean(data[:, 400:-400], axis=1)
    x = np.arange(len(y))
    maskL = sigma_clip(y[:140], masked=True, maxiters=None)
    maskR = sigma_clip(y[140:], masked=True, maxiters=None)
    modelL = np.polyval(np.polyfit(x[:140][~maskL.mask], y[:140][~maskL.mask], 1), x[:140])
    modelR = np.polyval(np.polyfit(x[140:][~maskR.mask], y[140:][~maskR.mask], 1), x[140:])
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
    mask = sigma_clip(sky - cont, masked=True, maxiters=None)
    for i in np.arange(5):
        nsky = sky * 1.
        mask.mask[1:] += mask.mask[:-1]
        mask.mask[:-1] += mask.mask[1:]
        nsky[mask.mask] = np.nan
        cont = convolve(nsky, G)
        while np.isnan(cont).sum():
            cont = interpolate_replace_nans(cont, G)
        mask = sigma_clip(sky - cont, masked=True, maxiters=None)
    return mask.mask, cont    

def correct_wavelength_to_sky(spectra, skylines, wave, thresh=3.):
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

def main():
    F = fits.open(op.join(args.directory, args.sciobs))
    args.log.info('Science observation: %s loaded' % (args.sciobs))
    skyflag = args.skyobs.lower() != 'none'
    if skyflag:
        S = fits.open(op.join(args.directory, args.skyobs))
        args.log.info('Sky observation: %s loaded' % (args.sciobs))
    else:
        S = None
    
    # Basic Info Dump
    wave = F[6].data[0]
    channel = args.sciobs.split('_')[-1][:-5]
    channel_dict = {'uv': 'BL', 'orange': 'BR', 'red': 'RL', 'farred': 'RR'}
    specinit = channel_dict[channel]
    darfile = op.join(DIRNAME, 'lrs2_config/dar_%s.dat' % specinit)
    skylinefile = op.join(DIRNAME, 'lrs2_config/%s_skylines.dat' % channel)
    T = Table.read(darfile, format='ascii.fixed_width_two_line')
    SkyLines = Table.read(skylinefile, format='ascii.fixed_width_two_line')
    
    SciSpectra = F[0].data
    if skyflag:
        SkySpectra = S[0].data
    if not args.dont_correct_wavelength_to_sky:
        CW = correct_wavelength_to_sky(SciSpectra, SkyLines, wave)
    
    fits.PrimaryHDU(CW).writeto('test.fits', overwrite=True)
#    x, y = (F[5].data[:, 0], F[5].data[:, 1])
#    wave_0 = np.mean(wave)
#    xoff = (np.interp(wave, T['wave'], T['x_0']) -
#            np.interp(wave_0, T['wave'], T['x_0']))
#    yoff = (np.interp(wave, T['wave'], T['y_0']) -
#            np.interp(wave_0, T['wave'], T['y_0']))

#    # Initial Models
#    xn, yn = (0., 0.)
#    sel = np.where(((x - xn)**2 + (y-yn)**2) > 5.0**2)[0]
#    v = biweight_location(sci_list[1][sel, :] / sky_list[1][sel, :],
#                          axis=(0,))
#    gal_image = biweight_location(sci_list[1] - v * sky_list[1], axis=(1,))
#    loc = np.argmax(gal_image)
#    args.log.info('Peak found at %0.2f, %0.2f' % (x[loc], y[loc]))
#    xn, yn = (x[loc], y[loc])
#    d = (x - xn)**2 + (y-yn)**2
#    thresh = np.percentile(d, 90)
#    sel = np.where(((x - xn)**2 + (y-yn)**2) > thresh)[0]
#    v = biweight_location(sci_list[1][sel, :] / sky_list[1][sel, :],
#                          axis=(0,))
#        XN = biweight_location(sci_list[1], axis=(1,))
#        YN = biweight_location(v*sky_list[1], axis=(1,))
#        data = XN / YN
#        mean, median, std = sigma_clipped_stats(info[0]/info[1])
#        sel = np.abs(data - median) < 3. * std
#        P = Polynomial2D(2)
#        fitter = LevMarLSQFitter()
#        fit = fitter(P, x[sel], y[sel], data[sel])
#        offset = fit(x, y)
#        sky_list[1] = sky_list[1] * offset[:, np.newaxis]
#        gal_image = biweight_location(sci_list[1] - v * sky_list[1], axis=(1,))
#        sky, temp, norm1, norm2 = solve_system(sci_list, sky_list, x, y, xoff,
#                                               yoff, gal_image)
#        skysub = sci_list[1] - sky
#        for S, name in zip([sky, skysub], ['sky', 'skysub']):
#            outname = '%s_%s_%s_%s_%s_cube.fits' % (scidate, sciobs, sciexp,
#                                                    specname, name)
#            zcube, zimage, xgrid, ygrid = make_frame(x, y, S, wave, T['wave'],
#                                                     xoff, yoff,
#                                                     wstart=wave_0-50.,
#                                                     wend=wave_0+50.)
#            write_cube(wave, xgrid, ygrid, zcube, outname)
#        outname = '%s_%s_%s_%s_%s.fits' % ('multi', scidate, sciobs, sciexp,
#                                           specname)
#        X = np.array([T['wave'], T['x_0'], T['y_0']])
#        f1 = create_header_objection(wave, sci_list[1], func=fits.PrimaryHDU)
#        f2 = create_header_objection(wave, sky)
#        f3 = create_header_objection(wave, skysub)
#        fits.HDUList([f1, f2, f3, fits.ImageHDU(sci_list[0]),
#                      fits.ImageHDU(wave), fits.ImageHDU(zimage),
#                      fits.ImageHDU(X)]).writeto(outname, overwrite=True)

main()
