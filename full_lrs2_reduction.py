# -*- coding: utf-8 -*-
"""
Created on Thu Oct  4 13:30:06 2018

@author: gregz
"""
import numpy as np
import os.path as op
import os
import glob
import sys
import warnings
import tarfile
import argparse as ap
import requests
import uuid

from datetime import datetime
from astropy.io.votable import parse_single_table

from astropy.io import fits
from astropy.table import Table
from utils import biweight_location
from scipy.signal import savgol_filter, medfilt2d
from scipy.ndimage.filters import percentile_filter
from scipy.interpolate import interp1d, interp2d, griddata
from input_utils import setup_logging
from astrometry import Astrometry
from astropy.stats import biweight_midvariance
from photutils import DAOStarFinder
from astropy.modeling.models import Moffat2D
from astropy.convolution import Gaussian2DKernel, interpolate_replace_nans
from sklearn.decomposition import PCA

#try:
#    from pyhetdex.het.telescope import HetpupilModel
#    hetpupil_installed = True
#except ImportError:
#    print('Cannot find HETpupilModel.  Please check pyhetdex installation.')
#    print('For now, using default 50m**2 for mirror illumination')
#    hetpupil_installed = False

parser = ap.ArgumentParser(add_help=True)

parser.add_argument("-d", "--date",
                    help='''Date for reduction''',
                    type=str, default='20181108')

parser.add_argument("-o", "--object",
                    help='''Object name, no input reduces all objects''',
                    type=str, default=None)

args = parser.parse_args(args=None)

blueinfo = [['BL', 'uv', '503_056_7001', [3640., 4645.], ['LL', 'LU'],
             [4350., 4375.], ['hg_b', 'cd-a_b', 'fear_r', 'cd_b']],
            ['BR', 'orange', '503_056_7001',
             [4635., 6950.], ['RU', 'RL'], [6270., 6470.],
             ['hg_b', 'cd-a_b', 'fear_r', 'cd_b']]]
redinfo = [['RL', 'red', '502_066_7002', [6450., 8400.], ['LL', 'LU'],
            [7225., 7425.], ['hg_r', 'cd-a_b', 'fear_r', 'cd_b']],
           ['RR', 'farred', '502_066_7002',
            [8275., 10500.], ['RU', 'RL'], [9280., 9530.],
            ['hg_r', 'cd-a_b', 'fear_r', 'cd_b']]]

fplane_file = '/work/03730/gregz/maverick/fplane.txt'
twi_date = args.date
sci_date = args.date

# FOR LRS2
instrument = 'lrs2'

dither_pattern = np.zeros((50, 2))

log = setup_logging('panacea_quicklook')

baseraw = '/work/03946/hetdex/maverick'


sci_path = op.join(baseraw, sci_date,  '%s', '%s%s', 'exp%s',
                   '%s', '2*_%sLL*sci.fits')
cmp_path = op.join(baseraw, twi_date,  '%s', '%s%s', 'exp*',
                   '%s', '2*_%sLL_cmp.fits')
twi_path = op.join(baseraw, twi_date,  '%s', '%s%s', 'exp*',
                   '%s', '2*_%sLL_twi.fits')
bias_path = op.join(baseraw, twi_date, '%s', '%s%s', 'exp*',
                    '%s', '2*_%sLL_zro.fits')


def get_script_path():
    return op.dirname(op.realpath(sys.argv[0]))


def get_ifucenfile(side, amp,
                   virusconfig='/work/03946/hetdex/maverick/virus_config',
                   skiprows=4):

    file_dict = {"uv": "LRS2_B_UV_mapping.txt",
                 "orange": "LRS2_B_OR_mapping.txt",
                 "red": "LRS2_R_NR_mapping.txt",
                 "farred": "LRS2_R_FR_mapping.txt"}

    ifucen = np.loadtxt(op.join(virusconfig, 'IFUcen_files',
                        file_dict[side]), usecols=[0, 1, 2], skiprows=skiprows)

    if amp == "LL":
        return ifucen[140:, 1:3][::-1, :]
    if amp == "LU":
        return ifucen[:140, 1:3][::-1, :]
    if amp == "RL":
        return ifucen[:140, 1:3][::-1, :]
    if amp == "RU":
        return ifucen[140:, 1:3][::-1, :]


def orient_image(image, amp, ampname):
    '''
    Orient the images from blue to red (left to right)
    Fibers are oriented to match configuration files
    '''
    if amp == "LU":
        image[:] = image[::-1, ::-1]
    if amp == "RL":
        image[:] = image[::-1, ::-1]
    if ampname is not None:
        if ampname == 'LR' or ampname == 'UL':
            image[:] = image[:, ::-1]
    return image


def make_avg_spec(wave, spec, binsize=35, per=50):
    ind = np.argsort(wave.ravel())
    T = 1
    for p in wave.shape:
        T *= p
    wchunks = np.array_split(wave.ravel()[ind],
                             T / binsize)
    schunks = np.array_split(spec.ravel()[ind],
                             T / binsize)
    nwave = np.array([np.mean(chunk) for chunk in wchunks])
    nspec = np.array([np.percentile(chunk, per) for chunk in schunks])
    nwave, nind = np.unique(nwave, return_index=True)
    return nwave, nspec[nind]


def base_reduction(filename):
    a = fits.open(filename)
    image = np.array(a[0].data, dtype=float)
    # overscan sub
    overscan_length = 32 * (image.shape[1] / 1064)
    O = biweight_location(image[:, -(overscan_length-2):])
    image[:] = image - O
    # trim image
    image = image[:, :-overscan_length]
    gain = a[0].header['GAIN']
    gain = np.where(gain > 0., gain, 0.85)
    rdnoise = a[0].header['RDNOISE']
    rdnoise = np.where(rdnoise > 0., rdnoise, 3.)
    amp = (a[0].header['CCDPOS'].replace(' ', '') +
           a[0].header['CCDHALF'].replace(' ', ''))
    try:
        ampname = a[0].header['AMPNAME']
    except:
        ampname = None
    a = orient_image(image, amp, ampname) * gain
    E = np.sqrt(rdnoise**2 + np.where(a > 0., a, 0.))
    return a, E


def power_law(x, c1, c2=.5, c3=.15, c4=1., sig=2.5):
        return c1 / (c2 + c3 * np.power(abs(x / sig), c4))


def get_powerlaw(image, trace, spec):
    YM, XM = np.indices(image.shape)
    inds = []
    for j in np.arange(trace.shape[0]):
        inds.append(np.where(np.abs(YM - trace[j, np.newaxis, :]).ravel() <
                             7.)[0])
    inds = np.hstack(inds)
    inds = np.unique(inds)
    inds = np.setdiff1d(np.arange(image.shape[0] * image.shape[1]), inds)
    y, x = np.unravel_index(inds, image.shape)
    xlim = [0, image.shape[1]]
    xy = np.median(spec, axis=1)
    bottom = np.percentile(xy, 15)
    sel = np.where(xy > (15. * bottom))[0]
    log.info('Number of fibers that need powerlaw modeling: %i' % len(sel))
    xp = np.hstack([np.arange(xlim[0], xlim[1], 128), image.shape[1]-1])
    ylim = [0, image.shape[0]]
    yz = np.hstack([np.arange(ylim[0], ylim[1], 128), image.shape[0]-1])
    plaw, XX, YY = ([], [], [])
    YM, XM = np.indices(trace.shape)
    for xi in xp:
        for yi in yz:
            d = np.sqrt((yi - trace)**2 + (xi - XM)**2)
            plaw.append(np.nansum(spec * power_law(d, 1.4e-5, c3=2., c4=1.0,
                                                   sig=1.5)))
            XX.append(xi)
            YY.append(yi)
    for xi in xp:
        for s in sel:
            y0 = int(trace[s, xi])
            for i in np.arange(-4, 6, 2):
                d = np.sqrt(((y0+i) - trace)**2 + (xi - XM)**2)
                plaw.append(np.nansum(spec * power_law(d, 1.4e-5, c3=2.,
                                                       c4=1.0, sig=1.5)))
                YY.append(y0+i)
                XX.append(xi)
    plaw, XX, YY = [np.hstack(j) for j in [plaw, XX, YY]]
    grid_x, grid_y = np.meshgrid(np.arange(image.shape[1]),
                                 np.arange(image.shape[0]))
    C = griddata(np.array([XX, YY]).swapaxes(0, 1), plaw, (grid_x, grid_y),
                 method='cubic', fill_value=0.0)
    norm = np.zeros((image.shape[1],))
    for b in np.arange(image.shape[1]):
        sel = (x == b)
        norm[b] = np.median(image[y[sel], x[sel]] / C[y[sel], x[sel]])
    return C * savgol_filter(norm, 141, 1)[np.newaxis, :], norm


def get_bigW(amp, array_wave, array_trace, image):
    bigW = np.zeros(image.shape)
    Y, X = np.indices(array_wave.shape)
    YY, XX = np.indices(image.shape)
    for x, at, aw, xx, yy in zip(np.array_split(X, 2, axis=0),
                                 np.array_split(array_trace, 2, axis=0),
                                 np.array_split(array_wave, 2, axis=0),
                                 np.array_split(XX, 2, axis=0),
                                 np.array_split(YY, 2, axis=0)):
        for j in np.arange(at.shape[1]):
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                p0 = np.polyfit(at[:, j], aw[:, j], 7)
            bigW[yy[:, j], j] = np.polyval(p0, yy[:, j])
    return bigW


def get_bigF(array_trace, image):
    bigF = np.zeros(image.shape)
    Y, X = np.indices(array_trace.shape)
    YY, XX = np.indices(image.shape)
    n, m = array_trace.shape
    F0 = array_trace[:, m/2]
    for j in np.arange(image.shape[1]):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            p0 = np.polyfit(array_trace[:, j], F0, 7)
        bigF[:, j] = np.polyval(p0, YY[:, j])
    return bigF


def get_twiflat_field(flt_path, amps, array_wave, array_trace, bigW,
                      common_wave, masterbias, specname):
    files1 = sorted(glob.glob(flt_path.replace('LL', amps[0])))
    files2 = sorted(glob.glob(flt_path.replace('LL', amps[1])))

    array_list = []
    for filename1, filename2 in zip(files1, files2):
        log.info('Prepping flat %s' % filename1)
        array_flt1, e1 = base_reduction(filename1)
        array_flt2, e2 = base_reduction(filename2)
        array_flt = np.vstack([array_flt1, array_flt2])
        array_flt[:] -= masterbias
        array_list.append(array_flt)
    array_list = np.array(array_list)
    if len(array_list) > 1:
        norm = np.median(array_list, axis=(1,2))
        array_flt = np.median(array_list / norm[:, np.newaxis, np.newaxis],
                              axis=0)
    else:
        array_flt = np.squeeze(np.array(array_list))
        array_flt[:] /= np.median(array_flt)
    
    log.info('Working on flat')
    x = np.arange(array_wave.shape[1])
    spectrum = array_trace * 0.
    for fiber in np.arange(array_wave.shape[0]):
        indl = np.floor(array_trace[fiber]).astype(int)
        indh = np.ceil(array_trace[fiber]).astype(int)
        spectrum[fiber] = array_flt[indl, x] / 2. + array_flt[indh, x] / 2.
    
    log.info('Getting powerlaw for side %s' % specname)
    plaw, norm = get_powerlaw(array_flt, array_trace, spectrum)
    array_flt[:] -= plaw
    array_flt[:] = np.where(array_flt < 0., 0., array_flt)
    fits.PrimaryHDU(plaw).writeto('test_plaw_%s.fits' % specname, overwrite=True)
    fits.PrimaryHDU(array_flt).writeto('test_spec_%s.fits' % specname, overwrite=True)
    smooth = savgol_filter(spectrum, 315, 1, axis=1)
    avg = biweight_location(smooth, axis=(0,))
    norm = biweight_location(smooth / avg, axis=(1,))
    nw, ns = make_avg_spec(array_wave, spectrum / norm[:, np.newaxis],
                           binsize=41, per=50)
    I = interp1d(nw, ns, kind='linear', fill_value='extrapolate')
    ftf = spectrum * 0.
    for fiber in np.arange(array_wave.shape[0]):
        model = I(array_wave[fiber])
        ftf[fiber] = savgol_filter(spectrum[fiber] / model, 151, 1)
    nw1, ns1 = make_avg_spec(array_wave, spectrum / ftf, binsize=41,
                             per=95)
    I = interp1d(nw1, ns1, kind='quadratic', fill_value='extrapolate')
    modelimage = I(bigW)
    flat = array_flt / modelimage
    flat[~np.isfinite(flat)] = 0.0
    flat[flat < 0.0] = 0.0
    return flat


def get_spectra(array_flt, array_trace):
    spectrum = array_trace * 0.
    x = np.arange(array_flt.shape[1])
    for fiber in np.arange(array_trace.shape[0]):
        indl = np.floor(array_trace[fiber]).astype(int)
        indh = np.ceil(array_trace[fiber]).astype(int)
        spectrum[fiber] = array_flt[indl, x] / 2. + array_flt[indh, x] / 2.
    return spectrum


def safe_division(num, denom, eps=1e-8, fillval=0.0):
    good = np.isfinite(denom) * (np.abs(denom) > eps)
    div = num * 0.
    if num.ndim == denom.ndim:
        div[good] = num[good] / denom[good]
        div[~good] = fillval
    else:
        div[:, good] = num[:, good] / denom[good]
        div[:, ~good] = fillval
    return div


def find_cosmics(Y, E, trace, thresh=8., ran=0):
    x = np.arange(trace.shape[1])
    C = Y * 0.
    for fiber in np.arange(trace.shape[0]):
        indl = np.floor(trace[fiber]).astype(int)
        T = np.zeros((4, trace.shape[1], 4))
        flag = True
        for ss, k in enumerate(np.arange(-1, 3)):
            try:
                T[0, :, ss] = Y[indl+k, x]
                T[1, :, ss] = E[indl+k, x]
                T[2, :, ss] = indl+k
                T[3, :, ss] = x
            except:
                flag = False
        if flag:
            m = np.median(T[0], axis=1)
            P = np.abs(T[0] - m[:, np.newaxis]) / T[1]
            C[T[2][P > thresh].astype(int), T[3][P > thresh].astype(int)] = 1.
    C = np.array(C, dtype=bool)
    log.info('Number of fiber pixels hit by cosmics: %i' % C.sum())
    return C


def weighted_extraction(image, error, flat, trace):
    E = safe_division(error, flat)
    E[E < 1e-8] = 1e9
    Y = safe_division(image, flat)
    nY = Y * 1.
    C = np.array(Y * 0., dtype=bool)
    for i in np.arange(1):
        cosmics = find_cosmics(nY, E, trace, ran=1)
        C = C + cosmics

    x = np.arange(trace.shape[1])
    spectrum = 0. * trace
    error_spec = 0. * trace
    Fimage = image * 0.
    for fiber in np.arange(trace.shape[0]):
        T = np.zeros((4, trace.shape[1], 4))
        indl = np.floor(trace[fiber]).astype(int)
        flag = False
        for ss, k in enumerate(np.arange(0, 2)):
            try:
                T[0, :, ss] = Y[indl+k, x]
                T[1, :, ss] = 1. / E[indl+k, x]**2
                T[2, :, ss] = ~C[indl+k, x]
                T[3, :, ss] = error[indl+k, x]**2
                Fimage[indl+k, x] = fiber + 1
            except:
                v = indl+k
                sel = np.where((v >= 0) * (v < Y.shape[0]))[0]
                T[0, sel, ss] = Y[v[sel], x[sel]]
                T[1, sel, ss] = 1. / E[v[sel], x[sel]]**2
                T[2, sel, ss] = ~C[v[sel], x[sel]]
                T[3, sel, ss] = error[v[sel], x[sel]]**2
                Fimage[v[sel], x[sel]] = fiber + 1
                flag = True
        if flag:
            if np.mean(indl) > (Y.shape[0]/2.):
                k = 2
            else:
                k = -1
            v = indl+k
            sel = np.where((v >= 0) * (v < len(x)))[0]
            a = np.sum(T[0, sel] * T[1, sel] * T[2, sel], axis=1)
            b = np.sum(T[1, sel] * T[2, sel], axis=1)
            spectrum[fiber, sel] = safe_division(a, b)
            sel1 = T[2, sel].sum(axis=1) < 2.
            spectrum[fiber][sel][sel1] = 0.0
            error_spec[fiber][sel][sel1] = 0.0
        else:
            a = np.sum(T[0] * T[1] * T[2], axis=1)
            b = np.sum(T[1] * T[2], axis=1)
            spectrum[fiber] = safe_division(a, b)
            a = np.sum(T[3] * T[1] * T[2], axis=1)
            b = np.sum(T[1] * T[2], axis=1)
            error_spec[fiber] = np.sqrt(safe_division(a, b))
            sel = T[2].sum(axis=1) < 2.
            spectrum[fiber][sel] = 0.0
            error_spec[fiber][sel] = 0.0
    return spectrum, error_spec, C, Y, Fimage


def get_trace_shift(sci_array, flat, array_trace, Yx):
    YM, XM = np.indices(flat.shape)
    inds = np.zeros((3, array_trace.shape[0], array_trace.shape[1]))
    XN = np.round(array_trace)
    inds[0] = XN - 1.
    inds[1] = XN + 0.
    inds[2] = XN + 1.
    inds = np.array(inds, dtype=int)
    Trace = array_trace * 0.
    FlatTrace = array_trace * 0.
    N = YM.max()
    x = np.arange(array_trace.shape[1])
    for i in np.arange(Trace.shape[0]):
        sel = YM[inds[0, i, :], x] >= 0.
        sel = sel * (YM[inds[2, i, :], x] < N)
        xmax = (YM[inds[1, i, sel], x[sel]] -
                (sci_array[inds[2, i, sel], x[sel]] -
                sci_array[inds[0, i, sel], x[sel]]) /
                (2. * (sci_array[inds[2, i, sel], x[sel]] -
                 2. * sci_array[inds[1, i, sel], x[sel]] +
                 sci_array[inds[0, i, sel], x[sel]])))
        Trace[i, sel] = xmax
        xmax = (YM[inds[1, i, sel], x[sel]] - (flat[inds[2, i, sel], x[sel]] -
                flat[inds[0, i, sel], x[sel]]) /
                (2. * (flat[inds[2, i, sel], x[sel]] - 2. *
                 flat[inds[1, i, sel], x[sel]] +
                 flat[inds[0, i, sel], x[sel]])))
        FlatTrace[i, sel] = xmax
    shifts = np.nanmedian(FlatTrace - Trace, axis=1)
    shifts = np.polyval(np.polyfit(np.nanmedian(FlatTrace, axis=1), shifts, 2),
                        Yx)
    return shifts


def modify_spectrum(spectrum, w, xloc, yloc):
    for i in np.arange(spectrum.shape[0]):
        sel = spectrum[i] == 0.
        I = interp1d(w[i][~sel], spectrum[i][~sel], kind='quadratic',
                     fill_value='extrapolate')
        dw = np.diff(w[i])
        dw = np.hstack([dw[0], dw])
        spectrum[i] = I(w[i]) / dw
    return spectrum
#    norm = np.zeros((spectrum.shape[0],))
#    bins = np.arange(w.min(), w.max(), 40)
#    x = bins + 20.
#    Z = np.zeros((spectrum.shape[0], len(bins)))
#    for j, bini in enumerate(bins):
#        for i in np.arange(spectrum.shape[0]):
#            sel = np.where((w[i] > bini) * (w[i] < bini + 100))[0]
#            norm[i] = np.mean(spectrum[i][sel])
#        smooth = norm / np.mean(norm)
#        Z[:, j] = smooth
#    ftf = spectrum * 0.
#    for i in np.arange(spectrum.shape[0]):
#        I = interp1d(x, Z[i, :], kind='quadratic', fill_value='extrapolate')
#        ftf[i] = I(w[i])
#    return spectrum / ftf


def extract_sci(sci_path, amps, flat, array_trace, array_wave, bigW,
                masterbias, pos):
    files1 = sorted(glob.glob(sci_path.replace('LL', amps[0])))
    files2 = sorted(glob.glob(sci_path.replace('LL', amps[1])))
    
    xloc, yloc = (pos[:, 0], pos[:, 1])
    array_list = []
    for filename1, filename2 in zip(files1, files2):
        log.info('Prepping sci %s' % filename1)
        array_flt1, e1 = base_reduction(filename1)
        array_flt2, e2 = base_reduction(filename2)
        array_flt = np.vstack([array_flt1, array_flt2])
        array_flt[:] -= masterbias
        array_list.append(array_flt)
    if len(array_list) > 1:
        sci_array = np.sum(array_list, axis=0)
    else:
        sci_array = np.squeeze(np.array(array_list))
    Xx = np.arange(flat.shape[1])
    Yx = np.arange(flat.shape[0])
    I = interp2d(Xx, Yx, flat, kind='cubic', bounds_error=False,
                 fill_value=0.0)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        shifts = get_trace_shift(sci_array, flat, array_trace, Yx)
    flat = I(Xx, Yx + shifts)
    log.info('Found shift for %s of %0.3f' % (files1[0], np.median(shifts)))
    array_list = []
    spec_list, error_list = ([], [])
    orig_list = []
    clist, flist, Flist = ([], [], [])
    for filename1, filename2 in zip(files1, files2):
        log.info('Fiber extraction sci %s' % filename1)
        array_flt1, e1 = base_reduction(filename1)
        array_flt2, e2 = base_reduction(filename2)
        array_flt = np.vstack([array_flt1, array_flt2])
        array_err = np.vstack([e1, e2])

        array_flt[:] -= masterbias
        log.info('Getting powerlaw for %s' % filename1)
        spectrum = get_spectra(array_flt, array_trace)
        plaw, norm = get_powerlaw(array_flt, array_trace, spectrum)
        array_flt[:] -= plaw
        array_list.append(array_flt)
        spectrum, error, c, fl, Fimage = weighted_extraction(array_flt,
                                                             array_err,
                                                             flat, array_trace)
        sel = ~np.isfinite(spectrum)
        spectrum[sel] = 0.0
        error[sel] = 0.0
        log.info('Number of 0.0 pixels in spectra: %i' %
                 (spectrum == 0.0).sum())
        speclist, errorlist = ([], [])
        spectrum = modify_spectrum(spectrum, array_wave, xloc, yloc)
        log.info('Number of 0.0 pixels in spectra: %i' %
                 (spectrum == 0.0).sum())
        for fiber in np.arange(array_wave.shape[0]):
            I = interp1d(array_wave[fiber], spectrum[fiber],
                         kind='linear', fill_value='extrapolate')
            coV = np.zeros((len(commonwave), spectrum.shape[1]))
            for i in np.arange(len(commonwave)):
                w = np.zeros((len(commonwave),))
                w[i] = 1.
                coV[i] = np.interp(array_wave[fiber], commonwave, w)
            error_interp = np.sqrt((coV * error[fiber]**2).sum(axis=1))
            speclist.append(I(commonwave))
            errorlist.append(error_interp)
        log.info('Finished rectifying spectra')
        spec_list.append(np.array(speclist))
        error_list.append(np.array(errorlist))
        orig_list.append(spectrum)
        clist.append(c)
        flist.append(fl)
        Flist.append(Fimage)
    return (np.array(array_list), np.array(spec_list), np.array(orig_list),
            np.array(clist, dtype=float), np.array(flist), np.array(Flist),
            np.array(error_list))


def get_masterbias(zro_path, amp):
    files = glob.glob(zro_path.replace('LL', amp))
    listzro = []
    for filename in files:
        a, error = base_reduction(filename)
        listzro.append(a)
    return np.median(listzro, axis=0)


def get_masterarc(arc_path, amp, arc_names, masterbias, specname):
    files = glob.glob(arc_path.replace('LL', amp))
    # if specname == 'farred':
    #     ofiles = glob.glob(arc_path.replace('LL', amp).replace('cmp', 'sci'))
    #     files = files + ofiles
    listarc, listarce = ([], [])
    for filename in files:
        f = fits.open(filename)
        # scicond = ('_056' in f[0].header['OBJECT']) and ('sci' in filename)
        if f[0].header['OBJECT'].lower() in arc_names:
            a, e = base_reduction(filename)
            a[:] -= masterbias
            listarc.append(a)
            listarce.append(e)
    listarc, listarce = [np.array(x) for x in [listarc, listarce]]
    # total_error = np.sqrt(np.sum(listarce**2, axis=0))
    return np.sum(listarc, axis=0)


def get_mastertwi(twi_path, amp, masterbias):
    files = glob.glob(twi_path.replace('LL', amp))
    listtwi = []
    for filename in files:
        a, e = base_reduction(filename)
        a[:] -= masterbias
        listtwi.append(a)
    twi_array = np.array(listtwi, dtype=float)
    norm = np.median(twi_array, axis=(1, 2))[:, np.newaxis, np.newaxis]
    return np.median(twi_array / norm, axis=0)


def get_trace_reference(specid, ifuslot, ifuid, amp, obsdate,
                        virusconfig='/work/03946/hetdex/'
                        'maverick/virus_config'):
    files = glob.glob(op.join(virusconfig, 'Fiber_Locations', '*',
                              'fiber_loc_%s_%s_%s_%s.txt' %
                              (specid, ifuslot, ifuid, amp)))
    dates = [op.basename(op.dirname(fn)) for fn in files]
    obsdate = datetime(int(obsdate[:4]), int(obsdate[4:6]),
                       int(obsdate[6:]))
    timediff = np.zeros((len(dates),))
    for i, datei in enumerate(dates):
        d = datetime(int(datei[:4]), int(datei[4:6]),
                     int(datei[6:]))
        timediff[i] = np.abs((obsdate - d).days)
    ref_file = np.loadtxt(files[np.argmin(timediff)])
    return ref_file


def get_trace(twilight, specid, ifuslot, ifuid, amp, obsdate):
    ref = get_trace_reference(specid, ifuslot, ifuid, amp, obsdate)
    N1 = (ref[:, 1] == 0.).sum()
    good = np.where(ref[:, 1] == 0.)[0]

    def get_trace_chunk(flat, XN):
        YM = np.arange(flat.shape[0])
        inds = np.zeros((3, len(XN)))
        inds[0] = XN - 1.
        inds[1] = XN + 0.
        inds[2] = XN + 1.
        inds = np.array(inds, dtype=int)
        Trace = (YM[inds[1]] - (flat[inds[2]] - flat[inds[0]]) /
                 (2. * (flat[inds[2]] - 2. * flat[inds[1]] + flat[inds[0]])))
        return Trace
    image = twilight
    N = 40
    xchunks = np.array([np.mean(x) for x in
                        np.array_split(np.arange(image.shape[1]), N)])
    chunks = np.array_split(image, N, axis=1)
    flats = [np.median(chunk, axis=1) for chunk in chunks]
    Trace = np.zeros((len(ref), len(chunks)))
    k = 0
    for flat, x in zip(flats, xchunks):
        diff_array = flat[1:] - flat[:-1]
        loc = np.where((diff_array[:-1] > 0.) * (diff_array[1:] < 0.))[0]
        peaks = flat[loc+1]
        loc = loc[peaks > 0.1 * np.median(peaks)]+1
        trace = get_trace_chunk(flat, loc)
        T = np.zeros((len(ref)))
        if len(trace) == N1:
            T[good] = trace
            for missing in np.where(ref[:, 1] == 1)[0]:
                gind = np.argmin(np.abs(missing - good))
                T[missing] = (T[good[gind]] + ref[missing, 0] -
                              ref[good[gind], 0])
        Trace[:, k] = T
        k += 1
    x = np.arange(twilight.shape[1])
    trace = np.zeros((Trace.shape[0], twilight.shape[1]))
    for i in np.arange(Trace.shape[0]):
        sel = Trace[i, :] > 0.
        trace[i] = np.polyval(np.polyfit(xchunks[sel], Trace[i, sel], 7), x)
    return trace, ref


def find_peaks(y, thresh=8.):
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
    return peak_loc, peaks/std, peaks


def robust_polyfit(x, y, order=3, niter=3):
    sel = y > 0.
    ymod = np.polyval(np.polyfit(x[sel], y[sel], order), x)
    for i in np.arange(niter):
        a = np.abs(y - ymod)
        mad = np.median(a)
        sel = a < 3. * mad
        if sel.sum() > (order+2):
            ymod = np.polyval(np.polyfit(x[sel], y[sel], order), x)
    return ymod


def count_matches(lines, loc, fib, cnt=5):
    found_lines = np.zeros((trace.shape[0], len(lines)))
    M = np.zeros((cnt, cnt))
    for k in np.arange(cnt):
        for j in np.arange(cnt):
            i1 = 0 + k
            i2 = -1 - j
            diff = [loc[fib][i1] - lines['col2'][0],
                    loc[fib][i2] - lines['col2'][-1]]
            m = (diff[1] - diff[0]) / (lines['col2'][-1] - lines['col2'][0])
            y = np.array(m * (lines['col2'] - lines['col2'][0]) +
                         diff[0] + lines['col2'])
            for i, line in enumerate(lines):
                col = y[i]
                v = np.abs(col - loc[fib])
                if np.min(v) < 3.:
                    found_lines[fib, i] = loc[fib][np.argmin(v)]
            M[k, j] = (found_lines[fib] > 0.).sum()
    return np.unravel_index(np.argmax(M), M.shape)


def get_wavelength_from_arc(image, trace, lines, side, amp):
    if side == 'uv':
        thresh = 3.  # 5
    if side == 'orange':
        thresh = 3.  # 8
    if side == 'red':
        thresh = 3.  # 50
    if side == 'farred':
        thresh = 3.  # 20
    spectrum = get_spectra(image, trace)
    fib = np.argmax(np.median(spectrum, axis=1))
    cont = percentile_filter(spectrum, 15, (1, 101))
    spectrum -= cont
    x = np.arange(trace.shape[1])
    loc = []
    ph, pr = ([], [])
    for i, spec in enumerate(spectrum):
        px, ps, py = find_peaks(spec, thresh=thresh)
        loc.append(px)
        ph.append(ps)
        pr.append(py)

    found_lines = np.zeros((trace.shape[0], len(lines)))
    ls = np.argsort(lines['col3'])[::-1]

    inds = np.argsort(pr[fib])[::-1]
    for ind in inds:
        off = loc[fib][ind] - lines['col2'][ls[0]]
        if np.abs(off) < 50.:
            break
    off = loc[fib][ind] - lines['col2'][ls[0]]
    found_lines[fib, ls[0]] = loc[fib][ind]
    y = lines['col2'] + off
    s = found_lines[fib] * 0.
    pp = s * 0.
    s[ls[0]] = 1.
    pp[ls[0]] = 0.
    for l in ls[1:]:
        guess = y[l]
        v = np.abs(guess - loc[fib])
        ER = lines['col3'][l] / lines['col3'][ls[0]]
        MR = pr[fib] / pr[fib][ind]
        EE = MR * np.sqrt(1./ph[fib]**2 + 1./ph[fib][ind])
        EE = np.max([EE, .1 * MR], axis=0)
        dist = v/2. + np.abs(ER - MR) / EE
        if np.min(dist) < 10.:
            ind1 = np.argmin(dist)
            found_lines[fib, l] = loc[fib][ind1]
            ll = np.where(found_lines[fib] > 0.)[0][0]
            lh = np.where(found_lines[fib] > 0.)[0][-1]
            diff = [found_lines[fib, ll] - lines['col2'][ll],
                    found_lines[fib, lh] - lines['col2'][lh]]
            m = ((diff[1] - diff[0]) /
                 (lines['col2'][lh] - lines['col2'][ll]))
            y = np.array(m * (lines['col2'] - lines['col2'][ll]) +
                         diff[0] + lines['col2'])
            s[l] = MR[ind1]
            pp[l] = dist[ind1]
    inds = np.where(found_lines[fib] > 0.)[0]
    for ind in inds:
        print(lines['col1'][ind], lines['col2'][ind], found_lines[fib][ind],
              lines['col3'][ind], s[ind], pp[ind])

    for i, line in enumerate(lines):
        if found_lines[fib, i] == 0.:
            continue
        for j in np.arange(0, fib)[::-1]:
            if len(loc[j]) < 1:
                continue
            k = j + 1
            v = found_lines[k, i]
            while (v == 0.) and (k < trace.shape[0]):
                k += 1
                v = found_lines[k, i]
            m = np.abs(loc[j] - v)
            if np.min(m) < 2.:
                found_lines[j, i] = loc[j][np.argmin(m)]
        for j in np.arange(fib+1, trace.shape[0]):
            if len(loc[j]) < 1:
                continue
            k = j - 1
            v = found_lines[k, i]
            while (v == 0.) and (k >= 0):
                k -= 1
                v = found_lines[k, i]
            m = np.abs(loc[j] - v)
            if np.min(m) < 2.:
                found_lines[j, i] = loc[j][np.argmin(m)]
    fits.PrimaryHDU(found_lines, 'fl_%s_%s.fits' % (side, amp), overwrite=True)
    for i, line in enumerate(lines):
        if np.sum(found_lines[:, i]) < (0.5 * trace.shape[0]):
            found_lines[:, i] = 0.0
            continue
        ind = np.array(found_lines[:, i], dtype=int)
        xt = trace[np.arange(trace.shape[0]), ind]
        yt = robust_polyfit(xt, found_lines[:, i])
        found_lines[:, i] = yt
    wave = trace * 0.0
    res = np.zeros((trace.shape[0],))
    for j in np.arange(trace.shape[0]):
        sel = found_lines[j, :] > 0.0
        if sel.sum() < 3:
            continue
        wave[j] = np.polyval(np.polyfit(found_lines[j, sel],
                             lines['col1'][sel], 3), x)
        res[j] = np.std(np.interp(found_lines[j, sel], x, wave[j]) -
                        lines['col1'][sel])
    missing = np.where(np.all(wave == 0., axis=1))[0]
    if len(missing):
        good = np.where(~np.all(wave == 0., axis=1))[0]
        for j in np.arange(trace.shape[1]):
            x = trace[good, j]
            y = wave[good, j]
            wave[missing, j] = np.polyval(np.polyfit(x, y, 3),
                                          trace[missing, j])
    log.info('Min, Max Wave: %0.2f, %0.2f' % (wave.min(), wave.max()))
    log.info('Mean Res, Median Res: %0.3f, %0.3f' % (np.mean(res),
                                                     np.median(res)))
    return wave


def get_objects(basefiles, attrs, full=False):
    s = []
    for fn in basefiles:
        F = fits.open(fn)
        s.append([])
        for att in attrs:
            s[-1].append(F[0].header[att])
        if full:
            area = get_mirror_illumination(fn)
            throughput = get_throughput(fn, s[-1][1])
            s[-1].append(area)
            s[-1].append(throughput)
    return s


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


def correct_ftf(rect, error):
    def outlier(y, y1, oi):
        m = np.abs(y[oi] - y1[oi])
        o = (y - y1) > 3. * np.median(m)
        return o

    def fit_sky_col(x, y):
        o = y == 0.
        low = np.percentile(y[~o], 16)
        mid = np.percentile(y[~o], 50)
        high = np.percentile(y[~o], 84)
        flag = False
        if (high - mid) > 2.0 * (mid - low):
            y1 = np.ones(x[~o].shape) * np.percentile(y[~o], 5)
            log.info('Object is too bright for fiber to fiber correction.')
            flag = True
        else:
            y1 = savgol_filter(y[~o], 31, 3)
        I = interp1d(x[~o], y1, kind='quadratic', fill_value='extrapolate')
        y1 = I(x)
        for i in np.arange(3):
            o += outlier(y, y1, ~o)
            y1 = savgol_filter(y[~o], 51, 3)
            I = interp1d(x[~o], y1, kind='quadratic', fill_value='extrapolate')
            y1 = I(x)
        return y1, o, flag

    x = np.arange(rect.shape[0])
    y = np.median(rect, axis=1)
    f, o, flag = fit_sky_col(x, y)
    ftf = f / np.median(f)
    if not flag:
        return rect / ftf[:, np.newaxis], error / ftf[:, np.newaxis]
    else:
        return rect, error

def sky_subtraction(rect, error, ncomponents=25):
    def outlier(y, y1, oi):
        m = np.abs(y[oi] - y1[oi])
        o = (y - y1) > 3. * np.median(m)
        return o

    def fit_sky_col(x, y):
        o = y == 0.
        low = np.percentile(y[~o], 16)
        mid = np.percentile(y[~o], 50)
        high = np.percentile(y[~o], 84)
        flag = False
        if (high - mid) > 2.0 * (mid - low):
            y1 = np.ones(x[~o].shape) * np.percentile(y[~o], 5)
            log.info('Object is too bright for regular sky subtraction.')
            log.info('Subtracting the 5th percentile')
            flag = True
        else:
            y1 = savgol_filter(y[~o], 31, 3)
        I = interp1d(x[~o], y1, kind='quadratic', fill_value='extrapolate')
        y1 = I(x)
        for i in np.arange(3):
            o += outlier(y, y1, ~o)
            y1 = savgol_filter(y[~o], 51, 3)
            I = interp1d(x[~o], y1, kind='quadratic', fill_value='extrapolate')
            y1 = I(x)
        return y1, o, flag

    x = np.arange(rect.shape[0])
    y = np.median(rect, axis=1)
    f, o, flag = fit_sky_col(x, y)
    o[:-1] += o[1:]
    o[1:] += o[:-1]
    if flag:
        sky = np.percentile(rect, 5, axis=0)[np.newaxis] * np.ones((280, 1))
        return sky
    md = np.median(rect[~o], axis=0)
    msub = rect[~o] - md
    pca = PCA(n_components=ncomponents)
    H = pca.fit_transform(msub.swapaxes(0, 1))
    res = rect * 0.
    for s in np.arange(rect.shape[0]):
        sol = np.linalg.lstsq(H, rect[s] - md)[0]
        res[s] = np.dot(H, sol)
    m = np.interp(x, x[~o], np.median(rect[~o] - md - res[~o], axis=1))
    sky = md + res + m[:, np.newaxis]
    return sky


def correct_fiber_to_fiber(data, xloc, yloc, seeing=1.5):
    D = np.sqrt((xloc[:, np.newaxis] - xloc)**2 +
                (yloc[:, np.newaxis] - yloc)**2)
    W = np.exp(-0.5 / (seeing/2.35)**2 * D**2)

    def outlier(y, y1, oi):
        m = np.abs(y[oi] - y1[oi])
        o = (y - y1) > 3. * np.median(m)
        return o

    o = data == 0.
    low = np.percentile(data[~o], 16)
    mid = np.percentile(data[~o], 50)
    high = np.percentile(data[~o], 84)
    if (high - mid) > 2.0 * (mid - low):
        log.info('Source is too bright to correct fiber to fiber')
        return np.ones(data.shape)
    smooth = (data[~o, np.newaxis] * W[~o]).sum(axis=0) / W[~o].sum(axis=0)
    o = outlier(data, smooth, data > 0.)
    for i in np.arange(3):
        smooth = (data[~o, np.newaxis] * W[~o]).sum(axis=0) / W[~o].sum(axis=0)
        o = outlier(data, smooth, ~o)
    return smooth


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
        zgrid[k, :, :] = ((data[sel, k][:, np.newaxis, np.newaxis] *
                           W[sel]).sum(axis=0) / W[sel].sum(axis=0) *
                          (scale**2 / area))
    wi = np.searchsorted(wave, wstart, side='left')
    we = np.searchsorted(wave, wend, side='right')

    zimage = biweight_location(zgrid[wi:we+1], axis=(0,))
    return zgrid, zimage, xgrid, ygrid


def write_cube(wave, xgrid, ygrid, zgrid, outname):
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


def find_source(image, xgrid, ygrid):
    std = np.sqrt(biweight_midvariance(image))
    daofind = DAOStarFinder(fwhm=4.0, threshold=5.*std)
    sources = daofind(image)
    print(sources)
    if len(sources) >= 1:
        ind = np.argmax(sources['flux'])
        xg = np.unique(xgrid)
        yg = np.unique(ygrid)
        xc = np.interp(sources[ind]['xcentroid'], np.arange(len(xg)), xg)
        yc = np.interp(sources[ind]['ycentroid'], np.arange(len(yg)), yg)
        return xc, yc
    else:
        return None


def extract_source(data, xc, yc, xoff, yoff, wave, xloc, yloc, error,
                   seeing=1.5):
    gamma = seeing / (np.sqrt(2**(1 / 3.5) - 1.) * 2.)
    PSF = Moffat2D(amplitude=1., x_0=0., y_0=0., alpha=3.5, gamma=gamma)
    spec = wave * 0.
    serror = wave * 0.
    for i in np.arange(len(wave)):
        x = xc + xoff[i]
        y = yc + yoff[i]
        PSF.x_0.value = x
        PSF.y_0.value = y
        W = PSF(xloc, yloc)
        W /= W.sum()
        spec[i] = (data[:, i] * W).sum() / (W**2).sum()
        serror[i] = np.sqrt((error[:, i]**2 * W).sum() / (W**2).sum())
    return spec, serror


def get_mirror_illumination(fn=None):
    ''' Use hetillum from illum_lib to calculate mirror illumination (cm^2) '''
    log.info('Getting mirror illumination')
    try:
        F = fits.open(fn)
        names = ['RHO_STRT', 'THE_STRT', 'PHI_STRT', 'X_STRT', 'Y_STRT']
        r, t, p, x, y = [F[0].header[name] for name in names]
        mirror_illum = float(os.popen('/home/00156/drory/illum_lib/hetillum -p'
                             ' -x "[%0.4f,%0.4f,%0.4f]" "[%0.4f,%0.4f]" 256' %
                                      (r, t, p, x, y)).read().split('\n')[0])
        area = mirror_illum * 51.4 * 1e4
    except:
        log.info('Using default mirror illumination value')
        area = 51.4 * 1e4
#    if hetpupil_installed:
#        log.info('Using HetpupilModel from pyhetdex')
#        try:
#            mirror_illum = HetpupilModel([fn], normalize=False)
#            area = mirror_illum.fill_factor[0] * 55. * 1e4
#        except:
#            log.info('Using default mirror illumination value')
#            area = 50. * 1e4
#    else:
#        log.info('Using default mirror illumination value')
#        area = 50. * 1e4
    log.info('Mirror illumination: %0.2f m^2' % (area/1e4))
    return area


def panstarrs_query(ra_deg, dec_deg, rad_deg, mindet=1,
                    maxsources=30000,
                    server=('https://archive.stsci.edu/panstarrs/search.php')):
    """
    Query Pan-STARRS DR1 @ MAST
    parameters: ra_deg, dec_deg, rad_deg: RA, Dec, field
                                          radius in degrees
                mindet: minimum number of detection (optional)
                maxsources: maximum number of sources
                server: servername
    returns: astropy.table object
    """
    r = requests.get(server, params={'RA': ra_deg, 'DEC': dec_deg,
                                     'SR': rad_deg, 'max_records': maxsources,
                                     'outputformat': 'VOTable',
                                     'ndetections': ('>%d' % mindet)})

    # write query data into local file
    name = str(uuid.uuid4()) + '.xml'
    outf = open(name, 'w')
    outf.write(r.text)
    outf.close()

    # parse local file into astropy.table object
    data = parse_single_table(name)
    os.remove(name)
    return data.to_table(use_names_over_ids=True)


def get_throughput(fn, exptime, path='/work/03946/hetdex/maverick'):
    return 1.
    attr = ['TARGTRA', 'TARGTDEC', 'GUIDLOOP', 'MJD', 'PSFMAG', 'STARMAGS']
    M = []
    path = op.join(path, args.date)
    f = op.basename(fn)
    DT = f.split('_')[0]
    y, m, d, h, mi, s = [int(x) for x in [DT[:4], DT[4:6], DT[6:8], DT[9:11],
                         DT[11:13], DT[13:15]]]
    d0 = datetime(y, m, d, h, mi, s)
    for gp in ['gc1', 'gc2']:
        tarfolder = op.join(path, gp, '%s.tar' % gp)
        T = tarfile.open(tarfolder, 'r')
        init_list = sorted([name for name in T.getnames()
                            if name[-5:] == '.fits'])

        final_list = []
        for t in init_list:
            DT = t.split('_')[0]
            y, m, d, h, mi, s = [int(x) for x in [DT[:4], DT[4:6], DT[6:8],
                                 DT[9:11], DT[11:13], DT[13:15]]]
            d = datetime(y, m, d, h, mi, s)
            p = (d - d0).seconds
            if (p > -10.) * (p < exptime+10.):
                final_list.append(t)
        for fn in final_list:
            fobj = T.extractfile(T.getmember(fn))
            f = fits.open(fobj)
            M.append([])
            if f[1].header['GUIDLOOP'] == 'ACTIVE':
                for att in attr:
                    M[-1].append(f[1].header[att])
    gmag = float(M[0][-1].split(',')[3])
    log.info("Guider header find g' mag: %0.2f" % gmag)
    if gmag < 0.:
        try:
            T1 = panstarrs_query(M[0]*15., M[1], 3. / 3600.)
        except:
            log.info('Could not get panstarrs match.')
            return 1.0
        if len(T1) == 0:
            log.info('No matches to panstarrs found.')
            return 1.0
        elif len(T1) > 1:
            log.info('Multiple matches within 3", using closest')
        gmag = T1['gMeanPSFMag'][0]
        if (gmag < 0.) + (gmag > 25.):
            log.info('Unreasonable g-mag from panstarrs: %0.2f' % gmag)
            return 1.0
    throughput = np.zeros((len(M),))
    for i, mi in enumerate(M):
        try:
            throughput[i] = 10**(-0.4 * (mi[4] - gmag))
        except:
            args.log.warning(mi)
    t = np.mean(throughput[throughput>0.0])
    log.info('Throughput for %s is %0.2f' % (path, t))
    return t


def check_if_standard(objname):
    standard_names = ['HD_19445', 'SA95-42', 'GD50', 'G191B2B', 'FEIGE_25',
                      'HILTNER_600', 'G193-74', 'PG0823+546', 'HD_84937',
                      'GD108', 'FEIGE_34', 'HD93521', 'GD140', 'HZ_21',
                      'FEIGE_66', 'FEIGE_67', 'G60-54', 'HZ_44', 'GRW+70_5824',
                      'BD+26+2606', 'BD+33_2642', 'G138-31', 'WOLF_1346',
                      'BD_+17_4708', 'FEIGE_110', 'GD248', 'HZ_4']
    for standard in standard_names:
        if standard.lower() in objname.lower():
            return True
    return False


def fit_response_cont(wv, sky, skip=5, fil_len=95, func=np.array):
    skym_s = 1. * sky
    sky_sm = savgol_filter(skym_s, fil_len, 1)
    allind = np.arange(len(wv), dtype=int)
    y = np.abs(np.diff(np.hstack([sky[0], sky])))
    sel = np.where(y > 1.5 * np.median(sky))[0]
    for j in np.arange(1, skip+1):
            sel = np.union1d(sel, sel + 1)
            sel = np.union1d(sel, sel - 1)
    sel = np.sort(np.unique(sel))
    sel = sel[skip:-skip]
    good = np.setdiff1d(allind, sel)
    skym_s = 1.*sky
    skym_s[sel] = np.interp(wv[sel], wv[good], sky_sm[good])
    sky_sm = savgol_filter(skym_s, fil_len, 1)
    for i in np.arange(5):
        mad = np.median(np.abs(sky - sky_sm))
        outlier = func(sky - sky_sm) < -1.5 * mad
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


def get_response(objname, commonwave, spec, specname):
    standard_names = ['HD_19445', 'SA95-42', 'GD50', 'G191B2B', 'FEIGE_25',
                      'HILTNER_600', 'G193-74', 'PG0823+546', 'HD_84937',
                      'GD108', 'FEIGE_34', 'HD93521', 'GD140', 'HZ_21',
                      'FEIGE_66', 'FEIGE_67', 'G60-54', 'HZ_44', 'GRW+70_5824',
                      'BD+26+2606', 'BD+33_2642', 'G138-31', 'WOLF_1346',
                      'BD_+17_4708', 'FEIGE_110', 'GD248', 'HZ_4']
    for standard in standard_names:
        if standard.lower() in objname.lower():
            filename = op.join('/work/03946/hetdex/maverick/virus_config/'
                               'standards',
                               'm' + standard.lower() + '.dat.txt')
            wave, standardmag = np.loadtxt(filename, usecols=(0, 1),
                                           unpack=True)
            fnu = 10**(0.4 * (-48.6 - standardmag))
            standard_flam = fnu * 2.99792e18 / wave**2
            standard_wave = wave
            flam = np.interp(commonwave, standard_wave, standard_flam)
            cont = fit_response_cont(commonwave, spec / flam, fil_len=11)
            return 1. / cont
    return None


def big_reduction(obj, bf, instrument, sci_obs, calinfo, amps, commonwave,
                  ifuslot, specname, standard=False, response=None):
    log.info('Extracting %s from %s' % (obj[0], bf))
    scifiles = sci_path % (instrument, instrument, sci_obs, '*',
                           instrument, ifuslot)
    images, rect, spec, cos, fl, Fi, E = extract_sci(scifiles, amps, calinfo[2],
                                              calinfo[1], calinfo[0], calinfo[3],
                                              calinfo[4], calinfo[5])
    cnt = 1
    wave_0 = np.mean(commonwave)
    darfile = op.join(DIRNAME, 'lrs2_config/dar_%s.dat' % specinit)
    T = Table.read(darfile, format='ascii.fixed_width_two_line')
    xoff = (np.interp(commonwave, T['wave'], T['x_0']) -
            np.interp(wave_0, T['wave'], T['x_0']))
    yoff = (np.interp(commonwave, T['wave'], T['y_0']) -
            np.interp(wave_0, T['wave'], T['y_0']))
    for im, r, s, c, fli, Fii, e in zip(images, rect, spec, cos, fl, Fi, E):
        fn = (sci_path % (instrument, instrument, sci_obs,
                          '%02d' % cnt, instrument, ifuslot))
        fn = glob.glob(fn)
        mini = get_objects(fn, ['OBJECT', 'EXPTIME'], full=True)
        log.info('Subtracting sky %s, exp%02d' % (obj[0], cnt))
        r[calinfo[-3][:, 1] == 1.] = 0.
        e[calinfo[-3][:, 1] == 1.] = 0.

        r /= mini[0][1]
        r /= mini[0][2]
        r /= mini[0][3]
        e /= mini[0][1]
        e /= mini[0][2]
        e /= mini[0][3]
        r, e = correct_ftf(r, e)
        sky = sky_subtraction(r, e)
        sky[calinfo[-3][:, 1] == 1.] = 0.
        skysub = r - sky
        if response is not None:
            r *= response
            e *= response
            sky *= response
            skysub *= response
        X = np.array([T['wave'], T['x_0'], T['y_0']])
        for S, name in zip([sky, skysub], ['sky', 'skysub']):
            outname = ('%s_%s_%s_%s_%s_cube.fits' % (args.date, sci_obs,
                       'exp%02d' % cnt, specname, name))
            zcube, zimage, xgrid, ygrid = make_frame(calinfo[5][:, 0],
                                                     calinfo[5][:, 1], S,
                                                     commonwave,
                                                     T['wave'],
                                                     xoff, yoff,
                                                     wstart=wave_0-50.,
                                                     wend=wave_0+50.)
            write_cube(commonwave, xgrid, ygrid, zcube, outname)
        loc = find_source(zimage, xgrid, ygrid)
        if loc is not None:
            log.info('Source found at %0.2f, %0.2f' % (loc[0], loc[1]))
            skyspec, errorskyspec = extract_source(sky, loc[0], loc[1], xoff, yoff,
                                     commonwave, calinfo[5][:, 0],
                                     calinfo[5][:, 1], e)
            skysubspec, errorskysubspec = extract_source(skysub, loc[0], loc[1], xoff, yoff,
                                        commonwave, calinfo[5][:, 0],
                                        calinfo[5][:, 1], e)
        else:
            skyspec = commonwave * 0.
            skysubspec = commonwave * 0.
            errorskyspec = commonwave * 0.
            errorskysubspec = commonwave * 0.
        if response is not None:
            f5 = fits.ImageHDU(np.vstack([commonwave, skysubspec, skyspec,
                                          errorskysubspec, errorskyspec, response]))
        else:
            f5 = fits.ImageHDU(np.vstack([commonwave, skysubspec, skyspec,
                                          errorskysubspec, errorskyspec]))
        outname = ('%s_%s_%s_%s_%s.fits' % ('multi', args.date, sci_obs,
                                            'exp%02d' % cnt, specname))
        cnt += 1
        f1 = create_header_objection(commonwave, r, func=fits.PrimaryHDU)
        f2 = create_header_objection(commonwave, sky)
        f3 = create_header_objection(commonwave, skysub)
        f4 = create_image_header(commonwave, xgrid, ygrid, zimage)
        f6 = create_header_objection(commonwave, e)

        fits.HDUList([f1, f2, f3, f6, f4, fits.ImageHDU(calinfo[5]), f5,
                      fits.ImageHDU(X), fits.ImageHDU(calinfo[3]),
                      fits.ImageHDU(im), fits.ImageHDU(fli), fits.ImageHDU(Fii),
                      fits.ImageHDU(c), fits.ImageHDU(s)]).writeto(outname, overwrite=True)
        if standard:
            return get_response(obj[0], commonwave, skysubspec, specname)


# LRS2-R
fiberpos, fiberspec = ([], [])
log.info('Beginning the long haul.')
allflatspec, allspec, allra, alldec, allx, ally, allsub = ([], [], [], [], [],
                                                           [], [])

DIRNAME = get_script_path()

for info in [redinfo[0], redinfo[1]]:
    specinit, specname, multi, lims, amps, slims, arc_names = info
    arc_lines = Table.read(op.join(DIRNAME, 'lrs2_config/lines_%s.dat' %
                                   specname), format='ascii')
    commonwave = np.linspace(lims[0], lims[1], 2064)
    specid, ifuslot, ifuid = multi.split('_')
    package = []
    flt_check_path = op.join(baseraw, args.date,  'lrs2', 'lrs20000*', 'exp01',
                             'lrs2', '2*_056LL_flt.fits')
    flt_files = sorted(glob.glob(flt_check_path))
    for fn in flt_files:
        o = fits.open(fn)[0].header['OBJECT']
        if specname in ['uv', 'orange']:
            if 'ldls' in o.lower():
                fltobs = op.basename(op.dirname(op.dirname(op.dirname(fn))))
        o = fits.open(fn)[0].header['OBJECT']
        if specname in ['red', 'farred']:
            if 'qth' in o.lower():
                fltobs = op.basename(op.dirname(op.dirname(op.dirname(fn))))
    twiflt_path = op.join(baseraw, twi_date,  '%s', '*', 'exp*',
                          '%s', '2*_%sLL_twi.fits')
    twibase = twiflt_path % (instrument, instrument, ifuslot)
    for amp in amps:
        amppos = get_ifucenfile(specname, amp)
        ##############
        # MASTERBIAS #
        ##############
        log.info('Getting Masterbias for ifuslot, %s, and amp, %s' %
                 (ifuslot, amp))
        zro_path = bias_path % (instrument, instrument, '00000*', instrument,
                                ifuslot)
        masterbias = get_masterbias(zro_path, amp)

        #####################
        # MASTERTWI [TRACE] #
        #####################
        log.info('Getting MasterFlat for ifuslot, %s, and amp, %s' %
                 (ifuslot, amp))
        masterflt = get_mastertwi(twibase, amp, masterbias)
        #twipath = twi_path % (instrument, instrument, '00000*', instrument,
        #                      ifuslot)
        #mastertwi = get_mastertwi(twibase, amp, masterbias)

        log.info('Getting Trace for ifuslot, %s, and amp, %s' %
                 (ifuslot, amp))
        trace, dead = get_trace(masterflt, specid, ifuslot, ifuid, amp,
                                twi_date)
        fits.PrimaryHDU(trace).writeto('test_trace.fits', overwrite=True)

        ##########################
        # MASTERARC [WAVELENGTH] #
        ##########################
        log.info('Getting MasterArc for ifuslot, %s, and amp, %s' %
                 (ifuslot, amp))
        lamp_path = cmp_path % (instrument, instrument, '00000*', instrument,
                                ifuslot)
        masterarc = get_masterarc(lamp_path, amp, arc_names, masterbias,
                                  specname)
        fits.PrimaryHDU(masterarc).writeto('wtf_%s_%s.fits' % (ifuslot, amp), overwrite=True)
        log.info('Getting Wavelength for ifuslot, %s, and amp, %s' %
                 (ifuslot, amp))
        wave = get_wavelength_from_arc(masterarc, trace, arc_lines, specname, amp)
        fits.PrimaryHDU(wave).writeto('test_wave.fits', overwrite=True)

        #################################
        # TWILIGHT FLAT [FIBER PROFILE] #
        #################################
        log.info('Getting bigW for ifuslot, %s, and amp, %s' %
                 (ifuslot, amp))
        bigW = get_bigW(amp, wave, trace, masterbias)
        package.append([wave, trace, bigW, masterbias, amppos, dead])
    # Normalize the two amps and correct the flat
    calinfo = [np.vstack([package[0][i], package[1][i]])
               for i in np.arange(len(package[0]))]
    calinfo[1][package[0][1].shape[0]:, :] += package[0][2].shape[0]
    log.info('Getting flat for ifuslot, %s, side, %s' % (ifuslot, specname))
    twiflat = get_twiflat_field(twibase, amps, calinfo[0], calinfo[1],
                                calinfo[2], commonwave, calinfo[3], specname)
    calinfo.insert(2, twiflat)
    flatspec = get_spectra(calinfo[2], calinfo[1])
    calinfo.append(flatspec)
    bigF = get_bigF(calinfo[1], calinfo[2])
    calinfo.append(bigF)
    #####################
    # SCIENCE REDUCTION #
    #####################
    basefiles = sorted(glob.glob(sci_path % (instrument, instrument, '0000*',
                                             '01', instrument, ifuslot)))
    all_sci_obs = [op.basename(op.dirname(op.dirname(op.dirname(fn))))[-7:]
                   for fn in basefiles]
    objects = get_objects(basefiles, ['OBJECT', 'EXPTIME'])

    response = None
    for sci_obs, obj, bf in zip(all_sci_obs, objects, basefiles):
        if check_if_standard(obj[0]) and (ifuslot in obj[0]):
            log.info('Getting Response Function from %s' % obj[0])
            response = big_reduction(obj, bf, instrument, sci_obs, calinfo,
                                     amps, commonwave, ifuslot, specname,
                                     standard=True)
    f = []
    names = ['wavelength', 'trace', 'flat', 'bigW', 'masterbias', 'xypos',
             'dead', 'flatspec', 'bigF']
    for i, cal in enumerate(calinfo):
        if i == 0:
            func = fits.PrimaryHDU
        else:
            func = fits.ImageHDU
        f.append(func(cal))
    if response is not None:
        f.append(fits.ImageHDU(np.array([commonwave, response], dtype=float)))
        names.append('response')
    for fi, n in zip(f, names):
        fi.header['EXTNAME'] = n
    fits.HDUList(f).writeto('cal_%s_%s.fits' % (args.date, specname),
                            overwrite=True)
    for sci_obs, obj, bf in zip(all_sci_obs, objects, basefiles):
        if args.object is None:
            big_reduction(obj, bf, instrument, sci_obs, calinfo, amps, commonwave,
                          ifuslot, specname, response=response)
        else:
            if args.object.lower() in obj[0].lower():
                big_reduction(obj, bf, instrument, sci_obs, calinfo, amps, commonwave,
                          ifuslot, specname, response=response)
            

#        header = fits.open(glob.glob(sci_path % (instrument, instrument, sci_obs, '01',
#                                         instrument,
#                                         ifuslots[0]))[0])[0].header
#        PA = float(header['PARANGLE'])
#        RA = float(header['TRAJRA'])
#        DEC = float(header['TRAJDEC'])
#        log.info('Observation at %0.4f %0.4f, PA: %0.3f' % (RA, DEC, PA))
#        A = Astrometry(RA, DEC, PA, 0., 0., fplane_file=fplane_file)
#        for i in np.arange(nexp):
#            log.info('Getting RA, Dec for exposure, %i, ifuslot, %s, and amp,'
#                     ' %s' % (i+1, ifuslot, amp))
#            ra, dec = A.get_ifupos_ra_dec(ifuslot,
#                                          amppos[:, 0] + dither_pattern[i, 0],
#                                          amppos[:, 1] + dither_pattern[i, 1])
#            allra.append(ra)
#            alldec.append(dec)
#            allx.append(A.fplane.by_ifuslot(ifuslot).y + amppos[:, 0] +
#                        dither_pattern[i, 0])
#            ally.append(A.fplane.by_ifuslot(ifuslot).x + amppos[:, 1] +
#                        dither_pattern[i, 1])
