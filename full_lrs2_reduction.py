# -*- coding: utf-8 -*-
"""
Created on Thu Oct  4 13:30:06 2018

@author: gregz
"""
import numpy as np
import fnmatch
import os.path as op
import os
import glob
import sys
import warnings
import tarfile
import argparse as ap
import requests
import uuid

from datetime import datetime, timedelta
from astropy.io.votable import parse_single_table
from fiber_utils import bspline_x0

from astropy.io import fits
from astropy.table import Table
from scipy.signal import savgol_filter
from distutils.dir_util import mkpath
from scipy.ndimage.filters import percentile_filter
from scipy.interpolate import interp1d, interp2d, griddata
from input_utils import setup_logging
from astrometry import Astrometry
from astropy.stats import biweight_midvariance, biweight_location
from astropy.modeling.models import Moffat2D, Polynomial2D, Gaussian2D
from astropy.modeling.fitting import LevMarLSQFitter
from sklearn.decomposition import PCA
from astropy.convolution import Gaussian1DKernel, Gaussian2DKernel, convolve


standard_names = ['HD_19445', 'SA95-42', 'GD50', 'G191B2B',
                  'HILTNER_600', 'G193-74', 'PG0823+546', 'HD_84937',
                  'GD108', 'FEIGE_34', 'HD93521', 'GD140', 'HZ_21',
                  'FEIGE_66', 'FEIGE_67', 'G60-54', 'HZ_44', 'GRW+70_5824',
                  'BD+26+2606', 'BD+33_2642', 'G138-31', 'WOLF_1346',
                  'BD_+17_4708', 'FEIGE_110', 'GD248', 'HZ_4',
                  'BD+40_4032', 'HILTNER_102',
                  'BD_+26_2606', 'GD_248', 'FEIGE_56', 'FEIGE_92',
                  'HZ_15', 'FEIGE_98', 'BD+08_2015', 'BD+25_3941',
                  'FEIGE_15', 'FEIGE_25', 'SA_95-42', 'BD+28_4211',
                  'HR6203']

log = setup_logging('panacea_quicklook')

parser = ap.ArgumentParser(add_help=True)

parser.add_argument("-d", "--date",
                    help='''Date for reduction''',
                    type=str, default='20181108')

parser.add_argument("-s", "--sides",
                    help='''"uv,orange,red,farred"''',
                    type=str, default="uv,orange,red,farred")

parser.add_argument("-o", "--object",
                    help='''Object name, no input reduces all objects''',
                    type=str, default=None)

parser.add_argument("-uf", "--use_flat",
                    help='''Use FLT instead of Twi''',
                    action="count", default=0)

parser.add_argument("-cf", "--correct_ftf",
                    help='''Correct fiber to fiber''',
                    action="count", default=0)

parser.add_argument("-md", "--model_dar",
                    help='''model DAR''',
                    action="count", default=0)

parser.add_argument("-cw", "--central_wave",
                    help='''Central Wavelength for collapsed Frame''',
                    type=float, default=None)

parser.add_argument("-wb", "--wavelength_bin",
                    help='''Wavelength Bin to collapse over (+/- bin size)''',
                    type=float, default=10.)

parser.add_argument("-sx", "--source_x",
                    help='''Source's x position at the central_wave''',
                    type=float, default=None)

parser.add_argument("-sy", "--source_y",
                    help='''Source's y position at the central_wave''',
                    type=float, default=None)

parser.add_argument("-ssd", "--standard_star_date",
                    help='''Standard Star Date for response function,
                    example: 20181101''',
                    type=str, default=None)

parser.add_argument("-sso", "--standard_star_obsid",
                    help='''Standard Star ObsID for response function,
                    example: 0000012''',
                    type=str, default=None)

parser.add_argument("-ad", "--arc_date",
                    help='''Arc Date for reduction''',
                    type=str, default=None)

parser.add_argument("-td", "--twi_date",
                    help='''Twilight Date for reduction''',
                    type=str, default=None)

parser.add_argument("-re", "--reduce_eng",
                    help='''Reduce Engineer Data''',
                    action="count", default=0)

parser.add_argument("-rf", "--reduce_flt",
                    help='''Reduce Flat Data''',
                    action="count", default=0)

parser.add_argument("-rd", "--reduce_drk",
                    help='''Reduce Dark Data''',
                    action="count", default=0)

args = parser.parse_args(args=None)

if args.standard_star_obsid is not None:
    args.standard_star_obsid = '%07d' % int(args.standard_star_obsid)
    if args.standard_star_date is None:
        log.error('Please include --standard_star_date DATE with call.')

for i in ['source_x', 'source_y']:
    for j in ['source_x', 'source_y', 'central_wave']:
        if i == j:
            continue
        if getattr(args, i) is not None:
            if getattr(args, j) is None:
                log.error('%s was set but not %s.' % (i, j))
                sys.exit(1)

args.sides = [x.replace(' ', '') for x in args.sides.split(',')]


blueinfo = [['BL', 'uv', '503_056_7001', [3640., 4645.], ['LL', 'LU'],
             [4350., 4375.], ['hg_b', 'cd-a_b', 'fear_r', 'cd_b', 'hg', 'cd', 'fear']],
            ['BR', 'orange', '503_056_7001',
             [4635., 6950.], ['RU', 'RL'], [6270., 6470.],
             ['hg_b', 'cd-a_b', 'fear_r', 'cd_b', 'hg', 'cd', 'fear']]]
redinfo = [['RL', 'red', '502_066_7002', [6450., 8400.], ['LL', 'LU'],
            [7225., 7425.], ['hg_r', 'cd-a_b', 'fear_r', 'cd_b', 'hg', 'cd', 'fear']],
           ['RR', 'farred', '502_066_7002',
            [8275., 10500.], ['RU', 'RL'], [9280., 9530.],
            ['hg_r', 'cd-a_b', 'fear_r', 'cd_b', 'hg', 'cd', 'fear']]]

listinfo = []
for side in args.sides:
    if side.lower() == 'uv':
        listinfo.append(blueinfo[0])
    if side.lower() == 'orange':
        listinfo.append(blueinfo[1])
    if side.lower() == 'red':
        listinfo.append(redinfo[0])
    if side.lower() == 'farred':
        listinfo.append(redinfo[1])


fplane_file = '/work/03730/gregz/maverick/fplane.txt'
twi_date = args.date
sci_date = args.date

# FOR LRS2
instrument = 'lrs2'

dither_pattern = np.zeros((50, 2))


baseraw = '/work/03946/hetdex/maverick'


sci_tar  = op.join(baseraw, sci_date,  '%s', '%s000*.tar')
sci_path = op.join(baseraw, sci_date,  '%s', '%s%s', 'exp%s',
                   '%s', '2*_%sLL*sci.fits')
cmp_path = op.join(baseraw, sci_date,  '%s', '%s%s', 'exp*',
                   '%s', '2*_%sLL_cmp.fits')
twi_path = op.join(baseraw, sci_date,  '%s', '%s%s', 'exp*',
                   '%s', '2*_%sLL_twi.fits')
bias_path = op.join(baseraw, sci_date, '%s', '%s%s', 'exp*',
                    '%s', '2*_%sLL_zro.fits')

if args.reduce_eng:
    sci_path = sci_path.replace('sci', 'eng')

if args.reduce_flt:
    sci_path = sci_path.replace('sci', 'flt')
    
if args.reduce_drk:
    sci_path = sci_path.replace('sci', 'drk')

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


def base_reduction(filename, tarname=None, get_header=False):
    if tarname is None:
        a = fits.open(filename)
    else:
        try:
            t = tarfile.open(tarname, 'r')
            a = fits.open(t.extractfile('/'.join(filename.split('/')[-4:])))
        except:
            log.warning('Could not open %s' % filename)
            return np.zeros((1032, 2064)), np.zeros((1032, 2064))
    image = np.array(a[0].data, dtype=float)
    # overscan sub
    overscan_length = int(32 * (image.shape[1] / 1064))
    O = biweight_location(image[:, -int(overscan_length-2):])
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
    header = a[0].header
    a = orient_image(image, amp, ampname) * gain
    E = np.sqrt(rdnoise**2 + np.where(a > 0., a, 0.))
    if tarname is not None:
        t.close()
    if get_header:
        return a, E, header
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
    xy = np.nanmedian(spec, axis=1)
    bottom = np.abs(np.nanpercentile(xy, 15))
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
    F0 = array_trace[:, int(m/2)]
    for j in np.arange(image.shape[1]):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            p0 = np.polyfit(array_trace[:, j], F0, 7)
        bigF[:, j] = np.polyval(p0, YY[:, j])
    return bigF

def get_tarname_from_filename(filename):
    tarname = op.dirname(op.dirname(op.dirname(filename))) + '.tar'
    return tarname

def get_twiflat_field(files, amps, array_wave, array_trace, bigW,
                      common_wave, masterbias, specname):
    files1 = [file.replace('LL', amps[0]) for file in files]
    files2 = [file.replace('LL', amps[1]) for file in files]
    tarnames = [get_tarname_from_filename(file) for file in files]
    array_list = []
    for filename1, filename2, tarname in zip(files1, files2, tarnames):
        log.info('Prepping flat %s' % filename1)
        array_flt1, e1 = base_reduction(filename1, tarname=tarname)
        array_flt2, e2 = base_reduction(filename2, tarname=tarname)
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
        try:
            spectrum[fiber] = array_flt[indl, x] / 2. + array_flt[indh, x] / 2.
        except:
            spectrum[fiber] = 0.
    
    log.info('Getting powerlaw for side %s' % specname)
    plaw, norm = get_powerlaw(array_flt, array_trace, spectrum)
    array_flt[:] -= plaw
    array_flt[:] = np.where(array_flt < 0., 0., array_flt)
    #fits.PrimaryHDU(plaw).writeto('test_plaw_%s.fits' % specname, overwrite=True)
    #fits.PrimaryHDU(array_flt).writeto('test_spec_%s.fits' % specname, overwrite=True)
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
                             per=50)
    B, c = bspline_x0(nw1, nknots=int(spectrum.shape[1]))
    sol = np.linalg.lstsq(c, ns1)[0]
    ns1 = np.dot(c, sol)
    I = interp1d(nw1, ns1, kind='quadratic', fill_value='extrapolate')
    modelimage = I(bigW)
    flat = array_flt / modelimage
    flat[~np.isfinite(flat)] = 0.0
    flat[flat < 0.0] = 0.0
    #fits.HDUList([fits.PrimaryHDU(array_flt), fits.ImageHDU(modelimage),
    #                 fits.ImageHDU(flat)]).writeto('flat_example_%s.fits' % specname, overwrite=True)
    return flat


def get_spectra(array_flt, array_trace):
    spectrum = array_trace * 0.
    x = np.arange(array_flt.shape[1])
    for fiber in np.arange(array_trace.shape[0]):
        indl = np.floor(array_trace[fiber]).astype(int)
        indh = np.ceil(array_trace[fiber]).astype(int)
        try:
            spectrum[fiber] = array_flt[indl, x] / 2. + array_flt[indh, x] / 2.
        except:
            log.warning('Index for getting spectrum out of bounds.')
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


def weighted_extraction(image, error, flat, trace, cthresh=8.):
    E = safe_division(error, flat)
    E[E < 1e-8] = 1e9
    Y = safe_division(image, flat)
    nY = Y * 1.
    C = np.array(Y * 0., dtype=bool)
    for i in np.arange(1):
        cosmics = find_cosmics(nY, E, trace, thresh=cthresh, ran=1)
        C = C + cosmics

    x = np.arange(trace.shape[1])
    spectrum = 0. * trace
    error_spec = 0. * trace
    Fimage = image * 0.
    for fiber in np.arange(trace.shape[0]):
        T = np.zeros((4, trace.shape[1], 4))
        indl = np.floor(trace[fiber]).astype(int)-1
        flag = False
        for ss, k in enumerate(np.arange(0, 4)):
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
    mid = int(Trace.shape[1] / 2)
    shifts = np.nanmedian((FlatTrace - Trace)[:, mid-200:mid+200], axis=1)
    shifts = np.polyval(np.polyfit(np.nanmedian(FlatTrace, axis=1), shifts, 2),
                        Yx)
    return shifts


def modify_spectrum(spectrum, error, w, xloc, yloc):
    dw = np.median(np.diff(w, axis=1), axis=0)
    dw = np.hstack([dw[0], dw])
    for i in np.arange(spectrum.shape[0]):
        sel = spectrum[i] == 0.
        I = interp1d(w[i][~sel], spectrum[i][~sel], kind='quadratic',
                     fill_value='extrapolate')
        spectrum[i] = I(w[i]) / dw
        error[i] /= dw
    bad = error == 0.
    for i in np.arange(1, 3):
        bad[:, :-i] += bad[:, i:]
        bad[:, i:] += bad[:, :-i]
    error[bad] = 0.
    return spectrum, error


def extract_sci(sci_path, amps, flat, array_trace, array_wave, bigW,
                masterbias, pos):
    files1 = get_filenames_from_tarfolder(get_tarname_from_filename(sci_path),
                                          sci_path.replace('LL', amps[0]))
    files2 = get_filenames_from_tarfolder(get_tarname_from_filename(sci_path),
                                          sci_path.replace('LL', amps[1]))
    for filen in files1:
        log.info('SECOND --- filename: %s' % filen)
    tarnames = [get_tarname_from_filename(file) for file in files1]
    xloc, yloc = (pos[:, 0], pos[:, 1])
    array_list, hdr_list = ([], [])
    for filename1, filename2, tarname in zip(files1, files2, tarnames):
        log.info('Prepping sci %s' % filename1)
        array_flt1, e1, header = base_reduction(filename1, tarname=tarname,
                                                get_header=True)
        array_flt2, e2 = base_reduction(filename2, tarname=tarname)
        array_flt = np.vstack([array_flt1, array_flt2])
        array_flt[:] -= masterbias
        array_list.append(array_flt)
        hdr_list.append(header)
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
    for filename1, filename2, tarname in zip(files1, files2, tarnames):
        log.info('Fiber extraction sci %s' % filename1)
        array_flt1, e1 = base_reduction(filename1, tarname=tarname)
        array_flt2, e2 = base_reduction(filename2, tarname=tarname)
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
        spectrum, error = modify_spectrum(spectrum, error, array_wave, xloc,
                                          yloc)
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
            sel = error[fiber] == 0.
            if sel.sum() > 0.:
                nsel = coV[:, sel].sum(axis=1) > 0.
                error_interp[nsel] = 0.
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
            np.array(error_list), hdr_list)


def get_masterbias(files, amp):
    files = [file.replace('LL', amp) for file in files]
    tarnames = [get_tarname_from_filename(file) for file in files]
    
    biassum = np.zeros((len(files), 1032, 2064))
    for j, filename in enumerate(files):
        tarname = tarnames[j]
        a, error = base_reduction(filename, tarname=tarname)
        biassum[j] = a
    return biweight_location(biassum, axis=0)


def get_masterarc(files, amp, arc_names, masterbias, specname, trace):
    files = [file.replace('LL', amp) for file in files]
    tarnames = [get_tarname_from_filename(file) for file in files]
    arcsum = np.zeros((1032, 2064))
    cnt = np.zeros((1032, 2064))
    for filename, tarname in zip(files, tarnames):
        t = tarfile.open(tarname, 'r')
        f = fits.open(t.extractfile('/'.join(filename.split('/')[-4:])))
        if f[0].header['OBJECT'].lower() in arc_names:
            a, e = base_reduction(filename, tarname=tarname)
            a[:] -= masterbias
            if np.median(a) < 3000.:
                #c = find_cosmics(a, e, trace, thresh=15., ran=0)
                arcsum += a
                cnt += 1.
    return arcsum / cnt


def get_mastertwi(files, amp, masterbias):
    files = [file.replace('LL', amp) for file in files]
    tarnames = [get_tarname_from_filename(file) for file in files]
    listtwi = []
    for filename, tarname in zip(files, tarnames):
        a, e = base_reduction(filename, tarname=tarname)
        a[:] -= masterbias
        if np.percentile(a, 75) > 100.:
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


def find_lines(spectrum, trace, nlines, thresh, fib, side=None):
    cont = percentile_filter(spectrum, 15, (1, 101))
    spectrum -= cont
    loc = []
    ph, pr = ([], [])
    lines = Table(nlines)
    for i, spec in enumerate(spectrum):
        px, ps, py = find_peaks(spec, thresh=thresh)
        sel = np.abs(px - 1032.) > 0.
        loc.append(px[sel])
        ph.append(ps[sel])
        pr.append(py[sel])
    
    if side == 'orange':
        names = ['Hg', 'Cd']
        v = []
        for name in names:
            selhg = lines['col4'] == name
            ma = np.argmax(arc_lines['col3'][selhg])
            sel = np.abs(loc[fib] - lines['col2'][selhg][ma]) < 50.
            v1 = np.max(pr[fib][sel])
            v2 = lines['col3'][selhg][ma]
            v.append([v1, v2])
        selhg = lines['col4'] == name
        if v[0][0] > v[1][0]:
            selhg = lines['col4'] == 'Hg'
            ma = np.argmax(arc_lines['col3'][selhg])
            mxv = lines['col3'][selhg][ma]
            lines['col3'][selhg] *= 1. / mxv
            selhg = (lines['col4'] == 'Cd') + (lines['col4'] == 'Ar')
            ma = np.argmax(arc_lines['col3'][selhg])
            mxv = lines['col3'][selhg][ma]
            nrt = v[1][0] / v[0][0]
            lines['col3'][selhg] *= nrt / mxv 
        else:
            selhg = (lines['col4'] == 'Cd') + (lines['col4'] == 'Ar')
            ma = np.argmax(arc_lines['col3'][selhg])
            mxv = lines['col3'][selhg][ma]
            lines['col3'][selhg] *= 1. / mxv
            selhg = lines['col4'] == 'Hg'
            ma = np.argmax(arc_lines['col3'][selhg])
            mxv = lines['col3'][selhg][ma]
            nrt = v[0][0] / v[1][0]
            lines['col3'][selhg] *= nrt / mxv
    found_lines = np.zeros((trace.shape[0], len(lines)))
    ls = np.argsort(lines['col3'])[::-1]
    if side == 'farred':
        distthresh = 15.
    else:
        distthresh = 50.
    sel = np.abs(loc[fib] - lines['col2'][ls[0]]) < distthresh
    ind = np.where(sel)[0][np.argmax(pr[fib][sel])]
    off = loc[fib][ind] - lines['col2'][ls[0]]

    found_lines[fib, ls[0]] = loc[fib][ind]
    y = lines['col2'] + off
    s = found_lines[fib] * 0.
    pp = s * 0.
    pph = s * 0.
    s[ls[0]] = 1.
    pp[ls[0]] = 0.
    for l in ls[1:]:
        guess = y[l]
        v = np.abs(guess - loc[fib])
        ER = lines['col3'][l] / lines['col3'][ls[0]]
        MR = pr[fib] / pr[fib][ind]
        EE = MR * np.sqrt(1./ph[fib]**2 + 1./ph[fib][ind])
        EE = np.max([EE, .1 * MR, 0.001 * np.ones(MR.shape)], axis=0)
        dist = v/2. + np.abs(ER - MR) / EE / 2.
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
            pph[l] = ph[fib][ind1]
    inds = np.where(found_lines[fib] > 0.)[0]
    delv = []
    for ind in inds:
        sel = np.where(found_lines[fib, ind] == found_lines[fib, inds])[0]
        if len(sel)>1:
            if np.any(pp[ind] > pp[inds[sel]]):
                delv.append(ind)
                found_lines[fib, ind] = 0.
    inds = np.delete(inds, delv)
    for ind in inds:
        print(lines['col1'][ind], lines['col2'][ind], found_lines[fib][ind],
              lines['col3'][ind], s[ind], pp[ind], pph[ind])    
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
            if np.min(m) < 2.0:
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
            if np.min(m) < 2.0:
                found_lines[j, i] = loc[j][np.argmin(m)]
    return found_lines


def get_wavelength_from_arc(image, trace, lines, side, amp, date, otherimage=None):
    spectrum = get_spectra(image, trace)
    fib = np.argmax(np.median(spectrum, axis=1))
    if side == 'uv' and date > 20161101:
        thresh = 5.  # 5
        spectrum2 = spectrum*0.
        fib = trace.shape[0] / 2
        for i in np.arange(trace.shape[0]):
            ll = int(np.max([0, i-4]))
            hl = int(np.min([trace.shape[0], i+5]))
            spectrum2[i] = np.median(spectrum[ll:hl], axis=0)
            G = Gaussian1DKernel(1.5)
            spectrum2[i] = convolve(spectrum2[i], G)
        spectrum = spectrum2 * 1.
        spectrum1 = get_spectra(otherimage, trace)
        spectrum2 = spectrum1*0.
        fib = trace.shape[0] / 2
        for i in np.arange(trace.shape[0]):
            ll = int(np.max([0, i-4]))
            hl = int(np.min([trace.shape[0], i+5]))
            spectrum2[i] = np.median(spectrum1[ll:hl], axis=0)
            G = Gaussian1DKernel(1.5)
            spectrum2[i] = convolve(spectrum2[i], G)
        spectrum1 = spectrum2 * 1.
        fib = int(trace.shape[0] / 2)
        found_lines = find_lines(spectrum, trace, lines, thresh, fib)
        found_lines1 = find_lines(spectrum1, trace, lines, thresh, fib)
        sel = (found_lines > 0.) * (found_lines1 > 0.)
        for i in np.arange(found_lines.shape[0]):
            found_lines[i] = (found_lines1[i] +
                              np.median(found_lines[i][sel[i]] -
                                        found_lines1[i][sel[i]]))
            found_lines[i][found_lines1[i] == 0.] = 0.
    else:
        thresh = 3.
        found_lines = find_lines(spectrum, trace, lines, thresh, fib, side)

    x = np.arange(trace.shape[1])
    found_lines[found_lines > 2060] = 0.
    found_lines1 = found_lines * 1.
    #fits.PrimaryHDU(found_lines).writeto('fl_%s_%s.fits' % (side, amp), overwrite=True)
    qft = np.zeros((len(lines),))
    for i, line in enumerate(lines):
        if np.sum(found_lines[:, i]) < (0.5 * trace.shape[0]):
            found_lines[:, i] = 0.0
            continue
        ind = np.array(found_lines[:, i], dtype=int)
        xt = trace[np.arange(trace.shape[0]), ind]
        yt = robust_polyfit(xt, found_lines[:, i])
        sel = found_lines1[:, i] > 0.
        qft[i] = np.std(found_lines1[sel, i] - yt[sel])
        if qft[i] > 0.25:
            found_lines[:, i] = 0.0
            continue
        found_lines[:, i] = yt
    wave = trace * 0.0
    res = np.zeros((trace.shape[0],))
    Res = np.zeros(found_lines.shape)
    for j in np.arange(trace.shape[0]):
        sel = found_lines[j, :] > 0.0
        if sel.sum() < 3:
            continue
        wave[j] = np.polyval(np.polyfit(found_lines[j, sel],
                             lines['col1'][sel], 3), x)
        res[j] = np.std(np.interp(found_lines[j, sel], x, wave[j]) -
                        lines['col1'][sel])
        Res[j, sel] = (np.interp(found_lines1[j, sel], x, wave[j]) -
                                lines['col1'][sel])
    #fits.PrimaryHDU(Res).writeto('Rl_%s_%s.fits' % (side, amp), overwrite=True)

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
        tarname = get_tarname_from_filename(fn)
        t = tarfile.open(tarname, 'r')
        F = fits.open(t.extractfile('/'.join(fn.split('/')[-4:])))
        s.append([])
        for att in attrs:
            s[-1].append(F[0].header[att])
        if full:
            area = get_mirror_illumination_guider(fn, s[-1][1])
            try:
                throughput = get_throughput(fn, s[-1][1])
            except:
                log.warning('Problem with throughput for %s' % fn)
                throughput = 1.0
            s[-1].append(area)
            s[-1].append(throughput)
        t.close()
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
            y1 = savgol_filter(y[~o], 31, 1)
        I = interp1d(x[~o], y1, kind='quadratic', fill_value='extrapolate')
        y1 = I(x)
        for i in np.arange(3):
            o += outlier(y, y1, ~o)
            y1 = savgol_filter(y[~o], 51, 1)
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


def sky_subtraction(rect, error, xloc, yloc):
    y = np.median(rect, axis=1)
    selg = y != 0.
    v = np.percentile(y[selg], 5)
    init_sel = selg * (y < v)
    init = np.percentile(rect[init_sel], 50, axis=0)
    df = np.diff(init)
    df = np.hstack([df[0], df])
    cont = np.abs(df) < np.percentile(np.abs(df), 25)
    G = Gaussian1DKernel(15)
    tempy = init * 1.
    tempy[~cont] = np.nan
    smooth_back = convolve(tempy, G, nan_treatment='interpolate',
                           preserve_nan=False)
    peak_loc, sn, v = find_peaks(init - smooth_back, thresh = 3)
    locs = np.round(peak_loc).astype(int)
    locs = np.sort(np.hstack([locs-2, locs-1, locs, locs+1, locs+2]))
    locs = locs[np.where(locs>=0)[0]]
    locs = locs[np.where(locs<2064)[0]]
    facs = np.arange(0.7, 1.3, 0.01)
    cnt = facs * 0.
    sol = np.ones((280,))
    for k in np.arange(280):
        good = np.zeros(rect[k].shape, dtype=bool)
        good[locs] = True
        good[error[k] == 0.] = False
        if good.sum() > 50.:
            for j, i in enumerate(facs):
                cnt[j] = (((rect[k] - i * init) * init)[good]).sum()
            ind = np.argsort(cnt)
            sol[k] = np.interp(0., cnt[ind], facs[ind])
    sol[np.abs(sol - 1.3) < 0.01] = 1.
    n1 = np.median(sol[:140])
    n2 = np.median(sol[140:])
    nsol = sol * 1.
    nsol[:140] = sol[:140] / n1
    nsol[140:] = sol[140:] / n2
    fitter = LevMarLSQFitter()
    P = Polynomial2D(2)
    good = nsol != 0.
    fit = fitter(P, xloc[good], yloc[good], nsol[good])
    off = np.abs(sol - fit(xloc, yloc))
    mad = np.median(off)
    good = (nsol != 0.) * (off <= 2. * mad)
    fit = fitter(P, xloc[good], yloc[good], nsol[good])
    model = fit(xloc, yloc)
    model[:140] *= n1
    model[140:] *= n2
    sky = init * model[:, np.newaxis]
    sky[~selg] = 0.
    res = 0. * sky
    skysub = rect - sky
    for j in np.arange(rect.shape[1]):
        E = error[:, j] * 1.
        E[E == 0.] = 1e9
        W = 1. / E**2
        W = np.sqrt(np.diag(W))
        A = model[:, np.newaxis]
        B = skysub[:, j]
        Aw = np.dot(W, A)
        Bw = np.dot(B, W)
        sol = np.linalg.lstsq(Aw, Bw)[0][0]
        sg = np.sign(sol)
        V = [np.abs(sol), 0.5 * np.median(E)]
        mult = V[np.argmin(V)] * sg
        res[:, j] = mult
    return sky + res


def make_frame(xloc, yloc, data, error, wave, dw, Dx, Dy, wstart=5700.,
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
        sel = np.isfinite(data[:, k]) * (error[:, k] != 0.)
        D = np.sqrt((xloc[:, np.newaxis, np.newaxis] - Dx[k] - xgrid)**2 +
                    (yloc[:, np.newaxis, np.newaxis] - Dy[k] - ygrid)**2)
        W = np.exp(-0.5 / (seeing/2.35)**2 * D**2)
        zgrid[k, :, :] = ((data[sel, k][:, np.newaxis, np.newaxis] *
                           W[sel]).sum(axis=0) / W[sel].sum(axis=0) *
                          (scale**2 / area))
    wi = np.searchsorted(wave, wstart, side='left')
    we = np.searchsorted(wave, wend, side='right')

    zimage = np.median(zgrid[wi:we+1], axis=(0,))
    return zgrid, zimage, xgrid, ygrid


def write_cube(wave, xgrid, ygrid, zgrid, outname, he):
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
    for key in he.keys():
        if key in hdu.header:
            continue
        if ('CCDSEC' in key) or ('DATASEC' in key):
            continue
        if ('BSCALE' in key) or ('BZERO' in key):
            continue
        try:
            hdu.header[key] = he[key]
        except:
            continue
    hdu.writeto(outname, overwrite=True)


def build_weight_matrix(x, y, sig=1.5):
    d = np.sqrt((x - x[:, np.newaxis])**2 + (y - y[:, np.newaxis])**2)
    G = np.exp(-0.5 * (d / sig)**2)
    G = G / G.sum(axis=0)[:, np.newaxis]
    return G.swapaxes(0, 1)


def mask_skylines_cosmics(wave, rect_spec, name, error):
    mask1 = rect_spec * 0.
    if op.exists(op.join(DIRNAME, 'lrs2_config', '%s_skylines.dat' % name)):
        T = Table.read(op.join(DIRNAME, 'lrs2_config', '%s_skylines.dat' % name),
                       format='ascii.fixed_width_two_line')
        for w in T['wavelength']:
            mask1[:, np.abs(wave - w) < 6.] = -1.
    mask2 = rect_spec * 0.
    mask2[error == 0.] = -1.
    mask2[1:, :] += mask2[:-1, :]
    mask2[:-1, :] += mask2[1:, :]
    if name == 'uv':
        mask1[79:83, 979:987] = -1.
    mask = (mask1 + mask2) < 0
    return mask


def get_all_cosmics(x, y, ispec, error):
    D = np.sqrt((x - x[:, np.newaxis])**2 + (y - y[:, np.newaxis])**2)
    for i in np.arange(D.shape[0]):
        D[i, :] = np.array(D[i, :] < 1.5, dtype=float)
    ispec[error==0.] = 0.
    T = ispec * 1.
    for i in np.arange(ispec.shape[1]):
        T[:, i] = np.dot(ispec[:, i], D)
    YY = ispec / T
    YY[np.isnan(YY)] = 0.
    return YY > 0.2

def convolve_spatially(x, y, spec, wave, name, error, ispec, sig_spatial=0.75,
                       sig_wave=1.5):
    W = build_weight_matrix(x, y, sig=sig_spatial)
    D = np.sqrt((x - x[:, np.newaxis])**2 + (y - y[:, np.newaxis])**2)
    for i in np.arange(D.shape[0]):
        D[i, :] = np.array(D[i, :] < 1.5, dtype=float)
    mask = mask_skylines_cosmics(wave, spec, name, error)
    Z = spec * 1.
    E = error**2
    Z[mask] = np.nan
    E[mask] = np.nan
    ispec[mask] = 0.
    T = ispec * 1.
    for i in np.arange(ispec.shape[1]):
        T[:, i] = np.dot(ispec[:, i], D)
    YY = ispec / T
    YY[np.isnan(YY)] = 0.
    Z[YY > 0.2] = np.nan
    E[YY > 0.2] = np.nan
    G = Gaussian1DKernel(sig_wave)
    for i in np.arange(spec.shape[0]):
        Z[i, :] = convolve(Z[i, :], G, nan_treatment='fill', fill_value=0.0)
        E[i, :] = convolve(E[i, :], G, nan_treatment='fill', fill_value=0.0)
    Z_copy = Z * 1.
    E_copy = np.sqrt(E)
    T = spec * 0.
    for i in np.arange(spec.shape[1]):
        Z[:, i] = np.dot(Z[:, i], W)
        E[:, i] = np.dot(E[:, i], W)
    E[:] = np.sqrt(E)
    Y = Z * 0.
    sel = E > 0.
    Y[sel] = Z[sel] / E[sel]
    Y[~np.isfinite(Y)] = 0.
    ind = np.unravel_index(np.nanargmax(Y[:, 50:-50],
                                        axis=None), Z[:, 50:-50].shape)
    l1 = ind[1] + 50 - 25
    l2 = ind[1] + 50 + 26
    return ind[1]+50, Z_copy[:, l1:l2], E_copy[:, l1:l2]


def find_source(dx, dy, skysub, commonwave, obj, specn, error,
                xoff, yoff, wave_0, ispec):
    D = np.sqrt((dx - dx[:, np.newaxis])**2 + (dy - dy[:, np.newaxis])**2)
    loc, sdimage, sderror = convolve_spatially(dx, dy, skysub, commonwave,
                                               specn, error, ispec*1.,
                                               sig_wave=1.5)
    BN = int(sdimage.shape[1] / 2)
    kSN = 0.
    for k in np.arange(BN):
        k = int(k)
        dimage = np.sum(sdimage[:, (BN-k):(BN+k+1)], axis=1)
        derror = np.sqrt(np.sum(sderror[:, (BN-k):(BN+k+1)]**2, axis=1))
        sn = dimage * 0.
        for i in np.arange(len(dimage)):
            sel = D[i, :] < 1.5
            S = np.sum(dimage[sel])
            N = np.sqrt(np.sum(derror[sel]**2))
            sn[i] = S / N
        SN = np.nanmax(sn)
        if kSN > SN:
            break
        else: 
            kSN = SN * 1.
    kind = 'Emission'
      
    if args.model_dar:
        D = get_standard_star_params(skysub, commonwave, calinfo[5][:, 0],
                                         calinfo[5][:, 1])
        
        xc, yc, xstd, ystd, xoff, yoff = D
        log.info('%s, %s: Source found at s/n: %0.2f' % (obj, specn, SN))
        return xc, yc, xstd, ystd, xoff, yoff
    if SN > 5.:
        ind = np.argmax(dimage)
        dist = np.sqrt((dx - dx[ind])**2 + (dy - dy[ind])**2)
        inds = dist < 1.5
        x_centroid = np.sum(dimage[inds] * dx[inds]) / np.sum(dimage[inds])
        y_centroid = np.sum(dimage[inds] * dy[inds]) / np.sum(dimage[inds])
        X = np.ones(commonwave.shape)
        G = Gaussian2D()
        fitter = LevMarLSQFitter()
        G.amplitude.value = dimage[ind]
        G.x_mean.value = x_centroid
        G.y_mean.value = y_centroid
        fit = fitter(G, dx[inds], dy[inds], dimage[inds])
        seeing = 2.35 * np.sqrt(fit.x_stddev * fit.y_stddev)
        if seeing < 0.75:
            log.info('%s, %s: %s source found at s/n: %0.2f but rejected for being too small, col: %i' % (obj, specn, kind, SN, loc))
            return None
        else:
            log.info('%s, %s: %s source found at s/n: %0.2f, with fwhm: %0.2f' % (obj, specn, kind, SN, seeing))

        xoff = xoff - xoff[loc]
        yoff = yoff - yoff[loc]
        return x_centroid, y_centroid, fit.x_stddev.value*X, fit.y_stddev.value * X, xoff, yoff
    else:
        log.info('%s, %s: No source found, s/n too low: %0.2f' % (obj, specn, SN))
        return None


def get_standard_star_params(data, commonwave, xloc, yloc):
    G = Gaussian2D()
    fitter = LevMarLSQFitter()
    wchunk = np.array([np.mean(chunk)
                       for chunk in np.array_split(commonwave, 11)])
    dchunk = [np.median(chunk, axis=1)
              for chunk in np.array_split(data, 11, axis=1)]
    xc, yc, xs, ys = [i * wchunk for i in [0., 0., 0., 0.]]
    for i in np.arange(11):
        y = dchunk[i]
        ind = np.argmax(y)
        dist = np.sqrt((xloc - xloc[ind])**2 + (yloc - yloc[ind])**2)
        inds = dist < 3.
        x_centroid = np.sum(y[inds] * xloc[inds]) / np.sum(y[inds])
        y_centroid = np.sum(y[inds] * yloc[inds]) / np.sum(y[inds])
        G.amplitude.value = y[ind]
        G.x_mean.value = x_centroid
        G.y_mean.value = y_centroid
        fit = fitter(G, xloc[inds], yloc[inds], y[inds])
        xc[i] = fit.x_mean.value * 1.
        yc[i] = fit.y_mean.value * 1.
        xs[i] = fit.x_stddev.value * 1.
        ys[i] = fit.y_stddev.value * 1.
    sel = xs > 0.
    xoff = np.polyval(np.polyfit(wchunk[sel], xc[sel], 2), commonwave)
    yoff = np.polyval(np.polyfit(wchunk[sel], yc[sel], 2), commonwave)
    xstd = np.median(xs[sel]) * np.ones(commonwave.shape)
    ystd = np.median(ys[sel]) * np.ones(commonwave.shape)
    N = len(commonwave)
    return xoff[N/2], yoff[N/2], xstd, ystd, xoff - xoff[N/2], yoff - yoff[N/2]


def get_bigarray(xloc, yloc):
    BigX = [xloc]
    BigY = [yloc]
    uy = np.unique(yloc)
    dy = np.mean(np.diff(uy))
    for i in np.arange(1, 10):
        ny = dy + np.max(np.hstack(BigY))
        x = np.hstack(BigX)[np.where(uy[-2] == np.hstack(BigY))[0]]
        y = ny * np.ones(x.shape)
        BigY.append(y)
        BigX.append(x)
        uy = np.unique(np.hstack(BigY))
        ny = -dy + np.min(np.hstack(BigY))
        x = np.hstack(BigX)[np.where(uy[1] == np.hstack(BigY))[0]]
        y = ny * np.ones(x.shape) 
        BigY.append(y)
        BigX.append(x)
        uy = np.unique(np.hstack(BigY))
    BigX = np.hstack(BigX)
    BigY = np.hstack(BigY)
    uy = np.unique(BigY)
    NX, NY = ([BigX], [BigY])
    for i in uy:
        sel = np.where(i == BigY)[0]
        dx = np.abs(np.mean(np.diff(BigX[sel])))
        xn = np.min(BigX[sel]) - np.arange(1, 10)*dx
        xn2 = np.max(BigX[sel]) + np.arange(1, 10)*dx
        yn = i * np.ones(xn.shape)
        yn2 = i * np.ones(xn2.shape)
        NX.append(xn)
        NX.append(xn2)
        NY.append(yn)
        NY.append(yn2)
    BigX = np.hstack(NX)
    BigY = np.hstack(NY)
    M = np.zeros(BigX.shape, dtype=bool)
    inds = np.zeros(xloc.shape, dtype=int)
    for i in np.arange(len(xloc)):
        ind = np.where((xloc[i] == BigX) * (yloc[i] == BigY))[0]
        M[ind] = True
        inds[i] = ind
    return BigX, BigY


def extract_source(data, xc, yc, xoff, yoff, wave, xloc, yloc, error,
                   xstd, ystd):
    PSF = Gaussian2D()
    spec = wave * 0.
    serror = wave * 0.
    weights = data * 0.
    mask = data * 0.
    bigX, bigY = get_bigarray(xloc, yloc)
    for i in np.arange(len(wave)):
        x = xc + xoff[i]
        y = yc + yoff[i]
        PSF.x_mean.value = x
        PSF.y_mean.value = y
        PSF.x_stddev = xstd[i]
        PSF.y_stddev = ystd[i]
        seeing = np.sqrt(ystd[i]*xstd[i])
        W = PSF(xloc, yloc)
        S = PSF(bigX, bigY).sum()
        W /= S
        M = (error[:, i] != 0.) * (np.sqrt((xloc-x)**2 + (yloc-y)**2) < seeing*2.5)
        weights[:, i] = W
        mask[:, i] = M
        spec[i] = (data[:, i] * W * M).sum() / (M * W**2).sum()
        serror[i] = np.sqrt((error[:, i]**2 * W * M).sum() / (M * W**2).sum())
    return spec, serror, weights, mask


def get_mirror_illumination(fn=None, default=51.4e4):
    ''' Use hetillum from illum_lib to calculate mirror illumination (cm^2) '''
    #log.info('Getting mirror illumination')
    try:
        F = fits.open(fn)
        names = ['RHO_STRT', 'THE_STRT', 'PHI_STRT', 'X_STRT', 'Y_STRT']
        r, t, p, x, y = [F[0].header[name] for name in names]
        # log.info('Rho, Theta, Phi, X, Y: %0.4f, %0.4f, %0.4f, %0.4f, %0.4f' %
        #          (r, t, p, x, y))
        mirror_illum = float(os.popen('/work/03730/gregz/maverick/illum_lib/hetillum -p'
                             ' -x "[%0.4f,%0.4f,%0.4f]" "[%0.4f,%0.4f]" 256' %
                                      (x, y, p, 0.042, 0.014)).read().split('\n')[0])
        area = mirror_illum * default
    except:
        log.info('Using default mirror illumination value')
        area = default
    # log.info('Mirror illumination: %0.2f m^2' % (area/1e4))
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

def truncate_list(lst):
    if len(lst) <= 5:
        return lst  # If the list already has 5 or fewer elements, return it as is.

    # Calculate indices for the three evenly spaced elements in the middle.
    step = (len(lst) - 1) / 4
    indices = [0, round(step), round(2 * step), round(3 * step), len(lst) - 1]

    # Select elements at the calculated indices
    return [lst[i] for i in indices]

def get_mirror_illumination_guider(fn, exptime, default=51.4e4,
                                   path='/work/03946/hetdex/maverick'):
    try:
        M = []
        path = op.join(path, args.date)
        f = op.basename(fn)
        DT = f.split('_')[0]
        y, m, d, h, mi, s = [int(x) for x in [DT[:4], DT[4:6], DT[6:8], DT[9:11],
                             DT[11:13], DT[13:15]]]
        d0 = datetime(y, m, d, h, mi, s)
        tarfolder = op.join(path, 'gc1', '*.tar')
        tarfolder = glob.glob(tarfolder)
        if len(tarfolder) == 0:
            area = 51.4e4
            log.info('Using default mirror illumination: %0.2f m^2' % (area/1e4))
            return area
        T = tarfile.open(tarfolder[0], 'r')
        init_list = sorted([name for name in T.getnames()
                            if name[-5:] == '.fits'])
        final_list = []
        for t in init_list:
            DT = op.basename(t).split('_')[0]
            y, m, d, h, mi, s = [int(x) for x in [DT[:4], DT[4:6], DT[6:8],
                                 DT[9:11], DT[11:13], DT[13:15]]]
            d = datetime(y, m, d, h, mi, s)
            p = (d - d0).seconds
            if (p > -10.) * (p < exptime+10.):
                final_list.append(t)
        final_list = truncate_list(final_list)
        for fn in final_list:
            fobj = T.extractfile(T.getmember(fn))
            M.append(get_mirror_illumination(fobj))
        M = np.array(M)
        sel = M != 51.4e4
        if sel.sum() > 0.:
            area = np.mean(M[sel])
        else:
            area = 51.4e4
        log.info('Final Mirror illumination: %0.2f m^2' % (area/1e4))
        return area
    except: 
        log.info('Using default mirror illumination: %0.2f m^2' % (default/1e4))
        return default


def get_throughput(fn, exptime, path='/work/03946/hetdex/maverick'):
    attr = ['GUIDLOOP', 'MJD', 'TRANSPAR']
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
            DT = op.basename(t).split('_')[0]
            y, m, d, h, mi, s = [int(x) for x in [DT[:4], DT[4:6], DT[6:8],
                                 DT[9:11], DT[11:13], DT[13:15]]]
            d = datetime(y, m, d, h, mi, s)
            p = (d - d0).seconds
            if (p > -10.) * (p < exptime+10.):
                final_list.append(t)
        final_list = truncate_list(final_list)
        for fnt in final_list:
            fobj = T.extractfile(T.getmember(fnt))
            f = fits.open(fobj)
            if f[1].header['GUIDLOOP'] == 'ACTIVE':
                M.append([])
                for att in attr:
                    M[-1].append(f[1].header[att])
    throughput = np.zeros((len(M),))
    for i, mi in enumerate(M):
        if len(mi) < 2:
            continue
        if mi[2] > 0.:
            throughput[i] = mi[2]
    t = np.mean(throughput[throughput>0.0])
    if np.isnan(t):
        log.warning('Could not find TRANSPAR measurments')
        t = 1.0
    if t > 1.1:
        log.info('The throughput for %s is %0.2f is too high, setting to 1.0' % (fn, t))
        t = 1.0
    if t < 0.1:
        log.info('The throughput for %s is %0.2f is too high, setting to 1.0' % (fn, t))
        t = 1.0
    log.info('Throughput for %s is %0.2f' % (fn, t))
    return t


def check_if_standard(objname):
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
    #scifiles = op.join(op.dirname(bf.replace('exp01', 'exp*')), '*%sLL*.fits' % ifuslot)
    scifiles = op.join(op.dirname(bf), '*%sLL*.fits' % ifuslot)
    images, rect, spec, cos, fl, Fi, E, header = extract_sci(scifiles, amps, calinfo[2],
                                              calinfo[1], calinfo[0], calinfo[3],
                                              calinfo[4], calinfo[5])
    cnt = 1
    if args.central_wave is None:
        wave_0 = np.mean(commonwave)
        wb = 50.
    else:
        wave_0 = args.central_wave
        wb = args.wavelength_bin
    darfile = op.join(DIRNAME, 'lrs2_config/dar_%s.dat' % specinit)
    T = Table.read(darfile, format='ascii.fixed_width_two_line')
    xoff = (np.interp(commonwave, T['wave'], T['x_0']) -
            np.interp(wave_0, T['wave'], T['x_0']))
    yoff = (np.interp(commonwave, T['wave'], T['y_0']) -
            np.interp(wave_0, T['wave'], T['y_0']))
    for im, r, s, c, fli, Fii, e, he in zip(images, rect, spec, cos,
                                            fl, Fi, E, header):
#        try:
        try:
            basename = 'LRS2/' + he['QPROG']
        except:
            if check_if_standard(obj[0]) and (ifuslot in obj[0]):
                basename = 'LRS2/STANDARDS'
            else:
                basename = 'LRS2/ORPHANS'
        mkpath(basename)
        
        pos = np.zeros((len(calinfo[5]), 6))
        pos[:, 0:2] = calinfo[5]
        try:
            PA = float(he['PARANGLE'])
            RA = float(he['TRAJRA'])
            DEC = float(he['TRAJDEC'])
            log.info('Observation at %0.4f %0.4f, PA: %0.3f' % (RA, DEC, PA))
            A = Astrometry(RA, DEC, PA, 0., 0., fplane_file=fplane_file)
            ra, dec = A.get_ifupos_ra_dec(ifuslot, calinfo[5][:, 0],
                                          calinfo[5][:, 1])

            fpx = A.fplane.by_ifuslot(ifuslot).y + calinfo[5][:, 0]
            fpy = A.fplane.by_ifuslot(ifuslot).x + calinfo[5][:, 1]

            pos[:, 2] = fpx
            pos[:, 3] = fpy
            pos[:, 4] = ra
            pos[:, 5] = dec
        except:
            log.warning('Astrometry Issue')
        fn = op.join(op.dirname(bf.replace('exp01', 'exp%02d' % cnt)), '*%sLL*.fits' % ifuslot)
        fn = get_filenames_from_tarfolder(get_tarname_from_filename(fn),
                                      fn)
        mini = get_objects(fn, ['OBJECT', 'EXPTIME'], full=True)
        log.info('Subtracting sky %s, exp%02d' % (obj[0], cnt))
        r[calinfo[6][:, 1] == 1.] = 0.
        e[calinfo[6][:, 1] == 1.] = 0.

        r /= mini[0][1]
        r /= mini[0][2]
        r /= mini[0][3]
        e /= mini[0][1]
        e /= mini[0][2]
        e /= mini[0][3]
        if args.correct_ftf:
            r, e = correct_ftf(r, e)
        bad = get_all_cosmics(pos[:, 0], pos[:, 1], r*1., e * 1.)
        e[bad] = 0.
        sky = sky_subtraction(r, e, pos[:, 0], pos[:, 1])
        sky[calinfo[6][:, 1] == 1.] = 0.
        skysub = r - sky
        if response is not None:
            r *= response
            e *= response
            sky *= response
            skysub *= response
        X = np.array([T['wave'], T['x_0'], T['y_0']])
        for S, name in zip([r, sky, skysub], ['obs', 'sky', 'skysub']):
            outname = ('%s_%s_%s_%s_%s_cube.fits' % (args.date, sci_obs,
                       'exp%02d' % cnt, specname, name))
            outname = op.join(basename, outname)
            
            zcube, zimage, xgrid, ygrid = make_frame(calinfo[5][:, 0],
                                                     calinfo[5][:, 1], S, e,
                                                     commonwave,
                                                     T['wave'],
                                                     xoff, yoff,
                                                     wstart=wave_0-wb,
                                                     wend=wave_0+wb)
            write_cube(commonwave, xgrid, ygrid, zcube, outname, he)
        loc = None
        if args.source_x is None:
            wi = np.searchsorted(commonwave, wave_0-wb, side='left')
            we = np.searchsorted(commonwave, wave_0+wb, side='right')
            dimage = np.median(skysub[:, wi:we+1], axis=1)
            derror = np.sqrt(np.sum(e[:, wi:we+1]**2, axis=1))*1.253 / np.sqrt(we-wi+1)
            loc1 = find_source(pos[:, 0], pos[:, 1],
                               skysub, commonwave, obj[0], specname, e,
                               xoff, yoff, wave_0, r)
            if loc1 is not None:
                loc = [0., 0., 0.]
                loc[0] = loc1[0]
                loc[1] = loc1[1]                
                loc[2] = 2.35 * np.mean(np.sqrt(loc1[2]*loc1[3]))
                xstd = loc1[2]
                ystd = loc1[3]
                xoff = loc1[4]
                yoff = loc1[4]
                log.info('Source seeing initially found to be: %0.2f' % loc[2])
        if args.source_x is not None:
            loc = [args.source_x, args.source_y, 1.5]
            xstd = np.ones(commonwave.shape) * 0.75
            ystd = np.ones(commonwave.shape) * 0.75
            
        if loc is not None:
            log.info('Source found at %0.2f, %0.2f' % (loc[0], loc[1]))
            skyspec, errorskyspec, w, m = extract_source(sky, loc[0], loc[1], xoff,
                                                   yoff, commonwave,
                                                   calinfo[5][:, 0],
                                                   calinfo[5][:, 1], e,
                                                   xstd, ystd)
            skysubspec, errorskysubspec, w, m = extract_source(skysub, loc[0],
                                                         loc[1], xoff, yoff,
                                                         commonwave,
                                                         calinfo[5][:, 0],
                                                         calinfo[5][:, 1], e,
                                                         xstd, ystd)
        else:
            skyspec = commonwave * 0.
            skysubspec = commonwave * 0.
            errorskyspec = commonwave * 0.
            errorskysubspec = commonwave * 0.
        if response is not None:
            f5 = np.vstack([commonwave, skysubspec, skyspec,
                            errorskysubspec, errorskyspec,
                            response])
        else:
            f5 = np.vstack([commonwave, skysubspec, skyspec,
                            errorskysubspec, errorskyspec,
                            np.ones(commonwave.shape)])

        f1 = create_header_objection(commonwave, r, func=fits.PrimaryHDU)
        f2 = create_header_objection(commonwave, sky)
        f3 = create_header_objection(commonwave, skysub)
        f4 = create_image_header(commonwave, xgrid, ygrid, zimage)
        f6 = create_header_objection(commonwave, e)
        for key in he.keys():
            if key in f1.header:
                continue
            if 'SEC' in key:
                continue
            if ('BSCALE' in key) or ('BZERO' in key):
                continue
            try:
                f1.header[key] = he[key]
            except:
                continue
        if loc is not None:
            f1.header['SOURCEX'] = loc[0]
            f1.header['SOURCEY'] = loc[1]
            f1.header['SEEING'] = loc[2]
            
        f1.header['MILLUM'] = mini[0][2]
        f1.header['THROUGHP'] = mini[0][3]
        names = ['observed_spectra', 'sky_spectra', 'skysub_spectra',
                 'error_spectra', 'collapsed_image', 'fiber_positions',
                 'extracted_spectrum', 'adr', 'bigw', 'image',
                 'flattened_image', 'trace', 'cosmics', 'unrectified_spectra']
        flist = [f1, f2, f3, f6, f4, fits.ImageHDU(pos), fits.ImageHDU(f5),
                 fits.ImageHDU(X), fits.ImageHDU(calinfo[3]),
                 fits.ImageHDU(im), fits.ImageHDU(fli), fits.ImageHDU(Fii),
                 fits.ImageHDU(c), fits.ImageHDU(s)]
        for fl, name in zip(flist, names):
            fl.header['EXTNAME'] = name
        outname = ('%s_%s_%s_%s_%s.fits' % ('multi', args.date, sci_obs,
                                            'exp%02d' % cnt, specname))
        outname = op.join(basename, outname)
        fits.HDUList(flist).writeto(outname, overwrite=True)
        outname = ('%s_%s_%s_%s_%s.fits' % ('spectrum', args.date, sci_obs,
                                            'exp%02d' % cnt, specname))
        outname = op.join(basename, outname)
        f1 = fits.PrimaryHDU(f5)
        for key in he.keys():
            if key in f1.header:
                continue
            if 'SEC' in key:
                continue
            if ('BSCALE' in key) or ('BZERO' in key):
                continue
            try:
                f1.header[key] = he[key]
            except:
                continue
        names = ['wavelength', 'F_lambda', 'Sky_lambda', 'e_F_lambda',
                 'e_Sky_lambda', 'response']
        f1.header['DWAVE'] = commonwave[1] - commonwave[0]
        f1.header['WAVE0'] = commonwave[0]
        f1.header['WAVESOL'] = 'WAVE0 + DWAVE * linspace(0, NAXIS1)'
        f1.header['WAVEUNIT'] = 'A'
        if loc is not None:
            f1.header['SOURCEX'] = loc[0]
            f1.header['SOURCEY'] = loc[1]
            f1.header['SEEING'] = loc[2]
            f1.header['MILLUM'] = mini[0][2]
            f1.header['THROUGHP'] = mini[0][3]

        if response is not None:
            f1.header['FLUXUNIT'] = 'ergs/s/cm2/A'
        else:
            f1.header['FLUXUNIT'] = 'e-/s/cm2/A'
        for i, name in enumerate(names):
            f1.header['ROW%i' % (i+1)] = name
        f1.writeto(outname, overwrite=True)
        if standard and ((skysubspec != 0.).sum() > 500):
            return get_response(obj[0], commonwave, skysubspec, specname)
        cnt += 1
#        except Exception as e: 
#            log.warning('Exposure %i Failed' % cnt)
#            cnt += 1

def get_filenames_from_tarfolder(tarfolder, path):
    T = tarfile.open(tarfolder, 'r')
    names = T.getnames()
    matches = fnmatch.filter(names, op.join('*', op.basename(path)))
    matches = [op.join(op.dirname(tarfolder), match) for match in matches]
    matches = sorted(matches)
    T.close()
    return matches

def get_cal_path(pathname, date, ndays=31):
    date_ = datetime(int(date[:4]), int(date[4:6]), int(date[6:]))
    # Controller change date
    date_controller_swap = datetime(2024, 7, 22)
    if date_ > date_controller_swap:
        flag_new = True
    else:
        flag_new = False
    filenames = []
    while len(filenames) == 0:
        datel = date_ - timedelta(days=int(ndays/2))
        for i in np.arange(ndays):
            ndate = datel + timedelta(days=int(i))
            if flag_new and (ndate <= date_controller_swap):
                continue
            daten = '%04d%02d%02d' % (ndate.year, ndate.month, ndate.day)
            npath = pathname.replace(date, daten)
            tarpath = get_tarname_from_filename(npath)
            for tarname in sorted(glob.glob(tarpath)):
                filenames.append(get_filenames_from_tarfolder(tarname, npath))
        flat_list = [item for sublist in filenames for item in sublist]
        filenames = sorted(flat_list)
        ndays += 1
    return filenames

def get_previous_night(daten):    
    datec_ = datetime(int(daten[:4]), int(daten[4:6]), int(daten[6:]))
    daten_ = datec_ - timedelta(days=1)
    daten = '%04d%02d%02d' % (daten_.year, daten_.month, daten_.day)
    return daten


# LRS2-R
fiberpos, fiberspec = ([], [])
log.info('Beginning the long haul.')
allflatspec, allspec, allra, alldec, allx, ally, allsub = ([], [], [], [], [],
                                                           [], [])

DIRNAME = get_script_path()

for info in listinfo:
    specinit, specname, multi, lims, amps, slims, arc_names = info
    if int(args.date) < 20161101:
        nnn = specname #'%s_old' % specname
    else:
        nnn = specname
    arc_lines = Table.read(op.join(DIRNAME, 'lrs2_config/lines_%s.dat' %
                                   nnn), format='ascii')
    commonwave = np.linspace(lims[0], lims[1], 2064)
    specid, ifuslot, ifuid = multi.split('_')
    package = []
    marc, mtwi, mflt = ([], [], [])
    twipath = twi_path % (instrument, instrument, '00000*', instrument,
                                ifuslot)
    twifiles = get_cal_path(twipath, args.date, ndays=15)
    flt_path = (twi_path.replace('twi', 'flt') %
                    (instrument, instrument, '00000*', instrument, ifuslot))
    fltfiles = get_cal_path(flt_path, args.date, ndays=15)
    if args.use_flat:
        twifiles = fltfiles
    for amp in amps:
        amppos = get_ifucenfile(specname, amp)
        ##############
        # MASTERBIAS #
        ##############
        log.info('Getting Masterbias for ifuslot, %s, and amp, %s' %
                 (ifuslot, amp))
        zro_path = bias_path % (instrument, instrument, '00000*', instrument,
                                ifuslot)
        zrofiles = get_cal_path(zro_path, args.date, ndays=2)
        masterbias = get_masterbias(zrofiles, amp)

        #####################
        # MASTERTWI [TRACE] #
        #####################
        log.info('Getting MasterFlat for ifuslot, %s, and amp, %s' %
                 (ifuslot, amp))
        log.info('Getting Trace for ifuslot, %s, and amp, %s' %
                 (ifuslot, amp))
        log.info('Number of twi files: %i' % len(twifiles))

        masterflt = get_mastertwi(twifiles, amp, masterbias)
        trace, dead = get_trace(masterflt, specid, ifuslot, ifuid, amp,
                                args.date)
        

        ##########################
        # MASTERARC [WAVELENGTH] #
        ##########################
        log.info('Getting MasterArc for ifuslot, %s, and amp, %s' %
                 (ifuslot, amp))
        lamp_path = cmp_path % (instrument, instrument, '00000*', instrument,
                                ifuslot)
        lampfiles = get_cal_path(lamp_path, args.date, ndays=15)
        log.info('Number of arc files: %i' % len(lampfiles))
        masterarc = get_masterarc(lampfiles, amp, arc_names, masterbias,
                                  specname, trace)
        
        lampfiles = get_cal_path(lamp_path.replace(args.date, '20181201'),
                                 '20181201', ndays=3)
        def_arc = get_masterarc(lampfiles, amp,
                                arc_names, masterbias, specname, trace)

        #fits.PrimaryHDU(masterarc).writeto('/work/03946/hetdex/maverick/run_lrs2/wtf_%s_%s.fits' % (ifuslot, amp), overwrite=True)
        log.info('Getting Wavelength for ifuslot, %s, and amp, %s' %
                 (ifuslot, amp))
        
        wave = get_wavelength_from_arc(masterarc, trace, arc_lines, specname,
                                       amp, int(args.date), otherimage=def_arc)
        #fits.PrimaryHDU(wave).writeto('test_wave.fits', overwrite=True)

        #################################
        # TWILIGHT FLAT [FIBER PROFILE] #
        #################################
        log.info('Getting bigW for ifuslot, %s, and amp, %s' %
                 (ifuslot, amp))
        bigW = get_bigW(amp, wave, trace, masterbias)
        package.append([wave, trace, bigW, masterbias, amppos, dead])
        log.info('Number of flt files: %i' % len(fltfiles))

        masterFlat = get_mastertwi(fltfiles, amp, masterbias)

        marc.append(masterarc)
        mtwi.append(masterflt)
        mflt.append(masterFlat)
    # Normalize the two amps and correct the flat
    calinfo = [np.vstack([package[0][i], package[1][i]])
               for i in np.arange(len(package[0]))]
    masterarc, masterflt, masterFlat = [np.vstack(x) for x in [marc, mtwi, mflt]]
    calinfo[1][package[0][1].shape[0]:, :] += package[0][2].shape[0]
    log.info('Getting flat for ifuslot, %s, side, %s' % (ifuslot, specname))
    twiflat = get_twiflat_field(twifiles, amps, calinfo[0], calinfo[1],
                                calinfo[2], commonwave, calinfo[3], specname)
    calinfo.insert(2, twiflat)
    flatspec = get_spectra(calinfo[2], calinfo[1])
    for mfile in [masterarc, masterflt, masterFlat]:
        masterarcerror = np.sqrt(3.**2 + np.where(mfile > 0., mfile, 0.))
        arcspec, ae, Cc, Yyy, Fff = weighted_extraction(mfile, masterarcerror,
                                                        calinfo[2], calinfo[1],
                                                        cthresh=500)
        sP = np.zeros((calinfo[0].shape[0], len(commonwave)))
        for fiber in np.arange(calinfo[0].shape[0]):
            I = interp1d(calinfo[0][fiber], arcspec[fiber],
                                 kind='linear', fill_value='extrapolate')
            sP[fiber] = I(commonwave)
        calinfo.append(sP)
    bigF = get_bigF(calinfo[1], calinfo[2])
    calinfo.append(bigF)
    #####################
    # SCIENCE REDUCTION #
    #####################
    response = None
    pathS = sci_path % (instrument, instrument, '0000*',
                                             '01', instrument, ifuslot)
    basefiles = []
    for tarname in glob.glob(get_tarname_from_filename(pathS)):
        basefiles.append(get_filenames_from_tarfolder(tarname, pathS))
    flat_list = [item for sublist in basefiles for item in sublist]
    basefiles = [f for f in sorted(flat_list) if "exp01" in f]

    all_sci_obs = [op.basename(op.dirname(op.dirname(op.dirname(fn))))[-7:]
                   for fn in basefiles]
    objects = get_objects(basefiles, ['OBJECT', 'EXPTIME'])
    if response is None:
        log.info('Getting average response')
        basename = 'LRS2/CALS'
        R = fits.open(op.join(DIRNAME,
                              'lrs2_config/response_%s.fits' % specname))
        response = R[0].data[1]*1.

    f = []
    names = ['wavelength', 'trace', 'flat', 'bigW', 'masterbias',
             'xypos', 'dead', 'arcspec', 'fltspec', 'Flatspec', 'bigF']
    for i, cal in enumerate(calinfo):
        if i == 0:
            func = fits.PrimaryHDU
        else:
            func = fits.ImageHDU
        f.append(func(cal))
    f.append(fits.ImageHDU(masterarc))
    names.append('masterarc')
    f.append(fits.ImageHDU(masterFlat))
    names.append('masterFlat')
    if response is not None:
        f.append(fits.ImageHDU(np.array([commonwave, response], dtype=float)))
        names.append('response')
    for fi, n in zip(f, names):
        fi.header['EXTNAME'] = n
    basename = 'LRS2/CALS'
    mkpath(basename)
    fits.HDUList(f).writeto(op.join(basename,
                            'cal_%s_%s.fits' % (args.date, specname)),
                            overwrite=True)
    for sci_obs, obj, bf in zip(all_sci_obs, objects, basefiles):
        log.info('Checkpoint --- Working on %s, %s' % (bf, specname))
        if args.object is None:
            big_reduction(obj, bf, instrument, sci_obs, calinfo, amps, commonwave,
                          ifuslot, specname, response=response)
        else:
            if args.object.lower() in obj[0].lower():
                big_reduction(obj, bf, instrument, sci_obs, calinfo, amps, commonwave,
                          ifuslot, specname, response=response)
            if check_if_standard(obj[0]) and (ifuslot in obj[0]):
                big_reduction(obj, bf, instrument, sci_obs, calinfo, amps, commonwave,
                          ifuslot, specname, response=response)
