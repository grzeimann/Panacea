# -*- coding: utf-8 -*-
"""
Created on Thu Oct  4 13:30:06 2018

@author: gregz
"""
import time
import numpy as np
import os.path as op
import glob
import warnings

from astropy.io import fits
from utils import biweight_location
from scipy.signal import savgol_filter, medfilt2d
from scipy.interpolate import interp1d, interp2d
from input_utils import setup_logging
from astrometry import Astrometry

dither_pattern = np.array([[0., 0.], [1.27, -0.73], [1.27, 0.73]])
virus_amps = ['LL', 'LU', 'RU', 'RL']
lrs2_amps = [['LL', 'LU'], ['RL', 'RU']]
fplane_file = '/work/03730/gregz/maverick/fplane.txt'
flt_obs = '%07d' % 15
twi_obs = '%07d' % 1
sci_obs = '%07d' % 13
twi_date = '20170205'
sci_date = twi_date
flt_date = twi_date

# FOR LRS2
instrument = 'lrs2'
AMPS = lrs2_amps[0]
dither_pattern = np.zeros((10, 2))

log = setup_logging('panacea_quicklook')

basered = '/work/03730/gregz/maverick'
baseraw = '/work/03946/hetdex/maverick'

twi_path = op.join(basered, 'reductions', twi_date, '%s', '%s%s', 'exp01',
                   '%s', 'multi_*_%s_*_LL.fits')

sci_path = op.join(baseraw, sci_date,  '%s', '%s%s', 'exp%s',
                   '%s', '2*_%sLL*.fits')
flt_path = op.join(baseraw, flt_date,  '%s', '%s%s', 'exp*',
                   '%s', '2*_%sLL*.fits')
sciflt_path = op.join(baseraw, twi_date,  '%s', '%s%s', 'exp*',
                      '%s', '2*_%sLL_twi.fits')
bias_path = op.join(baseraw, twi_date, '%s', '%s%s', 'exp*',
                    '%s', '2*_%sLL_zro.fits')


def get_cal_info(twi_path, amp):
    F = fits.open(glob.glob(twi_path.replace('LL', amp))[0])
    return (np.array(F['ifupos'].data, dtype=float),
            np.array(F['trace'].data, dtype=float),
            np.array(F['wavelength'].data, dtype=float))


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
    try:
        ampname = a[0].header['AMPNAME']
    except:
        ampname = None
    a = orient_image(image, amp, ampname)
    return a


def get_sciflat_field(flt_path, amp, array_wave, array_trace, common_wave,
                      masterbias, log):
    files = glob.glob(flt_path.replace('LL', amp))
    listflat = []
    array_flt = base_reduction(files[0])

    bigW = np.zeros(array_flt.shape)
    Y, X = np.indices(array_wave.shape)
    YY, XX = np.indices(array_flt.shape)
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
    listspec = []
    for filename in files:
        log.info('Working on sciflat %s' % filename)
        array_flt = base_reduction(filename) - masterbias
        x = np.arange(array_wave.shape[1])
        spectrum = array_trace * 0.
        for fiber in np.arange(array_wave.shape[0]):
            indl = np.floor(array_trace[fiber]).astype(int)
            indh = np.ceil(array_trace[fiber]).astype(int)
            spectrum[fiber] = array_flt[indl, x] / 2. + array_flt[indh, x] / 2.
        smooth = savgol_filter(spectrum, 315, 1, axis=1)
        avg = biweight_location(smooth, axis=(0,))
        norm = biweight_location(smooth / avg, axis=(1,))
        nw, ns = make_avg_spec(array_wave, spectrum / norm[:, np.newaxis],
                               binsize=41, per=95)
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
        listflat.append(flat)
        listspec.append(I(common_wave))
    flat = np.median(listflat, axis=(0,))
    flat[~np.isfinite(flat)] = 0.0
    flat[flat < 0.0] = 0.0

    return flat, bigW, np.nanmedian(listspec, axis=0)


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


def find_cosmics(Y, E, thresh=8.):
    A = medfilt2d(Y, (5, 1))
    S = safe_division((Y - A), E)
    P = S - medfilt2d(S, (1, 15))
    x, y = np.where(P > thresh)
    xx, yy = ([], [])
    for i in np.arange(-1, 2):
        for j in np.arange(-1, 2):
            sel = ((x + i) >= 0) * ((x + i) < Y.shape[0])
            sel2 = ((y + j) >= 0) * ((y + j) < Y.shape[1])
            sel = sel * sel2
            xx.append((x + i)[sel])
            yy.append((y + j)[sel])
    xx = np.hstack(xx)
    yy = np.hstack(yy)
    inds = np.ravel_multi_index([xx, yy], Y.shape)
    inds = np.unique(inds)
    C = np.zeros(Y.shape, dtype=bool).ravel()
    C[inds] = True
    C = C.reshape(Y.shape)
    log.info('Number of pixels affected by cosmics: %i' % len(x))
    log.info('Fraction of pixels affected by cosmics: %0.5f' %
             (1.*len(inds)/Y.shape[0]/Y.shape[1]))
    return C


def weighted_extraction(image, flat, trace):
    gain = 0.83
    rdnoise = 3.
    I = image * 1.
    I[I < 0.] = 0.
    E = np.sqrt(rdnoise**2 + gain * I) / gain
    E = safe_division(E, flat)
    E[E < 1e-8] = 1e9
    Y = safe_division(image, flat)
    cosmics = find_cosmics(Y, E)
    x = np.arange(trace.shape[1])
    spectrum = 0. * trace
    TT = np.zeros((trace.shape[0], 3, trace.shape[1], 4))
    for fiber in np.arange(trace.shape[0]):
        T = np.zeros((3, trace.shape[1], 4))
        indl = np.floor(trace[fiber]).astype(int)
        flag = False
        for ss, k in enumerate(np.arange(-1, 3)):
            try:
                T[0, :, ss] = Y[indl+k, x]
                T[1, :, ss] = 1. / E[indl+k, x]**2
                T[2, :, ss] = ~cosmics[indl+k, x]
            except:
                v = indl+k
                sel = np.where((v >= 0) * (v < Y.shape[0]))[0]
                T[0, sel, ss] = Y[v[sel], x[sel]]
                T[1, sel, ss] = 1. / E[v[sel], x[sel]]**2
                T[2, sel, ss] = ~cosmics[v[sel], x[sel]]
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
        else:
            a = np.sum(T[0] * T[1] * T[2], axis=1)
            b = np.sum(T[1] * T[2], axis=1)
            spectrum[fiber] = safe_division(a, b)
        TT[fiber] = T
    fits.PrimaryHDU(TT).writeto('wtf2.fits', overwrite=True)
    return spectrum


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
    shifts = np.polyval(np.polyfit(np.nanmedian(FlatTrace, axis=1), shifts, 1),
                        Yx)
    fits.HDUList([fits.PrimaryHDU(FlatTrace),
                  fits.ImageHDU(Trace)]).writeto('test_trace.fits',
                                                 overwrite=True)
    return shifts


def subtract_sci(sci_path, flat, array_trace, array_wave, bigW, masterbias):
    files = sorted(glob.glob(sci_path.replace('LL', amp)))
    array_list = []
    for filename in files:
        log.info('Skysubtracting sci %s' % filename)
        array_flt = base_reduction(filename) - masterbias
        array_list.append(array_flt)
    sci_array = np.sum(array_list, axis=0)
    Xx = np.arange(flat.shape[1])
    Yx = np.arange(flat.shape[0])
    I = interp2d(Xx, Yx, flat, kind='cubic', bounds_error=False,
                 fill_value=0.0)
    shifts = get_trace_shift(sci_array, flat, array_trace, Yx)
    flat = I(Xx, Yx + shifts)
    log.info('Found shift for %s of %0.3f' % (files[0], np.median(shifts)))
    array_list = []
    residual = []
    spec_list = []
    for filename in files:
        log.info('Skysubtracting sci %s' % filename)
        array_flt = base_reduction(filename) - masterbias
        array_list.append(array_flt)
        spectrum = weighted_extraction(array_flt, flat, array_trace)
        spectrum[~np.isfinite(spectrum)] = 0.0
        nw, ns = make_avg_spec(array_wave, spectrum, binsize=41)
        ns[~np.isfinite(ns)] = 0.0
        I = interp1d(nw, ns, kind='quadratic', fill_value='extrapolate')
        modelimage = I(bigW)
        residual.append((array_flt - modelimage*flat))
        speclist = []
        for fiber in np.arange(array_wave.shape[0]):
            dlam = np.diff(array_wave[fiber])
            dlam = np.hstack([dlam[0], dlam])
            I = interp1d(array_wave[fiber], spectrum[fiber] / dlam,
                         kind='quadratic', fill_value='extrapolate')
            speclist.append(I(commonwave))
        spec_list.append(np.array(spectrum))
    return np.array(array_list), np.array(residual), np.array(spec_list)


def get_masterbias(zro_path, amp):
    files = glob.glob(zro_path.replace('LL', amp))
    listzro = []
    for filename in files:
        a = base_reduction(filename)
        listzro.append(a)
    return np.median(listzro, axis=0)


# GET ALL VIRUS IFUSLOTS
twilist = glob.glob(twi_path % (instrument, instrument, twi_obs, instrument,
                                '*'))
ifuslots = [op.basename(x).split('_')[2] for x in twilist]
# LRS2-R
ifuslots = ['066']
fiberpos, fiberspec = ([], [])
log.info('Beginning the long haul.')
nexp = len(glob.glob(sci_path % (instrument, instrument, sci_obs, '*',
                                 instrument, ifuslots[0])))
header = fits.open(glob.glob(sci_path % (instrument, instrument, sci_obs, '01',
                                         instrument,
                                         ifuslots[0]))[0])[0].header
PA = float(header['PARANGLE'])
RA = float(header['TRAJRA'])
DEC = float(header['TRAJDEC'])
log.info('Observation at %0.4f %0.4f, PA: %0.3f' % (RA, DEC, PA))
A = Astrometry(RA, DEC, PA, 0., 0., fplane_file=fplane_file)
allflatspec, allspec, allra, alldec, allx, ally, allsub = ([], [], [], [], [],
                                                           [], [])

# Rectified wavelength
commonwave = np.linspace(6450, 8400, 3000)
N = len(ifuslots) * len(virus_amps)
t1 = time.time()
cnt = 0
cnt2 = 0
breakloop = False
for ifuslot in ifuslots:
    for amp in AMPS:
        log.info('Starting on ifuslot, %s, and amp, %s' % (ifuslot, amp))
        twibase = twi_path % (instrument, instrument, twi_obs, instrument,
                              ifuslot)
        amppos, trace, wave = get_cal_info(twibase, amp)
        if wave.ndim == 1:
            log.info('Insufficient cal data for ifuslot, %s, and amp, %s'
                     % (ifuslot, amp))
            continue
        log.info('Getting Masterbias for ifuslot, %s, and amp, %s' %
                 (ifuslot, amp))
        zro_path = bias_path % (instrument, instrument, '00000*', instrument,
                                ifuslot)
        masterbias = get_masterbias(zro_path, amp)
        twibase = sciflt_path % (instrument, instrument, '00000*', instrument,
                                 ifuslot)
        log.info('Getting SciFlat for ifuslot, %s, and amp, %s' %
                 (ifuslot, amp))
        twiflat, bigW, twispec = get_sciflat_field(twibase, amp, wave, trace,
                                                   commonwave, masterbias, log)
        allflatspec.append(twiflat)
        wave = np.array(wave, dtype=float)
        i1 = []
        scifiles = sci_path % (instrument, instrument, sci_obs, '*',
                               instrument, ifuslot)
        images, subimages, spec = subtract_sci(scifiles, twiflat, trace, wave,
                                               bigW, masterbias)
        allsub.append(images)
        allspec.append(spec)
        for i in np.arange(nexp):
            log.info('Getting RA, Dec for exposure, %i,  ifuslot, %s, and amp,'
                     ' %s' % (i+1, ifuslot, amp))
            ra, dec = A.get_ifupos_ra_dec(ifuslot,
                                          amppos[:, 0] + dither_pattern[i, 0],
                                          amppos[:, 1] + dither_pattern[i, 1])
            allra.append(ra)
            alldec.append(dec)
            allx.append(A.fplane.by_ifuslot(ifuslot).y + amppos[:, 0] +
                        dither_pattern[i, 0])
            ally.append(A.fplane.by_ifuslot(ifuslot).x + amppos[:, 1] +
                        dither_pattern[i, 1])

        t2 = time.time()
        cnt += 1
        time_per_amp = (t2 - t1) / cnt
        remaining_amps = (N - cnt)
        log.info('Time remaining: %0.2f' % (time_per_amp * remaining_amps))
    if breakloop:
        break

fitslist = [fits.PrimaryHDU(np.vstack(allspec)),
            fits.ImageHDU(np.array(allflatspec)),
            fits.ImageHDU(commonwave),
            fits.ImageHDU(np.array(allra)),
            fits.ImageHDU(np.array(alldec)),
            fits.ImageHDU(np.array(allx)),
            fits.ImageHDU(np.array(ally))]
fits.HDUList(fitslist).writeto('test_big.fits', overwrite=True)
flist1 = []
alls = np.vstack(allsub)
for j, resi in enumerate(alls):
    if j == 0:
        func = fits.PrimaryHDU
    else:
        func = fits.ImageHDU
    flist1.append(func(resi))
fits.HDUList(flist1).writeto('test_sub.fits', overwrite=True)
