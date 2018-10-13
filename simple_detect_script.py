# -*- coding: utf-8 -*-
"""
Created on Thu Oct  4 13:30:06 2018

@author: gregz
"""
from astropy.io import fits
from utils import biweight_location
import numpy as np
from scipy.signal import savgol_filter
import os.path as op
import glob
from scipy.interpolate import interp1d
from input_utils import setup_logging
import warnings
from astrometry import Astrometry

dither_pattern = np.array([[0., 0.], [1.27, -0.73], [1.27, 0.73]])
virus_amps = ['LL', 'LU', 'RU', 'RL']
lrs2_amps = [['LL', 'LU'], ['RL', 'RU']]
fplane_file = '/work/03730/gregz/maverick/fplane.txt'
flt_obs = '%07d' % 16
twi_obs = '%07d' % 02
sci_obs = '%07d' % 18
twi_date = '20181008'
sci_date = twi_date
flt_date = twi_date

log = setup_logging('panacea_quicklook')

basered = '/work/03730/gregz/maverick'
baseraw = '/work/03946/hetdex/maverick'

twi_path = op.join(basered, 'reductions', twi_date, '%s', '%s%s', 'exp01',
                   '%s', 'multi_*_%s_*_LL.fits')
sci_path = op.join(baseraw, sci_date,  '%s', '%s%s', 'exp%s',
                   '%s', '2*_%sLL*.fits')
flt_path = op.join(baseraw, flt_date,  '%s', '%s%s', 'exp*',
                   '%s', '2*_%sLL*.fits')


def get_cal_info(twi_path, amp):
    F = fits.open(glob.glob(twi_path.replace('LL', amp))[0])
    return F['ifupos'].data*1., F['trace'].data*1., F['wavelength'].data*1.


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


def make_avg_spec(wave, spec, binsize=35):
    ind = np.argsort(wave.ravel())
    T = 1
    for p in wave.shape:
        T *= p
    wchunks = np.array_split(wave.ravel()[ind],
                             T / binsize)
    schunks = np.array_split(spec.ravel()[ind],
                             T / binsize)
    nwave = np.array([np.mean(chunk) for chunk in wchunks])
    nspec = np.array([biweight_location(chunk) for chunk in schunks])
    nwave, nind = np.unique(nwave, return_index=True)
    return nwave, nspec[nind]


def base_reduction(filename):
    a = fits.open(filename)
    image = a[0].data
    # overscan sub
    overscan_length = 32 * (image.shape[1] / 1064)
    if image.shape[1] == 1064:
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


def get_flat_field(flt_path, amp, array_wave, array_trace, common_wave):
    files = glob.glob(flt_path.replace('LL', amp))
    listflt = []
    for filename in files:
        a = base_reduction(filename)
        listflt.append(a)
    array_flt = np.array(listflt)
    norm = np.median(array_flt, axis=(1, 2))
    array_flt = np.median(array_flt / norm[:, np.newaxis, np.newaxis], axis=0)
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
                           binsize=41)
    I = interp1d(nw, ns, kind='linear', fill_value='extrapolate')
    ftf = spectrum * 0.
    for fiber in np.arange(array_wave.shape[0]):
        model = I(array_wave[fiber])
        ftf[fiber] = savgol_filter(spectrum[fiber] / model, 151, 1)
    nw1, ns1 = make_avg_spec(array_wave, spectrum / ftf, binsize=41)
    Y, X = np.indices(array_wave.shape)
    YY, XX = np.indices(array_flt.shape)
    bigW = array_flt * 0.
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
    I = interp1d(nw1, ns1, kind='linear', fill_value='extrapolate')
    modelimage = I(bigW)
    flat = array_flt / modelimage
    return flat, bigW, I(common_wave)


# GET ALL VIRUS IFUSLOTS
twilist = glob.glob(twi_path % ('virus', 'virus', twi_obs, 'virus', '*'))
ifuslots = [op.basename(x).split('_')[2] for x in twilist]
fiberpos, fiberspec = ([], [])
log.info('Beginning the long haul.')
nexp = len(glob.glob(sci_path % ('virus', 'virus', flt_obs, '*', 'virus',
                                 ifuslots[0])))
header = fits.open(glob.glob(sci_path % ('virus', 'virus', sci_obs, '01',
                                         'virus', ifuslots[0]))[0])[0].header
PA = float(header['PARANGLE'])
RA = float(header['TRAJRA'])
DEC = float(header['TRAJDEC'])
log.info('Observation at %0.4f %0.4f, PA: %0.3f' % (RA, DEC, PA))
A = Astrometry(RA, DEC, PA, 0., 0., fplane_file=fplane_file)
allflats, allspec, allra, alldec = ([], [], [], [])

# Rectified wavelength
commonwave = np.linspace(3500, 5500, 1000)

for ifuslot in ifuslots:
    for amp in virus_amps:
        log.info('Starting on ifuslot, %s, and amp, %s' % (ifuslot, amp))
        twibase = twi_path % ('virus', 'virus', twi_obs, 'virus', ifuslot)
        amppos, trace, wave = get_cal_info(twibase, amp)
        fltbase = flt_path % ('virus', 'virus', flt_obs, 'virus', ifuslot)
        log.info('Getting Flat for ifuslot, %s, and amp, %s' % (ifuslot, amp))
        flat, bigW, flatspec = get_flat_field(fltbase, amp, wave, trace,
                                              commonwave)
        allflats.append(flatspec)
        for i in np.arange(nexp):
            log.info('Getting spectra for exposure, %i,  ifuslot, %s, and amp,'
                     ' %s' % (i+1, ifuslot, amp))
            ra, dec = A.get_ifupos_ra_dec(ifuslot,
                                          amppos[:, 0] + dither_pattern[i, 0],
                                          amppos[:, 1] + dither_pattern[i, 1])
            allra.append(ra)
            alldec.append(dec)
            scifile = glob.glob(sci_path % ('virus', 'virus', sci_obs,
                                            '%02d' % (i+1), 'virus',
                                            ifuslots[0]))[0].replace('LL', amp)
            image = base_reduction(scifile)
            spectrum = np.zeros(trace.shape[0], len(commonwave))
            temp = np.zeros((trace.shape[1], 6))
            temp2 = np.zeros((trace.shape[1], 6))
            x = np.arange(trace.shape[1])
            for fiber in np.arange(trace.shape[0]):
                indl = np.floor(trace[fiber]).astype(int)
                for k, j in enumerate(np.arange(-2, 4)):
                    temp[:, k] = image[indl, x]
                    temp2[:, k] = flat[indl, x]
                tempspec = (np.sum(temp * temp2, axis=1) /
                            np.sum(temp2, axis=1))
                I = interp1d(wave[fiber], spectrum[fiber], kind='quadratic',
                             fill_value='extrapolate')
                spectrum[fiber] = I(commonwave)
            allspec.append(spectrum)
fitslist = [fits.PrimaryHDU(np.array(allflats)),
            fits.ImageHDU(np.array(allspec)),
            fits.ImageHDU(commonwave),
            fits.ImageHDU(np.array(allra)),
            fits.ImageHDU(np.array(alldec))]
fits.HDUList(fitslist).writeto('test_big.fits', overwrite=True)
