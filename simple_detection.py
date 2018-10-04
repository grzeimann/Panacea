# -*- coding: utf-8 -*-
"""
Simple Detection Code
"""

import matplotlib
matplotlib.use('agg')
import glob
import sys
import numpy as np
import os.path as op
import fitsio
import argparse as ap

from astropy.io import fits
from astropy.convolution import convolve, Gaussian1DKernel, interpolate_replace_nans
from input_utils import setup_logging
from scipy.interpolate import interp1d
from astropy.stats import mad_std, biweight_location
from utils import biweight_midvariance
from astropy.table import Table
from scipy.signal import savgol_filter, medfilt


def get_script_path():
    return op.dirname(op.realpath(sys.argv[0]))

DIRNAME = get_script_path()


def grab_attribute(filename, args, attributes=[], amps=['LL', 'LU', 'RU',
                   'RL']):
    ''' grab specified attributes from multi* file '''
    basename = filename[:-8]
    s = [[] for a in attributes]
    norm = []
    for amp in amps:
        name = basename + '_%s.fits' % amp
        data = fitsio.read(name, 'twi_spectrum')
        xl = data.shape[1] / 3
        xh = 2 * data.shape[1] / 3
        norm.append(np.median(data[:, xl:xh], axis=1))
    for amp in amps:
        name = basename + '_%s.fits' % amp
        for i, attribute in enumerate(attributes):
            s[i].append(fitsio.read(name, attribute))
    X = [np.vstack(si) for si in s]
    X.append(np.hstack(norm))
    return X


def put_attribute(filename, args, data, attributes=[]):
    ''' put specified attributes into multi* file '''
    try:
        for i, attribute in enumerate(attributes):
            F = fitsio.FITS(filename, 'rw')
            F.write(data[i], extname=attribute+'_1')
    except IOError:
        for i, attribute in enumerate(attributes):
            args.log.warning('%s not found to add %s' % attribute)


def rectify(wave, spec, lims, mask=None, fac=1.0):
    N, D = wave.shape
    rect_wave = np.linspace(lims[0], lims[1], int(D*fac))
    rect_spec = np.zeros((N, len(rect_wave)))
    G = Gaussian1DKernel(1.5 * fac)
    for i in np.arange(N):
        dw = np.diff(wave[i])
        dw = np.hstack([dw[0], dw])
        if mask is None:
            x = wave[i]
            y = spec[i] / dw
        else:
            x = wave[i]
            y = (spec[i] / dw)
            y[mask[i]] = np.nan
            y = interpolate_replace_nans(y, G)
        I = interp1d(x, y, kind='quadratic',
                     bounds_error=False, fill_value=-999.)
        rect_spec[i, :] = I(rect_wave)
    return rect_wave, rect_spec


def build_weight_matrix(x, y, sig=1.5):
    d = np.sqrt((x - x[:, np.newaxis])**2 + (y - y[:, np.newaxis])**2)
    G = np.exp(-0.5 * (d / sig)**2)
    G = G / G.sum(axis=0)[:, np.newaxis]
    return G.swapaxes(0, 1)


def convolve_spatially(x, y, spec, wave, mask, sig_spatial=0.7, sig_wave=1.5):
    W = build_weight_matrix(x, y, sig=sig_spatial)
    Z = spec * 1.
    G = Gaussian1DKernel(sig_wave)
    for i in np.arange(spec.shape[0]):
        Z[i, :] = convolve(Z[i, :], G, nan_treatment='fill', fill_value=0.0)
    for i in np.arange(spec.shape[1]):
        Z[:, i] = np.dot(Z[:, i], W)
    return Z


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


def mask_cosmics(error, trace, cosmic_avoidance=4):
    cyind, cxind = np.where(error == -1)
    mask = np.array(trace * 0., dtype=bool)
    for xind, yind in zip(cxind, cyind):
        trace_a = trace[:, xind]
        fibs = np.where(np.abs(trace_a - yind) < cosmic_avoidance)[0]
        for fib in fibs:
            lx = (xind-cosmic_avoidance)
            hx = (xind+cosmic_avoidance) + 1
            mask[fib, lx:hx] = True
    return mask


def dummy_test(image):
    y = savgol_filter(image, 315, 1, axis=1)
    s = np.zeros((y.shape[1], 4))
    for i in np.arange(y.shape[1]):
        chunks = np.array_split(y[:, i], 4)
        avg = chunks[-1]
        n = [biweight_location(avg / chunk) for chunk in chunks]
        s[i, :] = np.array(n)

    chunks = np.array_split(y, 4, axis=0)
    norm = np.zeros((3*448, chunks[0].shape[1]))
    for k in np.arange(3):
        for i in np.arange(y.shape[1]):
            x = np.arange(chunks[0].shape[0])
            z = (chunks[k] * s[:, k] / chunks[-1])[:, i]
            test = z - medfilt(z, 51)
            threshold = 3. * np.nanmedian(np.abs(test))
            mask = np.abs(test) < threshold
            xchunk, zchunk, mchunk = [np.array_split(j, 4, axis=0) for j in [x, z, mask]]
            for xc, zc, mc in zip(xchunk, zchunk, mchunk):
                p = np.polyfit((xc/448.)[mc], zc[mc], 2)
                norm[xc+448*k, i] = np.polyval(p, xc/448.)
    return norm

def main():
    parser = ap.ArgumentParser(add_help=True)

    parser.add_argument("-f", "--filename",
                        help='''Filename that contains list of files''',
                        type=str, default=None)
    parser.add_argument("-ac", "--spatial_conv_size",
                        help='''Spatial Convolution Kernel Sigma (")''',
                        type=float, default=0.6)
    parser.add_argument("-ec", "--spectral_conv_size",
                        help='''Spectral Convolution Kernel Sigma (A)''',
                        type=float, default=2.5)
    parser.add_argument("-cc", "--spectral_cont_conv_size",
                        help='''Spectral Continuum Convolution Kernel Sigma (A)''',
                        type=float, default=25.)
    parser.add_argument("-t", "--threshold",
                        help='''Detection Threshold''',
                        type=float, default=5.0)
    parser.add_argument("-o", "--outdir",
                        help='''Out directory for detections''',
                        type=str, default='temp')
    args = parser.parse_args(args=None)
    args.log = setup_logging(logname='detection')

    filenames = [line.rstrip('\n').split()
                 for line in open(args.filename, 'r')]
    allwave, allspec, allifupos, allmask, alltwi, allmodel, allftf = ([], [], [], [], [], [], [])
    for filename in filenames:
        args.log.info('Reading in %s' % filename[0][:-8])
        dither = np.array([float(filename[2]), float(filename[3])])
        amps = ['LL', 'LU', 'RU', 'RL']
        attributes = ['wavelength', 'spectrum', 'twi_spectrum',
                      'ifupos', 'error', 'trace', 0, 'flat_image']
        w, s, f, i, e, t, T, m, n  = grab_attribute(filename[0], args,
                                             attributes=attributes, amps=amps)
        mask = mask_cosmics(e, t)
        norm = (n / np.median(n))[:, np.newaxis]
        allwave.append(w)
        allspec.append(s)#safe_division(s, f * norm))
        allifupos.append(i + dither)
        allmask.append(mask)
        alltwi.append(T)
        allmodel.append(m)
        allftf.append(f)
    allwave, allspec, allifupos, allmask, alltwi, allmodel, allftf = [np.array(np.vstack(x), dtype='float64')
                                            for x in [allwave, allspec,
                                                      allifupos, allmask, alltwi, allmodel, allftf]]
    args.log.info('Rectifying sky subtracted spectra')
    allmask = np.array(allmask, dtype=bool)
    rw, rs = rectify(allwave, allspec, [3500., 5500.], mask=allmask, fac=1.5)
    args.log.info('Convolving sky subtracted spectra for continuum')
    Zc = convolve_spatially(allifupos[:, 0], allifupos[:, 1], rs, rw, allmask,
                           sig_spatial=args.spatial_conv_size,
                           sig_wave=(args.spectral_cont_conv_size / (rw[1]-rw[0])))
    noise = biweight_midvariance(Zc, axis=(0, ))
    SNc = Zc / noise
    F1 = fits.PrimaryHDU(alltwi)
    args.log.info('Convolving sky subtracted spectra for emission')
    Ze = convolve_spatially(allifupos[:, 0], allifupos[:, 1], rs, rw, allmask,
                           sig_spatial=args.spatial_conv_size,
                           sig_wave=(args.spectral_conv_size / (rw[1]-rw[0])))
    noise = biweight_midvariance(Ze-Zc, axis=(0, ))
    SNe = (Ze-Zc) / noise
    F2 = fits.ImageHDU(allmodel)
    #norm = dummy_test(np.vstack([allspec, alltwi[:448,:]]))
    F3 = fits.ImageHDU(allftf)
    fits.HDUList([F1, F2, F3]).writeto('test.fits', overwrite=True)
    # peaks_fib, peaks_wave = np.where(SN > args.threshold)              
    
if __name__ == '__main__':
    main()
