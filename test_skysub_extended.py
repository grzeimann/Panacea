# -*- coding: utf-8 -*-
"""
Created on Fri Oct  5 14:53:10 2018

@author: gregz
"""

import os.path as op
import sys
from astropy.io import fits
from astropy.table import Table
from utils import biweight_location
import numpy as np
from scipy.interpolate import LSQBivariateSpline, interp1d
from astropy.convolution import Gaussian1DKernel, interpolate_replace_nans
from scipy.signal import medfilt
import argparse as ap
from input_utils import setup_logging
import warnings

def get_script_path():
    return op.dirname(op.realpath(sys.argv[0]))

DIRNAME = get_script_path()

blueinfo = [['BL', 'uv', 'multi_503_056_7001', [3640., 4640.], ['LL', 'LU']],
            ['BR', 'orange', 'multi_503_056_7001', [4660., 6950.], ['RU', 'RL']]]
redinfo = [['RL', 'red', 'multi_502_066_7002', [6450., 8400.], ['LL', 'LU']],
           ['RR', 'farred', 'multi_502_066_7002', [8275., 10500.], ['RU', 'RL']]]
#list_of_blue = [['20170121', 'lrs20000010', 'exp01', '20170121', 'lrs20000011', 'exp01'],
#                ['20170126', 'lrs20000010', 'exp01', '20170121', 'lrs20000011', 'exp01'],
#                ['20170126', 'lrs20000010', 'exp02', '20170121', 'lrs20000011', 'exp01'],
#                ['20170126', 'lrs20000010', 'exp03', '20170121', 'lrs20000011', 'exp01'],
#                ['20170126', 'lrs20000010', 'exp04', '20170121', 'lrs20000011', 'exp01'],
#                ['20170126', 'lrs20000010', 'exp05', '20170121', 'lrs20000011', 'exp01']]
#list_of_red  = [['20170121', 'lrs20000011', 'exp01', '20170121', 'lrs20000010', 'exp01']]
parser = ap.ArgumentParser(add_help=True)

parser.add_argument("-b", "--basedir",
                    help='''base directory for reductions''',
                    type=str, default=None)
parser.add_argument("-s", "--side",
                    help='''blue for LRS2-B and red for LRS2-R''',
                    type=str, default='blue')
parser.add_argument("-scd", "--scidateobsexp",
                    help='''Example: "20180112,lrs20000027,exp01"''',
                    type=str, default=None)
parser.add_argument("-skd", "--skydateobsexp",
                    help='''Example: "20180112,lrs20000027,exp01"''',
                    type=str, default=None)
args = parser.parse_args(args=None)
args.log = setup_logging('test_skysub')
if args.scidateobsexp is None:
    args.log.error('--scidateobsexp/-scd was not set.')
    sys.exit(1)
if args.skydateobsexp is None:
    args.log.error('--skydateobsexp/-skd was not set.')
    sys.exit(1)
if args.side == 'blue':
    list_of_blue = [args.scidateobsexp.split(',') +
                    args.skydateobsexp.split(',')]
if args.side == 'red':
    list_of_red = [args.scidateobsexp.split(',') +
                   args.skydateobsexp.split(',')]

basedir = op.join(args.basedir, '%s/lrs2/%s/%s/lrs2/%s')

skyline_file = op.join(DIRNAME, 'lrs2_config/%s_skylines.dat')


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
        D = np.sqrt((xloc[:, np.newaxis, np.newaxis] - Dx[k] - xgrid)**2 +
                    (yloc[:, np.newaxis, np.newaxis] - Dy[k] - ygrid)**2)
        W = np.exp(-0.5 / (seeing/2.35)**2 * D**2)
        N = W.sum(axis=0)
        zgrid[k, :, :] = ((data[:, k][:, np.newaxis, np.newaxis] *
                           W).sum(axis=0) / N / scale**2 / area)
    wi = np.searchsorted(wave, wstart, side='left')
    we = np.searchsorted(wave, wend, side='right')

    zimage = biweight_location(zgrid[wi:we+1], axis=(0,))
    return zgrid, zimage, xgrid, ygrid


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


def get_info(basefile, amps, lims, in_wave=None):
    ifup, spectrum, wave = ([], [], [])
    for amp in amps:
        sfile = basefile + '_%s.fits' % amp
        sci = fits.open(sfile)
        ifup.append(sci['ifupos'].data)
        if in_wave is None:
            wt = np.array(sci['wavelength'].data, dtype=float)
        else:
            N, D = in_wave.shape
            if (amp == 'LL') or (amp == 'RU'):
                wt = in_wave[:N/2, :]
            else:
                wt = in_wave[N/2:, :]
        rw, rs = rectify(wt, np.array(sci['spectrum'].data, dtype=float),
                         lims, fac=1.5)
        spectrum.append(rs)
        wave.append(rw)
    ifup, spectrum, wave = [np.vstack(x) for x in [ifup, spectrum, wave]]
    return ifup, spectrum, wave


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


def create_image_header(wave,  xgrid, ygrid, zgrid, func=fits.ImageHDU):
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
    mad1 = np.median(y1)
    mad2 = np.median(y2)
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
        C[:, 1] = sky_list[1][:, j]
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            sol = np.linalg.lstsq(C[sel], sci_list[1][sel, j])[0]
        norm1[j] = sol[0]
        norm2[j] = sol[1]
        newsci[:, j] = C[:, 0] * sol[0]
        newsky[:, j] = C[:, 1] * sol[1]
    return newsky, newsci, norm1, norm2


def main(reduc_info, info_list):
    scidate, sciobs, sciexp, skydate, skyobs, skyexp = reduc_info
    print('Working on %s %s %s' % (scidate, sciobs, sciexp))
    for side in info_list:
        specinit, specname, multi, lims, amps = side
        W = fits.open('/Users/gregz/cure/panacea/lrs2_config/'
                      '%s_wavelength.fits' % specname)
        sky_file = basedir % (skydate, skyobs, skyexp, multi)
        sci_file = basedir % (scidate, sciobs, sciexp, multi)
        darfile = '/Users/gregz/cure/panacea/lrs2_config/dar_%s.dat' % specinit
        T = Table.read(darfile, format='ascii.fixed_width_two_line')
        sci_list = get_info(sci_file, amps, lims, in_wave=W[0].data)
        sky_list = get_info(sky_file, amps, lims, in_wave=W[0].data)
        x, y = (sci_list[0][:, 0], sci_list[0][:, 1])
        wave = sci_list[2][0]
        wave_0 = np.mean(wave)
        xoff = (np.interp(wave, T['wave'], T['x_0']) -
                np.interp(wave_0, T['wave'], T['x_0']))
        yoff = (np.interp(wave, T['wave'], T['y_0']) -
                np.interp(wave_0, T['wave'], T['y_0']))

        # Initial Models
        xn, yn = (0., 0.)
        sel = np.where(((x - xn)**2 + (y-yn)**2) > 5.0**2)[0]
        v = biweight_location(sci_list[1][sel, :] / sky_list[1][sel, :],
                              axis=(0,))
        gal_image = biweight_location(sci_list[1] - v * sky_list[1], axis=(1,))
        sky, temp, norm1, norm2 = solve_system(sci_list, sky_list, x, y, xoff,
                                               yoff, gal_image)
        skysub = sci_list[1] - sky
        for S, name in zip([sky, skysub], ['sky', 'skysub']):
            outname = '%s_%s_%s_%s_%s_cube.fits' % (scidate, sciobs, sciexp,
                                                    specname, name)
            zcube, zimage, xgrid, ygrid = make_frame(x, y, S, wave, T['wave'],
                                                     xoff, yoff,
                                                     wstart=wave_0-50.,
                                                     wend=wave_0+50.)
            write_cube(wave, xgrid, ygrid, zcube, outname)
        outname = '%s_%s_%s_%s_%s.fits' % ('multi', scidate, sciobs, sciexp,
                                           specname)
        X = np.array([T['wave'], T['x_0'], T['y_0']])
        f1 = create_header_objection(wave, sci_list[1], func=fits.PrimaryHDU)
        f2 = create_header_objection(wave, sky)
        f3 = create_header_objection(wave, skysub)
        fits.HDUList([f1, f2, f3, fits.ImageHDU(sci_list[0]),
                      fits.ImageHDU(wave), fits.ImageHDU(zimage),
                      fits.ImageHDU(X)]).writeto(outname, overwrite=True)

if args.side == 'blue':
    for blue in list_of_blue:
        main(blue, blueinfo)

if args.side == 'red':
    for red in list_of_red:
        main(red, redinfo)
