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
from astropy.convolution import convolve
from scipy.signal import medfilt, savgol_filter
from skimage.feature import register_translation
import argparse as ap
from input_utils import setup_logging
import warnings
from astropy.modeling.models import Polynomial2D
from astropy.modeling.fitting import LevMarLSQFitter


get_newwave = True

def get_script_path():
    return op.dirname(op.realpath(sys.argv[0]))

DIRNAME = get_script_path()

blueinfo = [['BL', 'uv', 'multi_503_056_7001', [3640., 4640.], ['LL', 'LU'],
             [4350., 4375.]], ['BR', 'orange', 'multi_503_056_7001',
            [4660., 6950.], ['RU', 'RL'], [6270., 6470.]]]
redinfo = [['RL', 'red', 'multi_502_066_7002', [6450., 8400.], ['LL', 'LU'],
            [7225., 7425.]], ['RR', 'farred', 'multi_502_066_7002',
           [8275., 10500.], ['RU', 'RL'], [9280., 9530.]]]

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

targs = ["-b", "/Users/gregz/cure/reductions",
         "-s", "red", "-scd", "20181108,lrs20000025,exp01", "-skd",
         "20181108,lrs20000024,exp01"]

args = parser.parse_args(args=targs)
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


def align_wave_with_sky(wave, sky, l1, l2, error):
    skyline = get_skyline_file(op.join(DIRNAME,
                                       'lrs2_config/airglow_groups.dat'))
    sel = np.where((skyline[:, 0] > l1) * (skyline[:, 0] < l2))[0]
    kshift = np.zeros((sky.shape[0],))
    p = kshift * 0.
    G = Gaussian1DKernel(1.5)
    for i in np.arange(sky.shape[0]):
        xl = np.searchsorted(wave[i], l1)
        xh = np.searchsorted(wave[i], l2)
        p[i] = np.percentile(sky[i, xl:xh] / error[i, xl:xh], 98)
    if np.median(p) < 100.:
        args.log.info('Low S/N regime for sky fitting.')
    P = np.median(p)
    nchunks = np.min([sky.shape[0]/4, int(sky.shape[0] / (100. / P))])
    args.log.info('Using %i chunks b/c individual fiber s/n is %0.2f' %
                  (nchunks, P))
    nshift, fshift = ([], [])
    for wchunk, schunk, fchunk in zip(np.array_split(wave, nchunks),
                                      np.array_split(sky, nchunks),
                                      np.array_split(np.arange(
                                                     sky.shape[0]),
                                                     nchunks)):
        nw, ns = make_avg_spec(wchunk, schunk, binsize=(sky.shape[0] /
                                                        nchunks))
        xl = np.searchsorted(nw, l1)
        xh = np.searchsorted(nw, l2)
        m1 = sky.shape[1] / 2
        dw_pix = np.mean(wchunk[:, m1+1] - wchunk[:, m1])
        skymodel = make_skyline_model(nw, skyline[sel, :], 1.,
                                      kernel_size=2.1, dw_pix=dw_pix)
        cont = fit_continuum(nw, ns)
        y = convolve(ns-cont, G)
        dw = nw[xh] - nw[xh-1]
        shift = register_translation(skymodel[xl:xh, np.newaxis],
                                     y[xl:xh, np.newaxis],
                                     upsample_factor=100)
        nshift.append(shift[0][0]*dw)
        fshift.append(np.mean(fchunk))
#        import matplotlib.pyplot as plt
#        plt.figure
#        norm = np.max(y[xl:xh]) / np.max(skymodel[xl:xh])
#        plt.plot(nw + shift[0][0]*dw, y)
#        plt.plot(nw, skymodel*norm)
#        plt.ylim([-5, np.nanpercentile((ns-cont)[xl:xh], 99)*1.8])
#        plt.xlim([nw[xl], nw[xh]])
#        plt.show()
#        ans = raw_input('%0.2f' % shift[0][0])
#        if ans == 'q':
#            sys.exit(1)
    args.log.info(nshift)
    y = np.array(nshift)
    absy = np.abs(y - medfilt(y, 5))
    mad = np.median(absy)
    sel = absy < 2. * mad
    p = np.polyval(np.polyfit(np.array(fshift)[sel], np.array(nshift)[sel], 3),
                   np.arange(sky.shape[0]))
    args.log.info('Average wavelength offset: %0.3f' % np.median(p))
    return wave + p[:, np.newaxis]


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
        spectrum.append(sci['spectrum'].data)
        wave.append(wt)
    ifup, spectrum, wave = [np.vstack(x)
                                 for x in [ifup, spectrum, wave]]
    return [ifup, spectrum, wave]


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


def main(reduc_info, info_list):
    scidate, sciobs, sciexp, skydate, skyobs, skyexp = reduc_info
    print('Working on %s %s %s' % (scidate, sciobs, sciexp))
    for side in info_list:
        specinit, specname, multi, lims, amps, slims = side
        if specname == 'uv':
            W1 = fits.open('/Users/gregz/cure/panacea/lrs2_config/'
                           '%s_%s_wavelength.fits' % (specname, 'LL'))
            W2 = fits.open('/Users/gregz/cure/panacea/lrs2_config/'
                           '%s_%s_wavelength.fits' % (specname, 'LU'))
            W = np.vstack([W1[0].data, W2[0].data])
        elif specname == 'orange':
            W1 = fits.open('/Users/gregz/cure/panacea/lrs2_config/'
                           '%s_%s_wavelength.fits' % (specname, 'RU'))
            W2 = fits.open('/Users/gregz/cure/panacea/lrs2_config/'
                           '%s_%s_wavelength.fits' % (specname, 'RL'))
            W = np.vstack([W1[0].data, W2[0].data])
        else:
            W = fits.open('/Users/gregz/cure/panacea/lrs2_config/'
                          '%s_wavelength.fits' % specname)[0].data
        sky_file = basedir % (skydate, skyobs, skyexp, multi)
        sci_file = basedir % (scidate, sciobs, sciexp, multi)
        darfile = '/Users/gregz/cure/panacea/lrs2_config/dar_%s.dat' % specinit
        T = Table.read(darfile, format='ascii.fixed_width_two_line')
        sci_list = get_info(sci_file, amps, lims, in_wave=None)
        sky_list = get_info(sky_file, amps, lims, in_wave=None)
        sky = np.where(sky_list[1] < 0., 0., sky_list[1])
        error_sky = np.sqrt(2.*3**2 + 0.8 * sky)
        if get_newwave:
            newwave = sky_list[2] * 0.
            args.log.info('Getting new wavelength solution from sky for amp 1.')
            newwave[:140] = align_wave_with_sky(sky_list[2][:140],
                                                sky_list[1][:140], slims[0],
                                                slims[1], error_sky[:140])
            args.log.info('Getting new wavelength solution from sky for amp 2.')
            newwave[140:] = align_wave_with_sky(sky_list[2][140:],
                                                sky_list[1][140:], slims[0],
                                                slims[1], error_sky[:140])
        else:
            newwave = sky_list[2]
        rw, rs = rectify(newwave, np.array(sci_list[1], dtype=float),
                         lims, fac=1.5)
        sci_list[2] = rw*1.
        sci_list[1] = rs*1.
        rw, rs = rectify(newwave, np.array(sky_list[1], dtype=float),
                         lims, fac=1.5)
        sky_list[2] = rw*1.
        sky_list[1] = rs*1.
        x, y = (sci_list[0][:, 0], sci_list[0][:, 1])
        wave = sci_list[2]
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
        loc = np.argmax(gal_image)
        args.log.info('Peak found at %0.2f, %0.2f' % (x[loc], y[loc]))
        xn, yn = (x[loc], y[loc])
        d = (x - xn)**2 + (y-yn)**2
        thresh = np.percentile(d, 90)
        sel = np.where(((x - xn)**2 + (y-yn)**2) > thresh)[0]
        from astropy.stats import sigma_clipped_stats
        v = biweight_location(sci_list[1][sel, :] / sky_list[1][sel, :],
                              axis=(0,))
        XN = biweight_location(sci_list[1], axis=(1,))
        YN = biweight_location(v*sky_list[1], axis=(1,))
        data = XN / YN
        mean, median, std = sigma_clipped_stats(info[0]/info[1])
        sel = np.abs(data - median) < 3. * std
        P = Polynomial2D(2)
        fitter = LevMarLSQFitter()
        fit = fitter(P, x[sel], y[sel], data[sel])
        offset = fit(x, y)
        sky_list[1] = sky_list[1] * offset[:, np.newaxis]
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
        info = main(blue, blueinfo)

if args.side == 'red':
    for red in list_of_red:
        info = main(red, redinfo)
