# -*- coding: utf-8 -*-
"""
Created on Wed May 23 17:39:02 2018

@author: gregz
"""

import glob
import numpy as np
import os.path as op
import sys
from astropy.io import fits
from input_utils import setup_logging, setup_basic_parser
from scipy.interpolate import splev, splrep
from utils import biweight_location


def write_to_fits(hdu, outname):
        '''
        Writing fits file to outname
        '''
        try:
            hdu.writeto(outname, overwrite=True)
        except TypeError:
            hdu.writeto(outname, clobber=True)


def get_folder(args):
    folder = op.join(args.rootdir, args.date,
                     args.instrument, args.instrument +
                     "%07d" % (int(args.observation)),
                     'exp%02d' % (int(args.exposure_number)),
                     args.instrument)
    return folder


def build_cubename(args, prefix, basename):
    path = get_folder(args)
    return op.join(path, '%s%s_%s_sci.fits' %
                   (prefix, basename, args.ifuslot))


def get_spectrum(args, amp):
    folder = get_folder(args)
    filenames = sorted(glob.glob(op.join(folder, 'multi_*_%s_*_%s.fits' %
                                         (args.ifuslot, amp))))
    if len(filenames):
        F = fits.open(filenames[0])
    else:
        args.log.error('File does not exist: %s' % folder)
        sys.exit(1)
    try:
        wv, spec, ftf, x, y = (F['wavelength'].data, F[args.spectype].data,
                               F['fiber_to_fiber'].data,
                               F['ifupos'].data[:, 0], F['ifupos'].data[:, 1])
        basename = op.basename(F[0].header['rawfn']).split('_')[0]
    except:
        args.log.error('Could not get extensions: wavelength and %s' %
                       args.spectype)
        sys.exit(1)
    return wv, spec, ftf, x, y, basename


def rectify(wave, spec, rectified_dlam=1., minwave=None, maxwave=None):
    ''' Rectify spectra to same "rect_wave" '''
    dlam = np.zeros(wave.shape)
    dlam[:, 1:] = np.diff(wave, axis=1)
    dlam[:, 0] = dlam[:, 1]
    if rectified_dlam is None:
        rectified_dlam = np.nanmedian(dlam)

    rect_wave = np.arange(wave.min(), wave.max() + rectified_dlam,
                          rectified_dlam)
    if minwave is not None and maxwave is not None:
        wnew = np.arange(minwave, maxwave + rectified_dlam,
                         rectified_dlam)
    else:
        wnew = rect_wave * 1.
    rect_spec = np.zeros((spec.shape[0], len(wnew)))
    xs = np.linspace(0, 1, len(rect_wave))
    xn = np.interp(wnew, rect_wave, xs)
    for i in np.arange(spec.shape[0]):
        if np.all(spec[i] == 0):
            rect_spec[i, :] = 0.0
        else:
            y = spec[i] / dlam[i]
            xp = np.interp(wave[i], rect_wave, xs)
            tck = splrep(xp, y)
            rect_spec[i, :] = splev(xn, tck)
    rect_wave = wnew * 1.
    return rect_wave, rect_spec


def make_frame(xloc, yloc, data, wave, args, outname, ftf, scale=1.,
               seeing_fac=1.5):
    if args.test:
        data = 0. * data
        data[args.test_fiber] = 1.
        args.log.info('X, Y: %0.2f, %0.2f' % (xloc[args.test_fiber],
                                              yloc[args.test_fiber]))
    seeing = seeing_fac * scale
    a, b = data.shape
    x = np.arange(xloc.min()-scale,
                  xloc.max()+1*scale, scale)
    y = np.arange(yloc.min()-scale,
                  yloc.max()+1*scale, scale)
    xgrid, ygrid = np.meshgrid(x, y)
    zgrid = np.zeros((b,)+xgrid.shape)
    d = np.zeros((a,)+xgrid.shape)
    w = np.zeros((a,)+xgrid.shape)
    area = (1.5 / 2.)**2 * np.pi
    for i in np.arange(len(x)):
        for j in np.arange(len(y)):
            d[:, j, i] = np.sqrt((xloc - xgrid[j, i])**2 +
                                 (yloc - ygrid[j, i])**2)
            w[:, j, i] = np.exp(-1./2.*(d[:, j, i]/seeing)**2)

    avg = np.median(ftf, axis=1)
    sel = np.where((avg > 0.1) * np.isfinite(avg))[0]
    ws = w[sel, :, :].sum(axis=0)
    for k in xrange(b):
        zgrid[k, :, :] = ((data[sel, k][:, np.newaxis, np.newaxis] *
                           w[sel]).sum(axis=0) / ws * 1.9)# * scale**2 / area)
    wi = np.searchsorted(wave, args.wavestart, side='left')
    we = np.searchsorted(wave, args.waveend, side='right')

    zimage = biweight_location(zgrid[wi:we+1], axis=(0,))
    hdu = fits.PrimaryHDU(np.array(zimage, dtype='float32'))
    hdu.header['CRVAL1'] = x[0]
    hdu.header['CRVAL2'] = y[0]
    hdu.header['CRPIX1'] = 1
    hdu.header['CRPIX2'] = 1
    hdu.header['CTYPE1'] = 'pixel'
    hdu.header['CTYPE2'] = 'pixel'
    hdu.header['CDELT1'] = scale
    hdu.header['CDELT2'] = scale
    write_to_fits(hdu, outname)
    return zgrid, zimage, xgrid, ygrid


def main():
    parser = setup_basic_parser()
    parser.add_argument("-st", "--spectype",
                        help='''spectrum or sky_subtracted''',
                        type=str, default='sky_subtracted')
    parser.add_argument("-ws", "--wavestart",
                        help='''Start wavelength for collapse''',
                        type=float, default=4900)
    parser.add_argument("-we", "--waveend",
                        help='''End Wavelength for collapse''',
                        type=float, default=5350)
    parser.add_argument("-t", "--test",
                        help='''If test, then set test_fiber as well''',
                        action="count", default=0)
    parser.add_argument("-tf", "--test_fiber",
                        help='''fiber for testing centering''',
                        type=int, default=0)
    args = parser.parse_args(args=None)
    args.log = setup_logging(logname='collapsed')
    allwv, allspec, allx, ally, allftf = ([], [], [], [], [])
    for amp in ['LL', 'LU', 'RL', 'RU']:
        wv, spec, ftf, x, y, basename = get_spectrum(args, amp)
        allwv.append(wv)
        allspec.append(spec)
        allftf.append(ftf)
        allx.append(x)
        ally.append(y)
    allwv, allspec, allftf  = [np.vstack(i) for i in [allwv, allspec, allftf]]
    allx, ally = [np.hstack(i) for i in [allx, ally]]
    wave, spec = rectify(allwv, allspec, minwave=3500, maxwave=5500)
    if args.test:
        base = 'test_CoFeS'
    else:
        base = 'CoFeS'
    outname = build_cubename(args, base, basename)
    make_frame(allx, ally, spec, wave, args, outname, allftf)

if __name__ == '__main__':
    main()
