# -*- coding: utf-8 -*-
"""
Build new Fiber to Fiber
"""

import matplotlib
matplotlib.use('agg')
import glob
import numpy as np
import os.path as op
import splinelab
import fitsio


from astropy.io import fits
from distutils.dir_util import mkpath
from input_utils import setup_parser, set_daterange, setup_logging
from scipy.interpolate import splev, splrep
from utils import biweight_location
from bspline import Bspline


def bspline_matrix(x, nknots, norm=False):
    ''' Make the bspline knot matrix for linalg calculation later '''
    v = np.linspace(0, 1, nknots)
    k = splinelab.augknt(v, 3)
    B = Bspline(k, 3)
    if norm:
        x = (x - x.min()) / (x.max() - x.min() + 0.1)
    c = np.array([B(xp) for xp in x])
    return B, c


def check_if_type(date, obsid, args):
    ''' Test if header has IMAGETYP '''
    filenames = glob.glob(op.join(args.rootdir, date, args.instrument,
                                  args.instrument + obsid, 'exp01',
                                  args.instrument, 'multi_*_*_*_LL.fits'))
    try:
        kind = fits.open(filenames[0])[0].header['IMAGETYP']
    except:
        args.log.warn('No IMAGETYP in header for %s and observation %s'
                      % (date, obsid))
        return False
    if kind == args.type:
        return True
    else:
        return False


def build_filenames(date, obsid, args):
    '''
    Build directory structure and search for all the files in a given
    observation and exposure.
    '''
    if args.type == 'twi':
        expstr = '01'
    else:
        expstr = '*'
    filenames = glob.glob(op.join(args.rootdir, date, args.instrument,
                                  args.instrument + obsid, 'exp%s' % expstr,
                                  args.instrument, 'multi_*_*_*_LL.fits'))
    ifuslot_list = [op.basename(fn).split('_')[2] for fn in filenames]
    ifuslots = np.unique(ifuslot_list)
    exposure_list = [op.basename(op.dirname(op.dirname(fn)))[3:]
                     for fn in filenames]
    exposures = np.unique(exposure_list)
    return filenames, ifuslots, exposures, ifuslot_list, exposure_list


def grab_attribute(filename, args, attributes=[],
                   amps=['LL', 'LU', 'RU', 'RL']):
    ''' grab specified attributes from multi* file '''
    basename = filename[:-8]
    s = [[] for a in attributes]
    for amp in amps:
        name = basename + '_%s.fits' % amp
        try:
            for i, attribute in enumerate(attributes):
                s[i].append(fitsio.read(name, attribute))
        except IOError:
            args.log.warning('%s not found, filling with zeros' % name)
            for i, attribute in enumerate(attributes):
                s[i].append(np.zeros((112, 1032)))
        for i, attribute in enumerate(attributes):
            if s[i][-1].shape != (112, 1032):
                s[i][-1] = np.zeros((112, 1032))

    return [np.array(si) for si in s]


def put_attribute(filename, args, data, attributes=[]):
    ''' put specified attributes into multi* file '''
    try:
        for i, attribute in enumerate(attributes):
            F = fitsio.FITS(filename, 'rw')
            F.write(data[i], extname=attribute+'_1')
    except IOError:
        for i, attribute in enumerate(attributes):
            args.log.warning('%s not found to add %s' % attribute)


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


def main():
    parser = setup_parser()
    parser.add_argument("-t", "--type",
                        help='''Observation Type, twi or sci''',
                        type=str, default='twi')
    parser.add_argument("-n", "--nknots",
                        help='''Number of knots for bspline''',
                        type=int, default=7)
    parser.add_argument("-b", "--nbins",
                        help='''Number of bins to collapse data
                        for bspline fit''', type=int, default=40)
    args = parser.parse_args(args=None)
    args.log = setup_logging(logname='build_ftf')
    args = set_daterange(args)

    ifu_spline_dict = {}
    filename_dict = {}
    # HARDCODED SIZE FOR SPEED BUT MUST MATCH SIZE OF "rw" BELOW.
    B, C = bspline_matrix(np.linspace(0, 1, 2001), nknots=args.nknots)
    for datet in args.daterange:
        date = '%04d%02d%02d' % (datet.year, datet.month, datet.day)
        obsids = glob.glob(op.join(args.rootdir, date, args.instrument,
                                   args.instrument + '*'))
        obsids = [obsid[-7:] for obsid in obsids]
        for obsid in obsids:
            if not check_if_type(date, obsid, args):
                continue
            filenames, ifus, exps, i_list, e_list = build_filenames(date,
                                                                    obsid,
                                                                    args)
            for exposure in exps:
                file_list = [fn for fn, e in zip(filenames, e_list)
                             if e == exposure]
                ifuslot_list = [i for i, e in zip(i_list, e_list)
                                if e == exposure]
                ifuslot_amp = ['%s%s' % (ifu, amp) for ifu in ifuslot_list
                               for amp in ['LL', 'LU', 'RU', 'RL']]
                for ifua in ifuslot_amp:
                    if ifua not in ifu_spline_dict:
                        ifu_spline_dict[ifua] = []
                        filename_dict[ifua] = []
                args.log.info('Building Fiber to Fiber for %s, observation %s,'
                              ' exposure %s' % (date, obsid, exposure))
                allspec = []
                for filen, ifu in zip(file_list, ifuslot_list):
                    args.log.info('Reading in %s' % filen)
                    amps = ['LL', 'LU', 'RU', 'RL']
                    wave, spec = grab_attribute(filen, args,
                                                attributes=['wavelength',
                                                            'spectrum'],
                                                amps=amps)
                    for wv, sp, amp in zip(wave, spec, amps):
                        rw, rs = rectify(wv, sp, minwave=3500.,
                                         maxwave=5500.)
                        allspec.append(rs)
                        name = filen[:-8] + '_%s.fits' % amp
                        filename_dict['%s%s' % (ifu, amp)].append(name)
                allspec = np.array(allspec)
                avgspec = np.nanmedian(allspec, axis=(0, 1))
                X = np.arange(len(rw))
                XL = np.array_split(X, args.nbins)
                xloc = np.array([np.median(xl) for xl in XL])
                xloc = (xloc - 0.) / (len(rw) - 1.)
                B, c = bspline_matrix(xloc, nknots=args.nknots)
                for sp, ifua in zip(allspec, ifuslot_amp):
                    args.log.info('Working on ifuslot %s' % ifua)
                    div = sp / avgspec
                    splinecoeff = np.zeros((sp.shape[0], c.shape[1]))
                    div_list = np.array_split(div, args.nbins, axis=1)
                    mdiv = [np.nanmedian(d, axis=1) for d in div_list]
                    mdiv = np.array(mdiv).swapaxes(0, 1)
                    for i, fiber in enumerate(mdiv):
                        sel = np.where(np.isfinite(fiber))[0]
                        splinecoeff[i, :] = np.linalg.lstsq(c[sel, :],
                                                            fiber[sel])[0]
                    ifu_spline_dict[ifua].append(splinecoeff)

    ifu_ftf_dict = {}
    for ifu in ifu_spline_dict:
        fibers_spline_coeff = biweight_location(np.array(ifu_spline_dict[ifu]),
                                                axis=(0,))
        ftf = np.zeros(rs.shape)
        for i, fiber in enumerate(fibers_spline_coeff):
            ftf[i, :] = np.dot(C, fiber)
        ifu_ftf_dict[ifu] = ftf
        for filename in filename_dict[ifu]:
            args.log.info('Writing Fiber to Fiber to %s' % filename)
            put_attribute(filename, args, [ftf], attributes=['fiber_to_fiber'])

if __name__ == '__main__':
    main()
