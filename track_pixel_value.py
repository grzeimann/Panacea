# -*- coding: utf-8 -*-
"""
Created on Mon Jun 11 11:32:27 2018

@author: gregz
"""

import glob
import os.path as op
import numpy as np
import re

from amplifier import Amplifier
from distutils.dir_util import mkpath
from input_utils import setup_parser, set_daterange, setup_logging
from astropy.table import Table


def build_filenames(date, args):
    '''
    Build directory structure and search for unique observations, and return
    a single file for each observation to examine the header.
    '''
    basedir = op.join(args.rootdir, date, args.instrument,
                      args.instrument + '0000*', 'exp*', args.instrument)
    filenames_zro = sorted(glob.glob(op.join(basedir, '2*_zro.fits')))
    filenames_sci = sorted(glob.glob(op.join(basedir, '2*_sci.fits')))
    filenames = filenames_zro + filenames_sci
    dirnames = [op.dirname(fn) for fn in filenames]
    unique_dirnames, ind = np.unique(dirnames, return_index=True)
    return [(filenames[i][:-14], filenames[i][-8:-5]) for i in ind]


def Track_pixel_value(file_list, ifuslot, amp, args, date, yran=[10, 30],
                      xran=[900, 1200]):
    # Create empty lists for the left edge jump, right edge jump, and structure
    big_array = np.zeros(((xran[1] - xran[0]) * (yran[1] - yran[0]) + 2,
                          len(file_list)))
    amp_list = []
    for i, itm in enumerate(file_list):
        base, typ = itm
        fn = base + '%s%s_%s.fits' % (ifuslot, amp, typ)
        amp_list.append(Amplifier(fn, ''))
        amp_list[-1].subtract_overscan()
        amp_list[-1].trim_image()
        big_array[:, i+2] = amp_list[-1].image[yran[0]:(yran[1]+1),
                                               xran[0]:(xran[1]+1)].ravel()

    # Select only the bias frames that match the input amp, e.g., "RU"
    if not len(amp_list):
        args.log.warning('No zro or sci frames found for date range given')
        return None

    yind, xind = np.indices(amp_list[-1].image.shape)
    big_array[:, 0] = yind[yran[0]:(yran[1]+1), xran[0]:(xran[1]+1)].ravel()
    big_array[:, 1] = xind[yran[0]:(yran[1]+1), xran[0]:(xran[1]+1)].ravel()

    names = ['y', 'x']
    for ampn in amp_list:
        datetemp = re.split('[-,T]', ampn.header['DATE-OBS'])
        datev = datetemp[0] + datetemp[1] + datetemp[2]
        names.append('%s_%07d' % int(datev, ampn.header['OBSID']))
    T = Table(big_array, names=names)
    mkpath(op.join(args.folder, date))
    args.log.info('Writing pixels for %s, %s' % (amp_list[-1].specid, amp))
    T.write(op.join(args.folder, date, 'pixelvalues_%s_%s.dat' %
                    (amp_list[-1].specid, amp)),
            format='ascii.fixed_width_two_line')

parser = setup_parser()
parser.add_argument("-f", "--folder",
                    help='''Output folder''',
                    type=str, default='pixelvalues')

args = parser.parse_args(args=None)
args.log = setup_logging(logname='track_pixel_values')
args = set_daterange(args)
filenames = []
for date in args.daterange:
    date = '%04d%02d%02d' % (date.year, date.month, date.day)
    filenames = filenames + build_filenames(date, args)
for ifuslot in ['056']:  # ['056', '066']:
    for amp in ['RL']:  # ['LL', 'LU', 'RL', 'RU']:
        Track_pixel_value(filenames, ifuslot, amp, args,
                          '%04d%02d%02d' % (args.daterange[0].year,
                                            args.daterange[0].month,
                                            args.daterange[0].day))
