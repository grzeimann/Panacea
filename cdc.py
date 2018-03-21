# -*- coding: utf-8 -*-
"""
Created on Wed Mar 21 13:45:39 2018

@author: gregz
"""
from input_utils import setup_basic_parser, setup_logging
import os.path as op
import numpy as np
import glob


def build_filenames(args):
    '''
    Build directory structure and search for unique observations, and return
    a single file for each observation to examine the header.
    '''
    if args.exposure_number is None:
        expstr = 'exp*'
    else:
        expstr = 'exp%02d' % int(args.exposure_number)

    filename = op.join(args.rootdir, args.date, args.instrument,
                       '%s%07d' % (args.instrument, int(args.observation)),
                       expstr, args.instrument, 'multi*LL.fits')
    print(filename)
    filenames = glob.glob(filename)
    ifuslot_list = [op.basename(fn).split('_')[2] for fn in filenames]
    ifuslots = np.unique(ifuslot_list)
    exposure_list = [op.basename(op.dirname(op.dirname(fn)))[3:]
                     for fn in filenames]
    exposures = np.unique(exposure_list)

    return filenames, ifuslots, exposures, ifuslot_list, exposure_list

parser = setup_basic_parser()
args = parser.parse_args(args=None)
args.log = setup_logging(logname='amazeballs')
filenames, ifuslots, exposures, i_list, e_list = build_filenames(args)
args.log.info(exposures)
for exposure in exposures:
    file_list = [fn for fn, e in zip(filenames, e_list) if e == exposure]
    ifuslot_list = [i for i, e in zip(i_list, e_list) if e == exposure]
    args.log.info(ifuslot_list)

    
