# -*- coding: utf-8 -*-
"""
Generate panacea command
"""

import glob
import numpy as np
import os.path as op
import sys

from astropy.io import fits
from datetime import datetime as dt
from input_utils import setup_parser, set_daterange, setup_logging

standard_names = ['HD_19445', 'SA95-42', 'GD50', 'HZ_4', 'G191B2B', 'FEIGE_25',
                  'HILTNER_600', 'G193-74', 'PG0823+546', 'HD_84937',
                  'GD108', 'FEIGE_34', 'HD93521', 'GD140', 'HZ_21',
                  'FEIGE_66', 'FEIGE_67', 'G60-54', 'HZ_44', 'GRW+70_5824',
                  'BD+26+2606', 'BD+33_2642', 'G138-31', 'WOLF_1346',
                  'BD+17_4708', 'FEIGE_110', 'GD248']


def build_filenames(date, args):
    '''
    Build directory structure and search for unique observations, and return
    a single file for each observation to examine the header.
    '''
    filenames = glob.glob(op.join(args.rootdir, date, args.instrument,
                                  args.instrument + '*', 'exp01',
                                  args.instrument, '*.fits'))
    dirnames = [op.dirname(fn) for fn in filenames]
    unique_dirnames, ind = np.unique(dirnames, return_index=True)
    return [filenames[i] for i in ind]


def find_match(datet, list_name):
    timediff = []
    for twi in list_name:
        date_twi = twi.split('_')[0]
        date_twi = dt(int(date_twi[:4]), int(date_twi[4:6]), int(date_twi[6:]))
        timediff.append((datet - date_twi).days)
    closest_date = np.argmin(np.abs(timediff))
    diff = np.min(np.abs(timediff))
    return closest_date, diff

parser = setup_parser()
parser.add_argument("-t", "--target",
                    help='''Target Name''',
                    type=str, default=None)
parser.add_argument("-s", "--side",
                    help='''red or blue''',
                    type=str, default='red')
args = parser.parse_args(args=['-sd', '20170422', '-dl', '1', '-r',
                               '/Users/gregz/cure/lrs2_raw', '-in', 'lrs2',
                               '-t', 'GC', '-s', 'blue'])
                    
if args.side.lower() == 'blue':
    ifuslot = '056'
else:
    ifuslot = '066'
args.log = setup_logging(logname='build_panacea_command')
if args.target is None:
    args.log.error('Please provide "--target" argument on call')
    sys.exit(1)
args = set_daterange(args)
panacea_dict = {}
science_target_list = []
twi_list = []
standard_list = []
for datet in args.daterange:
    date = '%04d%02d%02d' % (datet.year, datet.month, datet.day)
    filenames = build_filenames(date, args)
    for filename in filenames:
        obsid = op.basename(op.dirname(op.dirname(op.dirname(filename)))).split(args.instrument)[1]
        keystring = date+'_'+obsid
        objectname = fits.open(filename)[0].header['OBJECT']
        if args.target.lower() in objectname.lower():
            science_target_list.append(keystring)
            print(keystring, objectname)
        if filename[-8:-5] == 'twi':
            twi_list.append(keystring)
        for standard in standard_names:
            if standard.lower() in objectname.lower():
                if ifuslot in objectname.lower():
                    standard_list.append(keystring)
        del objectname

twi_file = []
sci_file = []
std_file = []
std_post = []
no_repeats = []
for science_targ in science_target_list:
    date = science_targ.split('_')[0]
    obsid = science_targ.split('_')[1]
    datet = dt(int(date[:4]), int(date[4:6]), int(date[6:]))
    closest_date, diff = find_match(datet, twi_list)
    if len(standard_list):
        closest_date_st, diff_st = find_match(datet, standard_list)
        if diff_st < 1:
            standard_str = ('python panacea/panacea2.py -td %s -to %s -te 1 '
                            '--instr %s --instr_side %s --ifuslot %s -sd %s '
                            '-so %s -rs'
                            % (twi_list[closest_date].split('_')[0],
                               twi_list[closest_date].split('_')[1],
                               args.instrument, args.side, ifuslot,
                               standard_list[closest_date_st].split('_')[0],
                               standard_list[closest_date_st].split('_')[1]))
            std_file.append(standard_str)
            standard_str = ('python panacea/test_fit_lrs2.py --instr %s '
                            '--rootdir %s --side %s -d %s -o %s -e %d'
                            % (args.instrument, '/Users/gregz/cure/reductions',
                               args.side,
                               standard_list[closest_date_st].split('_')[0],
                               standard_list[closest_date_st].split('_')[1],
                               1))
            std_post.append(standard_str)
    if closest_date not in no_repeats:
        no_repeats.append(closest_date)
        twi_panacea_str = ('python panacea/panacea2.py -td %s -to %s -te 1 '
                           '--instr %s --instr_side %s --ifuslot %s -rt'
                           % (twi_list[closest_date].split('_')[0],
                              twi_list[closest_date].split('_')[1],
                              args.instrument, args.side, ifuslot))
        twi_file.append(twi_panacea_str)
    panacea_str = ('python panacea/panacea2.py -td %s -to %s -te 1 --instr %s '
                   '--instr_side %s --ifuslot %s -sd %s -so %s -rs'
                   % (twi_list[closest_date].split('_')[0],
                      twi_list[closest_date].split('_')[1],
                      args.instrument, args.side, ifuslot, date, obsid))
    sci_file.append(panacea_str)
std_file = np.unique(std_file)
std_post = np.unique(std_post)