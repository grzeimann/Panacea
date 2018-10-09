# -*- coding: utf-8 -*-
"""
Generate panacea command
"""

import glob
import numpy as np
import os.path as op
import sys
import slurmfile
import fnmatch

from astropy.io import fits
from astropy.table import Table
from datetime import datetime as dt
from input_utils import setup_parser, set_daterange, setup_logging

standard_names = ['HD_19445', 'SA95-42', 'GD50', 'HZ_4', 'G191B2B', 'FEIGE_25',
                  'HILTNER_600', 'G193-74', 'PG0823+546', 'HD_84937',
                  'GD108', 'FEIGE_34', 'HD93521', 'GD140', 'HZ_21',
                  'FEIGE_66', 'FEIGE_67', 'G60-54', 'HZ_44', 'GRW+70_5824',
                  'BD+26+2606', 'BD+33_2642', 'G138-31', 'WOLF_1346',
                  'BD_+17_4708', 'FEIGE_110', 'GD248']


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


def get_exposures(date, obs, args):
    exposures = sorted(glob.glob(op.join(args.rootdir, date, args.instrument,
                                         args.instrument + obs, 'exp*')))
    return [int(fn[-2:]) for fn in exposures]


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
                    help='''Target name or regular expression.  For example,
                    NGC9999*.  The asterix is necessary as the object names differ
                    from the target names as they also include ifuslot and track information.''',
                    type=str, default=None)
parser.add_argument("-s", "--side",
                    help='''red or blue''',
                    type=str, default='red')
parser.add_argument("-mt", "--exposure_time",
                    help='''Mininum Exposure Time''',
                    type=float, default=None)
parser.add_argument("-st", "--standards",
                    help='''set to true if just looking at standards''',
                    type=bool, default=False)
args = parser.parse_args(args=None)

if args.side.lower() == 'blue':
    ifuslot = '056'
    multi_name = 'multi_503_056_7001'
    sides = ['BL', 'BR']
else:
    ifuslot = '066'
    multi_name = 'multi_502_066_7002'
    sides = ['RL', 'RR']
args.log = setup_logging(logname='build_panacea_command')
if args.target is None:
    args.log.error('Please provide "--target" argument on call')
    sys.exit(1)
args = set_daterange(args)
panacea_dict = {}
science_target_list = []
twi_list = []
standard_list = []


object_table = [line.rstrip('\n').split() for line in open('/work/03730/gregz/maverick/object_list_all.dat')]
filenames = []
objects = []
dates = []
keystrings = []
datet = args.daterange[0]
DB = (datet.year, datet.month, datet.day)
datet1 = args.daterange[-1]
DE = (datet1.year, datet1.month, datet1.day)
for _object in object_table:
    if len(_object) == 3:
        filename = _object[0]
        objectname = _object[1]
        exptime = float(_object[2])
    else:
        filename = _object[0]
        exptime = float(_object[-1])
        objectname = ''
        for v in np.arange(1, len(_object)-1):
            objectname = objectname + _object[v]
    obsid = op.basename(op.dirname(op.dirname(op.dirname(filename)))).split(args.instrument)[1]
    date = op.basename(op.dirname(op.dirname(op.dirname(op.dirname(op.dirname(filename))))))
    keystring = date+'_'+obsid
    date_tup = (int(date[:4]), int(date[4:6]), int(date[6:]))
    if DB <= date_tup <= DE:
        if fnmatch.fnmatch(objectname.lower(), args.target.lower()) and filename[-8:-5] == 'sci':
            if args.exposure_time is not None:
                if exptime > args.exposure_time:    
                    science_target_list.append(keystring)
                    print('Science File Found: %s, %s, %0.1f' % (keystring, objectname, exptime))
            else:
                science_target_list.append(keystring)
                print('Science File Found: %s, %s, %0.1f' % (keystring, objectname, exptime))
        if filename[-8:-5] == 'twi':
            twi_list.append(keystring)
        for standard in standard_names:
            if standard.lower() in objectname.lower():
                if ifuslot in objectname.lower():
                    standard_list.append(keystring) 
    
twi_file = []
sci_file = []
com_file = []
std_file = []
std_post = []
no_repeats = []
if args.standards:
    target_list = standard_list
else:
    target_list = science_target_list

for science_targ in target_list:
    date = science_targ.split('_')[0]
    obsid = science_targ.split('_')[1]
    datet = dt(int(date[:4]), int(date[4:6]), int(date[6:]))
    closest_date, diff = find_match(datet, twi_list)
    if len(standard_list):
        closest_date_st, diff_st = find_match(datet, standard_list)
        if diff_st < 1:
            standard_str = ('python /work/03730/gregz/maverick/Panacea/panacea2.py -td %s -to %s -te 1 '
                            '--instr %s --instr_side %s --ifuslot %s -sd %s '
                            '-so %s -rs'
                            % (twi_list[closest_date].split('_')[0],
                               twi_list[closest_date].split('_')[1],
                               args.instrument, args.side, ifuslot,
                               standard_list[closest_date_st].split('_')[0],
                               standard_list[closest_date_st].split('_')[1]))
            std_file.append(standard_str)
            standard_str = ('python Panacea/response_lrs2.py --instr %s '
                            '--rootdir %s --side %s -d %s -o %s -e %d'
                            % (args.instrument, 'reductions', args.side,
                               standard_list[closest_date_st].split('_')[0],
                               standard_list[closest_date_st].split('_')[1],
                               1))
            std_post.append(standard_str)
    if closest_date not in no_repeats:
        no_repeats.append(closest_date)
        twi_panacea_str = ('python /work/03730/gregz/maverick/Panacea/panacea2.py -td %s -to %s -te 1 '
                           '--instr %s --instr_side %s --ifuslot %s -rt'
                           % (twi_list[closest_date].split('_')[0],
                              twi_list[closest_date].split('_')[1],
                              args.instrument, args.side, ifuslot))
        twi_file.append(twi_panacea_str)
    panacea_str = ('python /work/03730/gregz/maverick/Panacea/panacea2.py -td %s -to %s -te 1 --instr %s '
                   '--instr_side %s --ifuslot %s -sd %s -so %s -rs'
                   % (twi_list[closest_date].split('_')[0],
                      twi_list[closest_date].split('_')[1],
                      args.instrument, args.side, ifuslot, date, obsid))
    sci_file.append(panacea_str)
    exps = get_exposures(date, obsid, args)
    for exp in exps:
        panacea_str = ('python /work/03730/gregz/maverick/Panacea/combine_amp_reductions.py -f '
                       'reductions/%s/%s/%s%s/exp%02d/lrs2/%s -s %s -rc'
                       % (date, args.instrument, args.instrument, obsid,
                          exp, multi_name, args.side))
        com_file.append(panacea_str)

std_file = np.unique(std_file)
std_post = np.unique(std_post)
for f, basename in zip([twi_file, sci_file, std_file, std_post,
                        com_file],
                       ['rtwi', 'rsci', 'rstd', 'rresponse', 'rcom']):
    chunks = np.array_split(f, len(f) / 20 + 1)
    for j, chunk in enumerate(chunks):
        n = len(chunk)
        name = basename+'_%s_%i' % (args.side, j+1)
        f = open(name+'.slurm', 'w')
        s = slurmfile.slurmstring % (n, '%j', name)
        f.write(s)
        f.close()
        f = open(name, 'w')
        s = []
        for call in chunk:
            s.append(call)
        f.write('\n'.join(s))
        f.close()
        print('sbatch %s.slurm' % name)
