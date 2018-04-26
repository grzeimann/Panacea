# -*- coding: utf-8 -*-
"""
Created on Wed Jul 19 14:05:55 2017

@author: gregz
"""

import sys
import os
import re
import datetime
import logging
import glob

import argparse as ap
import os.path as op
import numpy as np
from utils import biweight_location
from datetime import datetime as dt
from astropy.io import fits
from scipy.stats import mode
import matplotlib.pyplot as plt
from matplotlib.dates import MonthLocator, DayLocator, DateFormatter

def parse_args(argv=None):
    # Arguments to parse include date range, inspection attribute, instrument,
    #           instrument_side
    parser = ap.ArgumentParser(description="DIAGNOSIS",
                            formatter_class=ap.RawTextHelpFormatter)
                            
    parser.add_argument("-sd","--start_date", 
                        help='''Start Date, e.g., 20170321''',
                        type=str, default=None)

    parser.add_argument("-ed","--end_date", 
                        help='''Start Date, e.g., 20170326''',
                        type=str, default=None)

    parser.add_argument("-dl","--date_length", 
                        help='''Days after/before start/end date, e.g., 10''',
                        type=int, default=None)

    parser.add_argument("-r","--rootdir", 
                        help='''Root Directory for Date''',
                        type=str, default='/work/03946/hetdex/maverick')

    parser.add_argument("-i","--ifuslot", 
                        help='''IFUSLOT, e.g., 076''',
                        type=str, default=None)

    parser.add_argument("-in","--instrument", 
                        help='''Instrument, e.g., virus''',
                        type=str, default='virus')

    parser.add_argument("-cb","--check_bias", 
                        help='''Check biases over date range''',
                        action="count", default=0)

    parser.add_argument("-cd","--check_dark", 
                        help='''Check darks over date range''',
                        action="count", default=0)

    parser.add_argument("-ct","--check_trace", 
                        help='''Check trace over date range''',
                        action="count", default=0)

    parser.add_argument("-cf","--check_fiberprofile", 
                        help='''Check fiber profile over date range''',
                        action="count", default=0)
    
    args = parser.parse_args(args=argv)
    
    args.log = setup_logging()
    
    if args.ifuslot is None:
        args.log.error('Please set the IFUSLOT you would like to diagnose.')
        sys.exit(1)
    
    args.instrument = args.instrument.lower()
    if args.instrument == 'virus':
        args.amps = ['LL','LU','RL','RU']
    if args.instrument == 'lrs2':
        args.amps = ['LL','LU','RL','RU']
    if args.instrument == 'virusw':
        args.amps = ['RU']

    args = set_daterange(args)                  
                              
    return args

def setup_logging():
    '''Set up a logger for shuffle with a name ``diagnosis``.

    Use a StreamHandler to write to stdout and set the level to DEBUG if
    verbose is set from the command line
    '''
    log = logging.getLogger('diagnosis')
    if not len(log.handlers):
        fmt = '[%(levelname)s - %(asctime)s] %(message)s'       
        fmt = logging.Formatter(fmt)
    
        level = logging.INFO

        handler = logging.StreamHandler()
        handler.setFormatter(fmt)
        handler.setLevel(level)
    
        log = logging.getLogger('diagnosis')
        log.setLevel(logging.DEBUG)
        log.addHandler(handler)
    return log

def set_daterange(args):
    dateatt = ['start_date', 'end_date']
    if args.date_length is None:
        if args.start_date is None:
            args.log.error('You must include two of the following: '
                           '"start_date", "end_date", or "date_length"')
            sys.exit(1)
        if args.end_date is None:
            args.log.error('You must include two of the following: '
                           '"start_date", "end_date", or "date_length"')
            sys.exit(1)
        dates = {}
        for da in dateatt:
            dates[da] = dt(int(getattr(args,da)[:4]),
                           int(getattr(args,da)[4:6]),
                           int(getattr(args,da)[6:]))
        
        args.daterange=[datetime.date.fromordinal(i) 
                        for i in range(dates[dateatt[0]].toordinal(), 
                                       dates[dateatt[1]].toordinal())]
    else:
        if args.start_date is not None and args.end_date is not None:
            args.log.warning('Using "start_date" and "date_length", '
                             'however, you specified "end_date" as well '
                             'which will not be used.')
            args.end_date = None
        if args.start_date is not None:
            base = datetime(int(args.start_date[:4]),int(args.start_date[4:6]),
                            int(args.start_date[:6]))
            args.daterange = [base + datetime.timedelta(days=x) 
                              for x in range(0, args.date_length)]

        if args.end_date is not None:
            base = datetime(int(args.end_date[:4]),int(args.end_date[4:6]),
                            int(args.end_date[:6]))
            args.daterange = [base - datetime.timedelta(days=x) 
                              for x in range(0, args.date_length)] 
                                  
    return args

def check_trace(args):
    pass

def check_fiberprofile(args):
    pass

def check_bias(args):
    for amp in args.amps:
        date_x = []
        dettemp = []
        ambtemp = []
        bias_y = []
        args.log.info('Looking at amp: %s' %amp)
        for fn in args.bias_list:
            if len(fn.split(amp+'_'))>1:
                try:    
                    F = fits.open(fn)
                    flag = True
                except:
                    args.log.warning('%s does exist, but must be corrupted.' %fn)
                    flag = False
                if flag:
                    flag = check_quality(args, F, fn)
                if flag:
                    biassec = re.split('[\[ \] \: \,]', 
                                       F[0].header['BIASSEC'])[1:-1]
                    biassec = [int(t)-((i+1)%2) for i,t in enumerate(biassec)]
                    data = F[0].data[biassec[2]:biassec[3], biassec[0]:biassec[1]]
                    dettemp.append(F[0].header['DETTEMP'])
                    ambtemp.append(F[0].header['AMBTEMP'])
                    bias_y.append(np.mean(data))
                    base = op.basename(fn)
                    Y, M, D, H, m, S = [int(s) for s in [base[:4],base[4:6],
                                                     base[6:8],base[9:11],
                                                     base[11:13],base[13:15]]]
                    date_x.append(datetime.datetime(Y, M, D, H, m, S))
                    del F[0].data
                    F.close()
        make_dateplot(date_x, bias_y, 
                      'bias_diagnositic_%s_%s.png'%(args.ifuslot,amp),
                      ylims=[900.,1100.])
        make_plot(dettemp, bias_y, 
                  'bias_diagnositic_dettemp_%s_%s.png'%(args.ifuslot,amp),
                  xlims=[-120.,-80.],
                  ylims=[900.,1100.])
        make_plot(ambtemp, bias_y, 
                  'bias_diagnositic_ambtemp_%s_%s.png'%(args.ifuslot,amp),
                  xlims=[0.,30.],
                  ylims=[900.,1100.])
                  
def check_dark(args):
    pass

def get_files(args):
    args.trace_list = []
    args.bias_list = []
    args.dark_list = []
    args.log.info('Looking at %i days' %len(args.daterange))
    for date in args.daterange:
        datestr = '%04d%02d%02d' %(date.year, date.month, date.day)
        files = glob.glob(op.join(args.rootdir, datestr, args.instrument,'*',
                          'exp01', args.instrument,'2*_%s*.fits' %args.ifuslot))
        for fn in files:
            imtype = fn[-8:-5]
            if imtype in ['sci','twi']:
                args.trace_list.append(fn)
            if imtype in ['zro']:
                args.bias_list.append(fn)
            if imtype in ['drk']:
                args.dark_list.append(fn)
                    
    return args
    
def ensure_no_stuckbits(F, args, fn): 
    bias = re.split('[\[ \] \: \,]', F[0].header['BIASSEC'])[1:-1]
    biassec = [int(t)-((i+1)%2) for i,t in enumerate(bias)]
    data = F[0].data[biassec[2]:biassec[3], biassec[0]:biassec[1]]
    mode_value = int(mode(data.ravel()).mode)
    missing_values = []
    for i in np.arange(mode_value-8, mode_value+9):
        if (data == i).sum() == 0:
            missing_values.append(i)
    for missing in missing_values:
        args.log.warning('The value %i is not represented '
                         'in the overscan region for %s' %(missing, fn))
    if len(missing_values):
        return False
    else:
        return True
    
def check_not_all_zeros(F, args, fn):
    bias = re.split('[\[ \] \: \,]', F[0].header['BIASSEC'])[1:-1]
    biassec = [int(t)-((i+1)%2) for i,t in enumerate(bias)]
    data = F[0].data[biassec[2]:biassec[3], biassec[0]:biassec[1]]
    if (data == 0).sum() == (data.shape[0]*data.shape[1]):
        return False
    else:
        return True

def check_quality(args, F, filename):
    flag = True
    if not check_not_all_zeros(F, args, op.basename(filename)):
        flag = False
        args.log.info('%s is all zeros' %filename)
    if flag:
        if not ensure_no_stuckbits(F, args, op.basename(filename)):
            flag = False
            args.log.info('%s has stuck bits' %filename)

    return flag    

def make_dateplot(x, y, outname, ylims=None):
    fig, ax = plt.subplots()
    fig.autofmt_xdate()
    ax.plot_date(x, y, ls='', marker='x')
    ax.xaxis.set_major_locator(MonthLocator())
    ax.xaxis.set_minor_locator(DayLocator())
    ax.xaxis.set_major_formatter(DateFormatter('%Y-%m-%d'))
    ax.fmt_xdata = DateFormatter('%Y-%m-%d %H:%M:%S')
    if ylims is not None:
        ax.set_ylim(ylims)
    fig.savefig(outname, dpi=200)
    plt.close(fig)

def make_plot(x, y, outname, xlims=None, ylims=None):
    fig, ax = plt.subplots()
    fig.autofmt_xdate()
    ax.plot(x, y, ls='', marker='x')
    if ylims is not None:
        ax.set_ylim(ylims)
    if xlims is not None:
        ax.set_xlim(xlims)
    fig.savefig(outname, dpi=200)
    plt.close(fig)
   
def main():
    '''
    Diagnose the dark current, bias levels, fiber trace, or fiber profile    
    '''
    args = parse_args()
    
    args = get_files(args)
    
    
    if args.check_bias:
        args.log.info('Investigating bias for %s' %args.ifuslot)
        check_bias(args)
        args.log.info('Done investigating bias for %s' %args.ifuslot)
    
    if args.check_dark:
        args = check_quality(args,'dark_list')
        check_dark(args)
    
    if args.check_trace:
        args = check_quality(args,'trace_list')
        check_trace(args)
        
    if args.check_fiberprofile:
        args = check_quality(args,'trace_list')
        check_fiberprofile(args)
    pass

if __name__ == '__main__':
    main() 
    
# Line up files to do all checks that I want to do for a given file, then
# discard file/amplifier.
   
# Twi checks, versus bias checks, versus dark checks, versus science checks
