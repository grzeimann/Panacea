#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Characterization script
------------------
Built for characterizing VIRUS instrument as well as LRS2 on HET

Incomplete Documentation

"""

import pandas as pd
import argparse as ap
from distutils.dir_util import mkpath
import textwrap
import glob
import os.path as op
from astropy.io import fits
from amplifier import Amplifier

def parse_args(argv=None):
    """Parse the command line arguments

    Parameters
    ----------
    argv : list of string
        arguments to parse; if ``None``, ``sys.argv`` is used

    Returns
    -------
    Namespace
        parsed arguments
    """
    description = textwrap.dedent('''Characterize - 
    
                     This script does ... (Fill in Later)
                     
                     ''')
                     
    parser = ap.ArgumentParser(description=description,
                            formatter_class=ap.RawTextHelpFormatter)
                        
    parser.add_argument("--specid", nargs='?', type=str, 
                        help='''List of SPECID's for processing. [REQUIRED]
                        Ex: "020,008".''', default = None)

    parser.add_argument("--instr", nargs='?', type=str, 
                        help='''Instrument to process. 
                        Default: "virus"
                        Ex: "camra" for lab data,
                            "lrs2" for lrs2.''', default = "camra")

    parser.add_argument("--output", nargs='?', type=str, 
                        help='''Output Directory
                        Default: \"characterized"''', 
                        default="characterized")
                        
    parser.add_argument("--rootdir", nargs='?', type=str, 
                        help='''Root Directory
                        Default: \"/work/03946/hetdex/maverick\"''', 
                        default="/work/03946/hetdex/maverick")

    parser.add_argument("-bd","--biadir_date", nargs='?', type=str,
                        help='''Bias Directory Date.    
                        Ex: \"20160412\"''', default=None)

    parser.add_argument("-bo","--biadir_obsid", nargs='?', type=str,
                        help='''Bias Directory ObsID.    
                        Ex: \"3\" or \"102\"''', default=None)
                        
    parser.add_argument("-be","--biadir_expnum", nargs='?', type=str,
                        help='''Science Directory exposure number.
                        Ex: \"1\" or \"05\"''', default=None) 

    parser.add_argument("-dd","--drkdir_date", nargs='?', type=str,
                        help='''Dark Directory Date.   
                        Ex: \"20160412\"''', default=None)

    parser.add_argument("-do","--drkdir_obsid", nargs='?', type=str,
                        help='''Dark Directory ObsID.  
                        Ex: \"3\" or \"102\"''', default=None)
                        
    parser.add_argument("-de","--drkdir_expnum", nargs='?', type=str,
                        help='''Science Directory exposure number.
                        Ex: \"1\" or \"05\"''', default=None) 

    parser.add_argument("-d","--debug", help='''Debug.''',
                        action="count", default=0)
                          
    args = parser.parse_args(args=argv)

    # Check that the arguments are filled
    if args.specid:
        args.specid = args.specid.replace(" ", "").split(',')
    else:
        msg = 'No SPECID was provided.'
        parser.error(msg)   

    labels = ['dir_date', 'dir_obsid', 'dir_expnum']
    observations=['bia', 'drk', 'pxf']
    for obs in observations:
        amp_list = []
        for label in labels[:2]:
            getattr(args, obs+label)
            if getattr(args, obs+label) is None:
                msg = '%s%s was not provided' %(obs, label)
                parser.error(msg) 
            else:
                setattr(args, obs+label, 
                        getattr(args, obs+label).replace(" ", "").split(','))
        if getattr(args, obs+labels[2]) is not None:
            setattr(args, obs+labels[2], 
                    getattr(args, obs+labels[2]).replace(" ", "").split(','))
        for date in getattr(args, obs+labels[0]):
            for obsid in getattr(args, obs+labels[1]):
                if getattr(args, obs+labels[2]) is not None:   
                    for expnum in getattr(args, obs+labels[2]):
                        folder = op.join(date, 
                                         args.instr, 
                                         "{:s}{:07d}".format(args.instr,
                                         int(obsid)), 
                                         "exp{:02d}".format(int(expnum)),
                                         args.instr)
                        files = sorted(glob.glob(op.join(args.rootdir, folder, 
                                                         '*')))
                        for fn in files:
                            amp_list.append(Amplifier(fn, '', name=obs))
                else:
                    folder = op.join(date, args.instr,
                                     "{:s}{:07d}".format(args.instr, 
                                                         int(obsid)))
                    files = sorted(glob.glob(op.join(args.rootdir, folder, '*', 
                                                     args.instr, '*')))
                    for fn in files:
                        amp_list.append(Amplifier(fn, '', name=obs))
        setattr(args,obs+'_list', amp_list)

    return args       


def measure_gain():
    # gain code from utility
    # print to file
    pass

def measure_readnoise():
    pass
    
def make_pixelflats():
    pass

def check_pixelflats():
    pass

def check_bias():
    pass

def check_darks():
    pass

def check_arcs():
    pass

def get_trace():
    pass

def relative_throughput():
    pass

def main():
    args = parse_args()
    
    pass

if __name__ == '__main__':
    main()  