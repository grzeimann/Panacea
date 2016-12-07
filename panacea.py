#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
IFU Reduction Code 
------------------
Built for the VIRUS instrument as well as LRS2 on HET

Fibers

"""

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

__all__ = ["Panacea"]

import os
import logging
import warnings
import textwrap

import argparse as ap
import numpy as np
import pandas as pd
import os.path as op

from amplifier import Amplifier

from six.moves import configparser

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
    description = textwrap.dedent('''Panacea - 
    
                     This script does ....
                     
                     ''')
                     
    parser = ap.ArgumentParser(description=description,
                            formatter_class=ap.RawTextHelpFormatter)
                        
    parser.add_argument("--specid", nargs='?', type=str, 
                        help='''List of SPECID's for processing. 
                        Ex: "020,008".''', default = None)

    parser.add_argument("--output", nargs='?', type=str, 
                        help='''Output Directory [REQUIRED]''', 
                        default=None)
                        
                        
    parser.add_argument("-s","--subsky", help='''Subtract sky (no arg nec.)''',
                        action="store_true")                        

    parser.add_argument("-rc","--remove_cosmics", 
                        help='''Make new error frame with cosmics set to "-1".''',
                        action="store_true")  

    parser.add_argument("-f","--fiberextract", 
                        help='''Fiberextract (no arg nec.)''',
                        action="store_true") 

    parser.add_argument("-m","--makecube", 
                        help='''Make 3d Cube (no arg nec.)''',
                        action="store_true") 
                        
    parser.add_argument("-d","--detect", 
                        help='''Run detect''',
                        action="store_true") 

    parser.add_argument("--cal_dir", nargs='?', type=str, 
                        help='''Calibration Directory [REQUIRED]''', 
                        default=None)

    parser.add_argument("-rmd","--remake_ditherfiles", 
                        help="Remake the dither files (if they already exist)",
                        action="store_true") 

    parser.add_argument("--instr", nargs='?', type=str, 
                        help='''Instrument to process. 
                        Default: "virus"
                        Ex: "camra" for lab data.''', default = "virus")
                        
    parser.add_argument("--rootdir", nargs='?', type=str, 
                        help='''Root Directory
                        Default: \"/work/03946/hetdex/maverick\"''', 
                        default="/work/03946/hetdex/maverick")

    parser.add_argument("--biasdir", type=str,
                        help= "Directory of biases to use",
                        default="/work/00115/gebhardt/maverick/cal/lib_bias")
                        
    parser.add_argument("-sd","--scidir_date", nargs='?', type=str,
                        help='''Science Directory Date.     [REQUIRED]
                        Ex: \"20160412\"''', default=None)

    parser.add_argument("-so","--scidir_obsid", nargs='?', type=str,
                        help='''Science Directory ObsID.    [REQUIRED]
                        Ex: \"3\" or \"102\"''', default=None)
                        
    parser.add_argument("-se","--scidir_expnum", nargs='?', type=str,
                        help='''Science Directory exposure number.
                        Ex: \"1\" or \"05\"''', default=None) 
                        
    parser.add_argument("--suboverscan_options", nargs="?", type=str, 
                        help='''Set subtract overscan options.
                        Default: \"-s -a -k 2.8 -t -z\"''', 
                        default="-s -a -k 2.8 -t -z")        
 
    parser.add_argument("--skysubtract_options", nargs="?", type=str, 
                        help='''Set sky subtraction options.
                        Default: \"-J -w 250 --output-both\".''', 
                        default="-J -w 250 -T 50 --output-both")

    parser.add_argument("--fiberextract_options", nargs="?", type=str, 
                        help='''Set fiber extract options.
                        Default: \"-n 1032 -W 3500,5500 -P\".''', 
                        default="-n 1032 -W 3500,5500 -P")

    parser.add_argument("--makecube_options", nargs="?", type=str, 
                        help='''Set makecube options.
                        Default: \"\".''', 
                        default="")
 
    parser.add_argument("--detect_options", nargs="?", type=str, 
                        help='''Set detect options.
                        Default: \" -d -S 2 --psf-size 6,6 -c 2.5 -C 2.5 ".''', 
                        default="-d -S 2 --psf-size 6,6 -c 2.5 -C 2.5 ")
                          
    args = parser.parse_args(args=argv)

    # Check that the arguments are filled
    if args.specid:
        args.specid = args.specid.replace(" ", "").split(',')
    else:
        args.specid = SPECID
        
    if args.run_insitu:
        args.run_insitu = True
    else:
        args.run_insitu = False
        
    if args.output is None:
        msg = 'No output directory was provided'
        parser.error(msg)    

    if args.cal_dir is None:
        msg = 'No calibration directory was provided'
        parser.error(msg)    
        
    if args.scidir_date is None:
        msg = 'No science directory date was provided'
        parser.error(msg) 
    else:
        args.scidir_date = args.scidir_date.replace(" ", "").split(',')

    if args.scidir_obsid is None:
        msg = 'No science directory ObsID was provided'
        parser.error(msg) 
    else:
        args.scidir_obsid = args.scidir_obsid.replace(" ", "").split(',')

    if args.scidir_expnum is not None:
        args.scidir_expnum = args.zrodir_expnum.replace(" ", "").split(',')
    
    args.sci_file_loc = []
    for date in args.scidir_date:
        for obsid in args.scidir_obsid:
            if args.scidir_expnum is not None:   
                for expnum in args.scidir_expnum:
                    args.sci_file_loc.append(op.join(args.rootdir, date, 
                                    args.instr, 
                                    "{:s}{:07d}".format(args.instr,int(obsid)), 
                                    "exp{:02d}".format(int(expnum))))
            else:
                args.sci_file_loc.append(op.join(args.rootdir, date, 
                                   args.instr,
                                   "{:s}{:07d}".format(args.instr,int(obsid))))
                
    return args 

class Panacea(object):
    # TODO Just be simple at first.
    """
    A reduction object 
    :param dim:
    """
    def __init__(self, filename, path, config_file):
        self.filename = filename
    
    # ORGANIZE THE FILES, THEN RUN THEM AND PUT THEM IN A LOGICAL PLACE
    # USE THE FILESYSTEM.
    
    
    
    def read_config(self, config_file):
        """Read configuration file
    
        Parameters
        ----------
        configfile : string
            name of the configuration file
    
        Returns
        -------
        config : :class:`~ConfigParser.ConfigParser` instance
        """
        self.config = configparser.ConfigParser(defaults={'xpa_method': 
                                                                      'local'})
        if not self.config.read(config_file):
            msg = ("File '{}' does not exist. Using the default one."
                   " You can get the full set of configuration files with the"
                   " command: ``shuffle_config``".format(config_file))
            warnings.warn(msg)
            config_def = os.path.join(os.path.dirname(__file__), 'configs',
                                      'shuffle.cfg')
            if not self.config.read(config_def):
                msg = ("ERROR: Failed to load the default configuration file.")
                raise ap.ArgumentTypeError(msg.format(config_file))
        
    def setup_logging(args):
        '''Set up a logger for shuffle with a name ``panacea``.
    
        Use a StreamHandler to write to stdout and set the level to DEBUG if
        verbose is set from the command line
        '''
        fmt = '[%(levelname)s - %(asctime)s] %(message)s'
        if args.verbose == 0:
            level = logging.WARNING
        elif args.verbose == 1:
            level = logging.INFO
        else:
            level = logging.DEBUG
            fmt = '[%(levelname)s - %(filename)s - %(asctime)s] %(message)s'
        fmt = logging.Formatter(fmt)
    
        handler = logging.StreamHandler()
        handler.setFormatter(fmt)
        handler.setLevel(level)
    
        log = logging.getLogger('panacea')
        log.setLevel(logging.DEBUG)
        log.addHandler(handler)


        

       