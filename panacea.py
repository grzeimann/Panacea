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

import argparse as ap
import numpy as np
import pandas as pd
import os.path as op
import cPickle as pickle
import fiber_utils as FU

from utils import biweight_location, biweight_midvariance
from astropy.io import fits, ascii
from fibers import Fibers

from six.moves import configparser

class Panacea(object):
    """
    A reduction object 
    :param dim:
    """
    def __init__(self, filename, path, config_file):
        self.filename = filename
        self.path = path
        self.config = self.read_config(config_file)
        self.overscan = None
        self.fitsfile = fits.open(filename)
        self.get_trimsec_biassec()
        self.check_for_fibers()
        self.image = self.fitsfile[0].data
        self.check_overscan()
        self.amp = None
        
        
    def actions(self):        
        #TODO COMBINE THESE

        self.orient()
        self.divide_pixelflats()

        self.define_fibers()        
        self.get_trace()
        self.get_wavelength()
        self.get_fiberpars()
        
        self.extract_spectra()
        self.get_ra_dec()
        self.find_sources()
    
    def check_for_fibers(self):
        fiber_fn = op.join(self.path, 'fibers_%s' %self.filename[:-5])
        if op.exists(fiber_fn):
            with open(fiber_fn, 'r') as f:
                self.fibers = pickle.read(f)
                self.overscan = self.fibers[0].overscan
        else:
            self.fibers = None
    
    
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

    def orient(self):
        '''
        Orient the images from blue to red (left to right)
        Fibers are oriented to match configuration files
        '''
        if self.amp == "LU":
            self.image[:] = self.image[::-1,::-1]
        if self.amp == "RL":
            self.image[:] = self.image[::-1,::-1]

# scale biases and darks based on sky expectation?
# gather all data necessary for a given measurement
# don't reproduce unnecessary info, just access it, and let the cpu do the work
# Use SAO xpaset ds9 to look at things as an option? Masterbias, subtracted frames?
            
        
    def check_overscan(self, recalculate=False):
        #TODO Make default overscan value: None
        if (self.overscan is None) or recalculate:
            self.overscan = biweight_location(self.image[self.biassec[2]:self.biassec[3],
                                                         self.biassec[0]:self.biassec[1]])
        

       