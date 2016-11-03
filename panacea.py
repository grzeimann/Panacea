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

import numpy as np
import pandas as pd
from utils import biweight_location, biweight_midvariance
from astropy.io import fits, ascii
import logging
import fiber_utils

class Panacea(object):
    """
    A reduction object 
    :param dim:
    """
    def __init__(self, filename, args):
        self.filename = filename
        self.image = fits.open(filename)[0].data
        self.amp = None
        self.args = args
        
        
    def actions(self):        
        #TODO COMBINE THESE
        
        self.join_amps()
        self.divide_pixelflats()

        self.define_fibers()        
        self.get_trace()
        self.get_wavelength()
        self.get_fiberpars()
        
        self.extract_spectra()
        self.get_ra_dec()
        self.find_sources()
        
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

# IF then structure
# scale biases and darks based on sky expectation?
# gather all data necessary for a given measurement
# check if evaluated, if not then evaluate
# don't reproduce unnecessary info, just access it, and let the cpu do the work
# Use SAO xpaset ds9 to look at things as an option? Masterbias, subtracted frames?
            
        
    def check_overscan(fiber, image, recalculate=False):
        #TODO define image
        #TODO Make default overscan value: None
        if fiber.overscan is None:
            fiber.overscan = biweight_location(image[fiber.biassec[2]:fiber.biassec[3],fiber.biassec[0]:fiber.biassec[1]])
            # TODO place overscan somewhere
        elif recalculate:
            fiber.overscan = biweight_location(image[by1:by2,bx1:bx2])
            # TODO place overscan somewhere
       