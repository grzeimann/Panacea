#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
IFU Reduction Code 
------------------
Built for the VIRUS instrument as well as LRS2 on HET


"""

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

__all__ = ["Panacea"]

import numpy as np
import pandas as pd
from utils import biweight_location, biweight_midvariance
from astropy.io import fits, ascii

class Panacea(object):
    """
    A reduction object 
    :param dim:
    """
    def __init__(self, args=[], kwargs={}):
        self.args = args
        self.kwargs = kwargs
        
    def model(self):
        model = bias + (dark * self.time + fiber_weight * spectrum) * gain 
        
    def actions(self):
        self.remove_overscan()
        
        #TODO COMBINE THESE
        self.remove_dark()
        self.remove_residual_bias()
        
        self.join_amps()
        self.divide_pixelflats()

        self.define_fibers()        
        self.get_trace()
        self.get_wavelength()
        self.get_fiberpars()
        
        self.extract_spectra()
        self.get_ra_dec()
        self.find_sources()
        

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

    def check_dark(self, recalculate=False):
        #TODO define image
        #TODO Make default overscan value: None
        if fiber.dark_mult is None:
            fiber.dark_mult = biweight_location(image[by1:by2,bx1:bx2])
            # TODO place overscan somewhere
        elif recalculate:
            fiber.dark_mult = biweight_location(image[by1:by2,bx1:bx2])
            # TODO place overscan somewhere            