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
        


            
        
    def check_overscan(self, recalculate=False):
        #TODO define image
        #TODO Make default overscan value: None
        if overscan_value is None:
            overscan = biweight_location(image[by1:by2,bx1:bx2])
            # TODO place overscan somewhere
        elif recalculate:
            overscan = biweight_location(image[by1:by2,bx1:bx2])
            # TODO place overscan somewhere

    def check_dark(self, recalculate=False):
        #TODO define image
        #TODO Make default overscan value: None
        if dark_mult_value is None:
            overscan = biweight_location(image[by1:by2,bx1:bx2])
            # TODO place overscan somewhere
        elif recalculate:
            overscan = biweight_location(image[by1:by2,bx1:bx2])
            # TODO place overscan somewhere            