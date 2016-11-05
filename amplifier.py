#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Amplifier Class
-----------
To be used in conjuction with IFU reduction code, Panacea


"""

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

from utils import biweight_location
from astropy.io import fits
import os.path as op
import sys
import re

__all__ = ["Amplifier"]

class Amplifier:
    def __init__(self, filename, path):
        ''' 
        Initialize class
        ----------------
        '''
        if not op.exists(filename):
            #TODO log warning message
            sys.exit(1)

        F = fits.open(filename)
        self.filename = filename
        self.basename = op.basename(filename)[:-5]
        self.path            
        self.N,self.D = F[0].data.shape
        self.overscan_value = None
        self.gain = F[0].header['GAIN']
        self.rdnoise = F[0].header['RDNOISE']
        self.amp = (F[0].header['CCDPOS'].replace(' ', '') 
                    + F[0].header['CCDHALF'].replace(' ', ''))
        self.trimsec = re.split('[\[ \] \: \,]', F[0].header['TRIMSEC'])[1:-1]
        self.biassec = re.split('[\[ \] \: \,]', F[0].header['BIASSEC'])[1:-1]
        self.fibers = None
        
    def check_overscan(self, recalculate=False):
        if (self.overscan_value is None) or recalculate:
            self.overscan_value = biweight_location(self.image[
                                              self.biassec[2]:self.biassec[3],
                                              self.biassec[0]:self.biassec[1]])
   
    def save(self):
        fn = op.join(path, 'amp_%.pkl' % basename)
        with open(fn, 'wb') as f:
           pickle.dump(f, self)
       