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
import numpy as np
import os
import sys
import re
import glob
import cPickle as pickle
from fiber_utils import get_trace_from_image, fit_fibermodel_nonparametric
from fiber import Fiber

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
        self.filename = op.abspath(filename)
        self.basename = op.basename(filename)[:-5]
        self.path = path            
        self.N,self.D = F[0].data.shape
        self.overscan_value = None
        self.gain = F[0].header['GAIN']
        self.rdnoise = F[0].header['RDNOISE']
        self.amp = (F[0].header['CCDPOS'].replace(' ', '') 
                    + F[0].header['CCDHALF'].replace(' ', ''))
        trim = re.split('[\[ \] \: \,]', F[0].header['TRIMSEC'])[1:-1]
        self.trimsec = [int(t)-((i+1)%2) for i,t in enumerate(trim)]        
        bias = re.split('[\[ \] \: \,]', F[0].header['BIASSEC'])[1:-1]
        self.biassec = [int(t)-((i+1)%2) for i,t in enumerate(bias)]        
        self.fiber_fn = None
        self.fibers = []
        self.image = None
        
    def check_overscan(self, image, recalculate=False):
        if (self.overscan_value is None) or recalculate:
            self.overscan_value = biweight_location(image[
                                              self.biassec[2]:self.biassec[3],
                                              self.biassec[0]:self.biassec[1]])
   
    def save(self):
        fn = op.join(self.path, 'amp_%s.pkl' % self.basename)
        if not op.exists(self.path):
            os.mkdir(self.path)
        with open(fn, 'wb') as f:
           pickle.dump(self, f)

    def load_fibers(self):
        fn = op.join(self.path, 'fiber_*_%s.pkl' % self.basename)
        glob.glob(fn)
        if not op.exists(self.path):
            os.mkdir(self.path)
        with open(fn, 'wb') as f:
           pickle.dump(self, f)
    
    def orient(self, image):
        '''
        Orient the images from blue to red (left to right)
        Fibers are oriented to match configuration files
        '''
        if self.amp == "LU":
            image[:] = image[::-1,::-1]
        if self.amp == "RL":
            image[:] = image[::-1,::-1]
        return image
       
    def get_image(self):
        image = fits.open(self.filename)[0].data
        self.check_overscan(image)
        self.image = self.orient(image[self.trimsec[2]:self.trimsec[3], 
                                       self.trimsec[0]:self.trimsec[1]])
    
    def get_trace(self):
        if self.image is None:
            self.get_image()
        allfibers, xc = get_trace_from_image(self.image)
        try:
            np.vstack(allfibers)
        except ValueError:
            print("The traces didn't find the same number of fibers for each "
                  "column. Going to more complicated mode ...")
            sys.exit(1)
            # TODO make smarter for matching to known number of fibers
        for i, fibery in enumerate(allfibers):
            append_flag = False
            try:
                F = self.fibers[i]
            except IndexError:    
                F = Fiber(self.D, i+1)
                append_flag = True
            F.trace_x[xc] = xc
            F.trace_y[xc] = fibery
            F.fit_trace_poly()
            F.eval_trace_poly()
            if append_flag:
                self.fibers.append(F)
                
    def get_fibermodel(self, poly_order=3):
        if self.image is None:
            self.get_image()
        if not self.fibers:
            self.get_trace()
        sol, xcol, binx = fit_fibermodel_nonparametric(self.image, self.fibers)
        nfibs, ncols, nbins = sol.shape
        for fiber in self.fibers:
            xf = []
            for i in xrange(nbins):
                np.polyfit()
            fiber.binx = binx
            fiber.fibmodel = 
            fiber.fibmodel_poly_order = poly_order
        
        
        
        
       