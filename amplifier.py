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
from fiber_utils import get_norm_nonparametric
from fiber import Fiber

__all__ = ["Amplifier"]

class Amplifier:
    def __init__(self, filename, path, calpath=None, debug=False):
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
        if not op.exists(self.path):
            os.mkdir(self.path)
            
        self.N, self.D = F[0].data.shape
        self.D -= 64
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
        self.type = F[0].header['IMAGETYP'].replace(' ', '')
        self.calpath = calpath
        self.debug = debug
        
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
        if not self.fibers:
            fn = op.join(self.path, 'fiber_*_%s.pkl' % self.basename)
            files = sorted(glob.glob(fn))
            for fiber_fn in files:
                if op.exists(fiber_fn):
                    with open(fiber_fn, 'r') as f:
                        self.fibers.append(pickle.load(f))

      
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
    
    def get_trace(self, fdist=2.):
        if self.image is None:
            self.get_image()
        if self.type == 'twi':
            allfibers, xc = get_trace_from_image(self.image, debug=self.debug)
            self.allfibers = allfibers
            brcol = np.argmin(np.abs(xc-self.D*.47))
            standardcol = allfibers[brcol]
            cols1 = xc[brcol::-1]
            cols2 = xc[(brcol+1)::1]
            # Initialize fiber traces
            for i in xrange(len(standardcol)):
                append_flag = False
                try: 
                    F = self.fibers[i]
                except IndexError:    
                    F = Fiber(self.D, i+1, self.path, self.filename)
                    append_flag = True
                if append_flag:
                    self.fibers.append(F)
                self.fibers[i].init_trace_info()
                self.fibers[i].trace_x[brcol] = brcol
                self.fibers[i].trace_y[brcol] = standardcol[i]
            for c in cols1:
                loc = np.where(xc==c)[0]
                for i in xrange(len(standardcol)):
                    yvals = allfibers[loc]
                    xloc = np.argmin(np.abs(self.fibers[i].trace_x - c))
                    floc = np.argmin(np.abs(self.fibers[i].trace_y[xloc] 
                                            - yvals))
                    if (np.abs(self.fibers[i].trace_y[xloc] 
                                                       - yvals[floc]) < fdist):
                        self.fibers[i].trace_x[c] = c
                        self.fibers[i].trace_y[c] = yvals[floc]
            for c in cols2:
                loc = np.where(xc==c)[0]
                for i in xrange(len(standardcol)):
                    yvals = allfibers[loc]
                    xloc = np.argmin(np.abs(self.fibers[i].trace_x - c))
                    floc = np.argmin(np.abs(self.fibers[i].trace_y[xloc] 
                                            - yvals))
                    if (np.abs(self.fibers[i].trace_y[xloc] 
                                                       - yvals[floc]) < fdist):
                        self.fibers[i].trace_x[c] = c
                        self.fibers[i].trace_y[c] = yvals[floc]
            for fiber in self.fibers:
                fiber.fit_trace_poly()
                fiber.eval_trace_poly()

        if self.type == 'sci':
            fn = op.join(self.calpath, 'fiber_*_%s.pkl' % self.basename)
            files = sorted(glob.glob(fn))
            for i, fiber_fn in enumerate(files):
                append_flag = False
                try:
                    F = self.fibers[i]
                except IndexError:    
                    F = Fiber(self.D, i+1, self.path, self.filename)
                    append_flag = True
                with open(fiber_fn, 'r') as f:
                    F1 = pickle.load(f)
                F.trace = F1.trace * 1.
                if append_flag:
                    self.fibers.append(F)
                
                
    def get_fibermodel(self, poly_order=3, use_default=False):
        if self.image is None:
            self.get_image()
        if not self.fibers:
            self.get_trace()
        if self.type == 'twi':
            sol, xcol, binx = fit_fibermodel_nonparametric(self.image, 
                                                           self.fibers,
                                                           debug=self.debug,
                                                       use_default=use_default)
            nfibs, ncols, nbins = sol.shape
            for i, fiber in enumerate(self.fibers):
                fiber.init_fibmodel_info(nbins)
                fiber.fibmodel_poly_order = poly_order
                fiber.fibmodel_x[xcol] = xcol
                fiber.fibmodel_y[xcol,:] = sol[i,:,:]
                fiber.binx = binx
                fiber.fit_fibmodel_poly()
                fiber.eval_fibmodel_poly()
        if self.type == 'sci':
            fn = op.join(self.calpath, 'fiber_*_%s.pkl' % self.basename)
            files = sorted(glob.glob(fn))
            for i, fiber_fn in enumerate(files):
                append_flag = False
                try:
                    F = self.fibers[i]
                except IndexError:    
                    F = Fiber(self.D, i+1, self.path, self.filename)
                    append_flag = True
                with open(fiber_fn, 'r') as f:
                    F1 = pickle.load(f)
                F.fibmodel = F1.fibmodel * 1.
                if append_flag:
                    self.fibers.append(F)        

    def fiberextract(self, poly_order=3, use_default=False):
        if self.image is None:
            self.get_image()
        if not self.fibers:
            self.get_trace()
        if not self.fibers[0].fibmodel:
            self.get_fibermodel(poly_order=poly_order, use_default=use_default)
        norm = get_norm_nonparametric(self.image, self.fibers, 
                                      debug=self.debug)
        for i, fiber in enumerate(self.fibers):
            fiber.spectrum = norm[i,:]
        
               
       