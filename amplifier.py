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
from fiber_utils import get_norm_nonparametric, check_fiber_trace
from fiber_utils import calculate_wavelength
from fiber import Fiber
from astropy.convolution import Gaussian1DKernel, convolve
from scipy.signal import medfilt

__all__ = ["Amplifier"]

class Amplifier:
    def __init__(self, filename, path, refit=False, calpath=None, debug=False):
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
        self.specid = '%03d' %F[0].header['SPECID']
        self.ifuid = F[0].header['IFUID'].replace(' ', '')
        self.ifuslot ='%03d' %F[0].header['IFUSLOT']
        self.refit = refit
        
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
        
    
    def save_fibers(self):
        for fiber in self.fibers:
            fiber.save(self.specid, self.ifuslot, self.ifuid, self.amp)
       
    def get_image(self):
        image = np.array(fits.open(self.filename)[0].data, dtype=int)
        self.check_overscan(image)
        image[:] = image - self.overscan_value
        self.image = self.orient(image[self.trimsec[2]:self.trimsec[3], 
                                       self.trimsec[0]:self.trimsec[1]])
    
    def get_trace(self, fdist=2., check_trace=True):
        if self.image is None:
            self.get_image()
        if self.type == 'twi' and self.refit:
            allfibers, xc = get_trace_from_image(self.image, interp_window=2.5,
                                                 debug=self.debug)
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
                    yvals = allfibers[int(loc)]
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
                    yvals = allfibers[int(loc)]
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

        else:
            fn = op.join(self.calpath,'fiber_*_%s_%s_%s_%s.pkl' %(self.specid, 
                                                                  self.ifuslot,
                                                                  self.ifuid,
                                                                  self.amp))
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
                F.trace_polyvals = F1.trace_polyvals * 1.
                F.eval_trace_poly()
                if append_flag:
                    self.fibers.append(F)
        if check_trace:
            outfile = op.join(self.path,'trace_%s.pdf' %self.basename)
            check_fiber_trace(self.image, self.fibers, outfile)
            
                
    def get_fibermodel(self, poly_order=3, use_default=False, 
                       make_plots=False):
        if self.image is None:
            self.get_image()
        if not self.fibers:
            self.get_trace()
        if self.type == 'twi' and self.refit:
            sol, xcol, binx = fit_fibermodel_nonparametric(self.image, 
                                                           self.fibers,
                                                           debug=self.debug,
                                                       use_default=use_default,
                                                       plot=make_plots,
                                                       outfolder=self.path,
                                                       fiber_group=8)
            nfibs, ncols, nbins = sol.shape
            for i, fiber in enumerate(self.fibers):
                fiber.init_fibmodel_info(nbins)
                fiber.fibmodel_poly_order = poly_order
                fiber.fibmodel_x[xcol] = xcol
                fiber.fibmodel_y[xcol,:] = sol[i,:,:]
                fiber.binx = binx
                fiber.fit_fibmodel_poly()
                fiber.eval_fibmodel_poly()
        else:
            fn = op.join(self.calpath,'fiber_*_%s_%s_%s_%s.pkl' %(self.specid, 
                                                                  self.ifuslot,
                                                                  self.ifuid,
                                                                  self.amp))
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
                F.fibmodel_polyvals = F1.fibmodel_polyvals * 1.
                F.binx = F1.binx
                F.eval_fibmodel_poly()
                if append_flag:
                    self.fibers.append(F)        


    def fiberextract(self, poly_order=3, use_default_profile=False):
        if self.image is None:
            self.get_image()
        if not self.fibers:
            self.get_trace()
        if not self.fibers[0].fibmodel_polyvals:
            self.get_fibermodel(poly_order=poly_order, 
                                use_default=use_default_profile)
        else:
            for fiber in self.fibers:
                fiber.eval_fibmodel_poly()
        norm = get_norm_nonparametric(self.image, self.fibers, 
                                      debug=self.debug)
        for i, fiber in enumerate(self.fibers):
            fiber.spectrum = norm[i,:]
    
    
    def get_wavelength_solution(self, poly_order=3, wave_order=3, 
                                use_default_profile=False, init_lims=None):
        if self.image is None:
            self.get_image()
        if not self.fibers:
            self.get_trace()
        if not self.fibers[0].fibmodel_polyvals:
            self.get_fibermodel(poly_order=poly_order, 
                                use_default=use_default_profile)
        else:
            for fiber in self.fibers:
                fiber.eval_fibmodel_poly() 
        if not self.fibers[0].spectrum:
            norm = get_norm_nonparametric(self.image, self.fibers, 
                                          debug=self.debug)
            for i, fiber in enumerate(self.fibers):
                fiber.spectrum = norm[i,:]
        if not init_lims:
            print("Please provide initial wavelength endpoint guess")
            sys.exit(1)
        solar_peaks = np.loadtxt('/Users/gregz/cure/panacea/solar_lines.dat')
        solar_spec = np.loadtxt('/Users/gregz/cure/virus_early/virus_config/solar_spec/medium_sun.spec')
        gauss = Gaussian1DKernel(13)
        conv = convolve(solar_spec[:,1],gauss)
        solar_spec[:,1] = medfilt(conv,301)
        sel = ((solar_peaks[:,1]>1.05) * (solar_peaks[:,0]>(init_lims[0]-100.))
                * (solar_peaks[:,0]<(init_lims[1]+100.)))
                
        for i, fiber in enumerate(self.fibers):
            #if i==0:
                fiber.wavelength, fiber.wave_polyvals = calculate_wavelength(
                                             np.arange(self.D), fiber.spectrum,
                                             solar_peaks[sel,:], solar_spec, 
                                             init_lims=init_lims, 
                                             debug=self.debug)
            #else:
            #    fiber.wavelength, fiber.wave_polyvals = calculate_wavelength(
            #                                 np.arange(self.D), fiber.spectrum,
            #                                 solar_peaks[sel,:], init_lims=init_lims, 
            #                                 debug=self.debug, init_sol=self.fibers[i-1].wave_polyvals)
        
