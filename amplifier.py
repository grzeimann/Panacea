#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Amplifier Class
-----------
To be used in conjuction with IFU reduction code, Panacea


"""

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

from distutils.dir_util import mkpath
from utils import biweight_location, biweight_filter
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
from fiber_utils import calculate_wavelength_chi2, get_model_image
from fiber_utils import check_fiber_profile, check_wavelength_fit
from fiber import Fiber

__all__ = ["Amplifier"]

class Amplifier:
    def __init__(self, filename, path, refit=False, calpath=None, debug=False,
                 darkpath="/Users/gregz/cure/virus_early/virus_config/lib_dark",
                 biaspath="/Users/gregz/cure/virus_early/virus_config/lib_bias",
                 virusconfig="/Users/gregz/cure/virus_early/virus_config/",
                 dark_mult=1., bias_mult=0.):
        ''' 
        Initialize class
        ----------------
        '''
        if not op.exists(filename):
            #TODO log warning message
            print("File Does Not Exist: %s" %filename)
            sys.exit(1)

        F = fits.open(filename)
        self.header = F[0].header
        self.filename = op.abspath(filename)
        self.basename = op.basename(filename)[:-5]
        self.virusconfig = virusconfig
        self.path = path
        if not op.exists(self.path):
            mkpath(self.path)
            
        self.N, self.D = F[0].data.shape
        if self.D == 1064:
            self.D -= 32
        if self.D == 2128:
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
        self.darkpath = darkpath
        self.biaspath = biaspath
        self.dark_mult = dark_mult
        self.bias_mult = bias_mult
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
            fn = op.join(self.path, 'fiber_*_%s_%s_%s_%s.pkl' %(self.specid, 
                                                                  self.ifuslot,
                                                                  self.ifuid,
                                                                  self.amp))
            files = sorted(glob.glob(fn))
            for fiber_fn in files:
                if op.exists(fiber_fn):
                    with open(fiber_fn, 'r') as f:
                        self.fibers.append(pickle.load(f))
                        self.fibers[-1].eval_fibmodel_poly()
                        self.fibers[-1].eval_wave_poly()
                        
    def load_cal_property(self, prop):
        fn = op.join(self.calpath,'fiber_*_%s_%s_%s_%s.pkl' %(self.specid, 
                                                                  self.ifuslot,
                                                                  self.ifuid,
                                                                  self.amp))
        if isinstance(prop, basestring):
            prop = [prop]
        files = sorted(glob.glob(fn))
        for i, fiber_fn in enumerate(files):
            fcheck = op.join(self.calpath,'fiber_%03d_%s_%s_%s_%s.pkl' %(i+1, 
                                                                  self.specid, 
                                                                  self.ifuslot,
                                                                  self.ifuid,
                                                                  self.amp))
            if fcheck != fiber_fn:
                print("Mismatch in loading fiber files from cal directory.")
                print("The %i fiber did not match %s" %(i,fiber_fn))
                sys.exit(1)
            append_flag = False
            try:
                F = self.fibers[i]
            except IndexError:    
                F = Fiber(self.D, i+1, self.path, self.filename)
                append_flag = True
            with open(fiber_fn, 'r') as f:
                F1 = pickle.load(f)
            for pro in prop:
                setattr(F, pro, 1. * getattr(F1, pro))
            if append_flag:
                self.fibers.append(F)
    
    def load_all_cal(self):
        self.load_cal_property(['trace','fibmodel_polyvals',
                                'fibmodel_x','fibmodel_y','binx',
                                'wave_polyvals','fiber_to_fiber'])
        for fiber in self.fibers:
            fiber.eval_fibmodel_poly()
            fiber.eval_wave_poly()
      
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
        image = np.array(fits.open(self.filename)[0].data, dtype=float)
        self.check_overscan(image)
        image[:] = image - self.overscan_value
        image = image[self.trimsec[2]:self.trimsec[3], 
                                       self.trimsec[0]:self.trimsec[1]]
        darkimage = np.array(fits.open(op.join(self.darkpath, 
                                                'masterdark_%s_%s.fits' 
                                                %(self.specid, self.amp)))[0].data, 
                              dtype=float)
        biasimage = np.array(fits.open(op.join(self.biaspath, 
                                                'masterbias_%s_%s.fits' 
                                                %(self.specid, self.amp)))[0].data, 
                              dtype=float)                  
        image[:] = ((image - self.dark_mult * darkimage 
                           - self.bias_mult * biasimage)
                    *self.gain)
        self.image = self.orient(image)
    
    def find_shift(self):
        fn = op.join(self.calpath,'fiber_*_%s_%s_%s_%s.pkl' %(self.specid, 
                                                                  self.ifuslot,
                                                                  self.ifuid,
                                                                  self.amp))
        files = sorted(glob.glob(fn))
        shift = []
        for i, fiber_fn in enumerate(files):
            fcheck = op.join(self.calpath,'fiber_%03d_%s_%s_%s_%s.pkl' %(i+1, 
                                                                  self.specid, 
                                                                  self.ifuslot,
                                                                  self.ifuid,
                                                                  self.amp))
            if fcheck != fiber_fn:
                print("Mismatch in loading fiber files from cal directory.")
                print("The %i fiber did not match %s" %(i,fiber_fn))
                sys.exit(1)
            try:
                F = self.fibers[i]
            except IndexError:    
                print("No trace measured yet, so no shift measured")
                return None
            with open(fiber_fn, 'r') as f:
                F1 = pickle.load(f)
            F1.eval_trace_poly()
            col = self.D/2
            width = 20
            shift.append(biweight_location(F.trace[(col-width):(col+width+1)] 
                              - F1.trace[(col-width):(col+width+1)]))
        return biweight_location(np.array(shift))
                        
       
    def get_trace(self, fdist=2., check_trace=True, calculate_shift=False, 
                  trace_poly_order=3, guess_diff=8.7):
        if self.image is None:
            self.get_image()
        if self.type == 'twi' or self.refit:
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
                self.fibers[i].trace_poly_order = trace_poly_order
                self.fibers[i].init_trace_info()
                self.fibers[i].trace_x[brcol] = brcol
                self.fibers[i].trace_y[brcol] = standardcol[i]
            for c in cols1:
                loc = np.where(xc==c)[0]
                for i in xrange(len(standardcol)):
                    yvals = allfibers[int(loc)]
                    if not yvals.size:
                        continue
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
                    if not yvals.size:
                        continue
                    xloc = np.argmin(np.abs(self.fibers[i].trace_x - c))
                    floc = np.argmin(np.abs(self.fibers[i].trace_y[xloc] 
                                            - yvals))
                    if (np.abs(self.fibers[i].trace_y[xloc] 
                                                       - yvals[floc]) < fdist):
                        self.fibers[i].trace_x[c] = c
                        self.fibers[i].trace_y[c] = yvals[floc]
            for fiber in self.fibers:
                fiber.fit_trace_poly()
                if np.sum(fiber.trace_x != fiber.flag)>(0.9*self.D):
                    fiber.eval_trace_poly()
            for fib, fiber in enumerate(self.fibers):
                if np.sum(fiber.trace_x != fiber.flag)<=(0.9*self.D):
                    sel = fiber.trace_x != fiber.flag
                    k=1
                    done = False
                    while done==False:
                        if fib<(len(self.fibers)-1):
                            if (self.fibers[fib+k].trace == 0).sum()==0:
                                if np.sum(sel)>10:
                                    dif = biweight_location(fiber.trace_y[sel] 
                                               - self.fibers[fib+k].trace_y[sel])
                                else:
                                    dif = guess_diff
                                fiber.trace = self.fibers[fib+k].trace+dif
                                done = True
                            else:
                                k+=1
                        else:
                            if (self.fibers[fib-k].trace == 0).sum()==0:
                                if np.sum(sel)>10:
                                    dif = biweight_location(fiber.trace_y[sel] 
                                               - self.fibers[fib-k].trace_y[sel])
                                else:
                                    dif = guess_diff
                                fiber.trace = self.fibers[fib-k].trace+dif
                                done = True
                            else:
                                k+=1
                        

            if calculate_shift:
                self.net_trace_shift = self.find_shift()
                if self.net_trace_shift is not None:
                    for fiber in self.fibers:
                        fiber.trace_polyvals[-1] += self.net_trace_shift
                        fiber.eval_trace_poly()
            
        else:
            self.load_cal_property('trace_polyvals')
            for fiber in self.fibers:
                fiber.eval_trace_poly()
                
        if check_trace:
            outfile = op.join(self.path,'trace_%s.png' %self.basename)
            check_fiber_trace(self.image, self.fibers, outfile)
            
                
    def get_fibermodel(self, fibmodel_poly_order=3, trace_poly_order=3, 
                       use_default=False, bins=15, 
                       make_ind_plots=False, calculate_shift=False, 
                       check_fibermodel=False):
        if self.image is None:
            self.get_image()
        if not self.fibers:
            self.get_trace(calculate_shift=calculate_shift,
                           trace_poly_order=trace_poly_order)
        if self.type == 'twi' or self.refit:
            sol, xcol, binx = fit_fibermodel_nonparametric(self.image, 
                                                           self.fibers,
                                                           debug=self.debug,
                                                       use_default=use_default,
                                                       plot=make_ind_plots,
                                                       outfolder=self.path,
                                                       fiber_group=8,
                                                       bins=bins)
            nfibs, ncols, nbins = sol.shape
            for i, fiber in enumerate(self.fibers):
                fiber.fibmodel_poly_order = fibmodel_poly_order
                fiber.fibmodel_x = xcol
                fiber.fibmodel_y = sol[i,:,:]
                fiber.binx = binx
                fiber.fit_fibmodel_poly()
                fiber.eval_fibmodel_poly()
        else:
            self.load_cal_property(['fibmodel_polyvals','binx'])
            for fiber in self.fibers:
                fiber.eval_fibmodel_poly()
                
        if check_fibermodel:
            if self.fibers[0].spectrum is None:
                norm = get_norm_nonparametric(self.image, self.fibers, 
                                              debug=self.debug)
                for i, fiber in enumerate(self.fibers):
                    fiber.spectrum = norm[i,:]
            outfile = op.join(self.path,'fibmodel_%s.png' %self.basename)
            check_fiber_profile(self.image, self.fibers, outfile)


    def fiberextract(self, fibmodel_poly_order=3, use_default_profile=False, 
                     calculate_shift=False, trace_poly_order=3):
        if self.image is None:
            self.get_image()
        if not self.fibers:
            self.get_trace(calculate_shift=calculate_shift,
                           trace_poly_order=trace_poly_order)
        if self.fibers[0].fibmodel_polyvals is None:
            self.get_fibermodel(fibmodel_poly_order=fibmodel_poly_order, 
                                use_default=use_default_profile)
        else:
            for fiber in self.fibers:
                fiber.eval_fibmodel_poly()
        norm = get_norm_nonparametric(self.image, self.fibers, 
                                      debug=self.debug)
        for i, fiber in enumerate(self.fibers):
            fiber.spectrum = norm[i,:]
    
    
    def get_wavelength_solution(self, fibmodel_poly_order=3, trace_poly_order=3, 
                                wave_order=3, use_default_profile=False, 
                                init_lims=None, interactive=False, 
                                calculate_shift=False, check_wave=False,
                                check_fibermodel=False, default_fib=0,
                                filt_size_sky=51, filt_size_ind=21):
                                    
        solar_spec = np.loadtxt(op.join(self.virusconfig,'solar_spec/virus_temp.txt'))
        if self.image is None:
            self.get_image()
        if not self.fibers:
            self.get_trace(calculate_shift=calculate_shift,
                           trace_poly_order=trace_poly_order)
        if self.fibers[0].fibmodel_polyvals is None:
            self.get_fibermodel(fibmodel_poly_order=fibmodel_poly_order, 
                                use_default=use_default_profile,
                                check_fibermodel=check_fibermodel)
        else:
            for fiber in self.fibers:
                fiber.eval_fibmodel_poly() 
        if self.type == 'twi' or self.refit:
            if self.fibers[0].spectrum is None:
                norm = get_norm_nonparametric(self.image, self.fibers, 
                                              debug=self.debug)
                for i, fiber in enumerate(self.fibers):
                    fiber.spectrum = norm[i,:]
            if init_lims is None:
                print("Please provide initial wavelength endpoint guess")
                sys.exit(1)
            for k in xrange(2):
                content=False
                cnting = 0
                while content==False:
                    fiber = self.fibers[default_fib]
                    fiber.wavelength, fiber.wave_polyvals = calculate_wavelength_chi2(
                                                     np.arange(self.D), fiber.spectrum, solar_spec, 
                                                     init_lims=init_lims, 
                                                     debug=self.debug, 
                                                     interactive=interactive)
                    # Boundary check:
                    if (np.abs(fiber.wavelength.min()-init_lims[0])>100. 
                          or np.abs(fiber.wavelength.max()-init_lims[1])>100.):
                        default_fib = np.random.choice(112) 
                        cnting+=1
                    else:
                        content=True
                    if cnting>9:
                        print("Ran through 9 loops to find good solution, but couldn't")
                        sys.exit(1)
                fc = np.arange(len(self.fibers))
                if default_fib==0:
                    fibs1 = []
                else:
                    fibs1 = fc[(default_fib-1)::-1]
                if default_fib==(len(self.fibers)-1):
                    fibs2 = []
                else:
                    fibs2 = fc[(default_fib+1)::1]
                for fib in fibs1:
                    if self.debug:
                        print("Working on Fiber %i" %fib)
                    fiber = self.fibers[fib]
                    fiber.wavelength, fiber.wave_polyvals = calculate_wavelength_chi2(
                                                     np.arange(self.D), fiber.spectrum, solar_spec, 
                                                     init_lims=init_lims, 
                                                     debug=self.debug, 
                                                     interactive=False, init_sol=self.fibers[fib+1].wave_polyvals)
                for fib in fibs2:
                    if self.debug:
                        print("Working on Fiber %i" %fib)
                    fiber = self.fibers[fib]
                    fiber.wavelength, fiber.wave_polyvals = calculate_wavelength_chi2(
                                                     np.arange(self.D), fiber.spectrum, solar_spec, 
                                                     init_lims=init_lims, 
                                                     debug=self.debug, 
                                                     interactive=False, init_sol=self.fibers[fib-1].wave_polyvals)
                if k==0:
                    self.get_master_sky(filt_size_sky=filt_size_sky, filt_size_ind=filt_size_ind,
                                norm=True)
                    solar_spec = np.zeros((len(self.masterwave),2))
                    solar_spec[:,0] = self.masterwave
                    solar_spec[:,1] = self.mastersmooth
        else:
            self.load_cal_property(['wave_polyvals'])
            for fiber in self.fibers:
                fiber.eval_wave_poly()
        if check_wave:
            if self.fibers[0].spectrum is None:
                norm = get_norm_nonparametric(self.image, self.fibers, 
                                                  debug=self.debug)
                for i, fiber in enumerate(self.fibers):
                    fiber.spectrum = norm[i,:]
            outfile = op.join(self.path,'wavesolution_%s.png' %self.basename)
            check_wavelength_fit(self.fibers, solar_spec, outfile)
                
    def sky_subtraction(self, fibmodel_poly_order=3, trace_poly_order=3, 
                        wave_order=3, use_default_profile=False, 
                        init_lims=None, interactive=False, 
                        filt_size_ind=21, filt_size_agg=51, 
                        filt_size_final=51, filt_size_sky=51, 
                        calculate_shift=False):
        if self.image is None:
            self.get_image()
        if not self.fibers:
            self.get_trace(calculate_shift=calculate_shift,
                           trace_poly_order=trace_poly_order)
        if self.fibers[0].fibmodel_polyvals is None:
            self.get_fibermodel(fibmodel_poly_order=fibmodel_poly_order, 
                                use_default=use_default_profile)
        else:
            for fiber in self.fibers:
                fiber.eval_fibmodel_poly()
        if self.fibers[0].fiber_to_fiber is None:
            self.load_cal_property(['fiber_to_fiber'])
        if self.fibers[0].spectrum is None:
            self.fiberextract( trace_poly_order= trace_poly_order, 
                              use_default_profile=use_default_profile)
        if self.fibers[0].wave_polyvals is None:
            self.get_wavelength_solution(wave_order=wave_order, 
                                         use_default_profile=use_default_profile, 
                                         init_lims=init_lims,
                                         interactive=interactive)
        else:
            for fiber in self.fibers:
                fiber.eval_wave_poly()
        if self.fibers[0].fiber_to_fiber is None:
            self.get_fiber_to_fiber(wave_order=wave_order, 
                                    use_default_profile=use_default_profile, 
                                    init_lims=init_lims,
                                    interactive=interactive,
                                    filt_size_ind=filt_size_ind,
                                    filt_size_agg=filt_size_agg,
                                    filt_size_final=filt_size_final)
        self.get_master_sky(filt_size_sky=filt_size_sky, filt_size_ind=filt_size_ind,
                            sky=True)
        for fib, fiber in enumerate(self.fibers):
            fiber.sky_spectrum = (fiber.fiber_to_fiber 
                     * np.interp(fiber.wavelength, self.masterwave, self.mastersky))
        self.skyframe = get_model_image(self.image, self.fibers, 'sky_spectrum',
                                        debug=self.debug)
        self.clean_image = self.image - self.skyframe
        
        
    def get_master_sky(self, filt_size_sky=51, filt_size_ind=21, sky=False,
                       norm=False):
        masterwave = []
        if sky:
            masterspec = []
        if norm:
            mastersmooth=[]
        for fib, fiber in enumerate(self.fibers):
            if norm:
                y = biweight_filter(fiber.spectrum, filt_size_ind)
                mastersmooth.append(y/fiber.spectrum)
            masterwave.append(fiber.wavelength)
            if sky:            
                masterspec.append(fiber.spectrum/fiber.fiber_to_fiber)
        masterwave = np.hstack(masterwave)
        ind = np.argsort(masterwave)
        masterwave[:] = masterwave[ind]
        self.masterwave = masterwave
        if sky:
            masterspec = np.hstack(masterspec)
            masterspec[:] = masterspec[ind]
            self.mastersky = biweight_filter(masterspec, filt_size_sky)
        if norm:
            mastersmooth = np.hstack(mastersmooth)
            mastersmooth[:] = mastersmooth[ind]
            self.mastersmooth = biweight_filter(mastersmooth, filt_size_sky)
                    
    def get_fiber_to_fiber(self, fibmodel_poly_order=3, trace_poly_order=3,
                           calculate_shift=False,
                           wave_order=3, use_default_profile=False, 
                           init_lims=None, interactive=False, filt_size_ind=21, 
                           filt_size_agg=51, filt_size_final=51, check_wave=False,
                           check_fibermodel=False):
        if self.image is None:
            self.get_image()
        if not self.fibers:
            self.get_trace(calculate_shift=calculate_shift,
                           trace_poly_order=trace_poly_order)
        if self.fibers[0].fibmodel_polyvals is None:
            self.get_fibermodel(fibmodel_poly_order=fibmodel_poly_order, 
                                use_default=use_default_profile, 
                                check_fibermodel=check_fibermodel)
        else:
            for fiber in self.fibers:
                fiber.eval_fibmodel_poly()
        if self.type == 'twi' or self.refit:
            if self.fibers[0].spectrum is None:
                self.fiberextract()
            if self.fibers[0].wave_polyvals is None:
                self.get_wavelength_solution(wave_order=wave_order, 
                                             use_default_profile=use_default_profile, 
                                             init_lims=init_lims,
                                             interactive=interactive,
                                             check_wave=check_wave)
            else:
                for fiber in self.fibers:
                    fiber.eval_wave_poly()
            masterwave = []
            masterspec = []
            for fib, fiber in enumerate(self.fibers):
                masterwave.append(fiber.wavelength)
                masterspec.append(fiber.spectrum)
            masterwave = np.hstack(masterwave)
            smoothspec = np.hstack(masterspec)
            ind = np.argsort(masterwave)
            masterwave[:] = masterwave[ind]
            smoothspec[:] = smoothspec[ind]
            self.averagespec = biweight_filter(smoothspec, filt_size_agg)
            for fib, fiber in enumerate(self.fibers):
                fiber.fiber_to_fiber = biweight_filter(masterspec[fib] 
                        / np.interp(fiber.wavelength, masterwave, self.averagespec),filt_size_final)

        else:
            self.load_cal_property(['fiber_to_fiber'])         
        

        