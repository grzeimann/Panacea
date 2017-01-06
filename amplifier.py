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
import cosmics


__all__ = ["Amplifier"]

class Amplifier:
    def __init__(self, filename, path, refit=False, calpath=None, skypath=None,
                 debug=False, darkpath=None, biaspath=None, virusconfig=None,
                 dark_mult=1., bias_mult=0., use_pixelflat=True, 
                 specname=None, fdist=2., check_trace=True, 
                 calculate_shift=False, trace_poly_order=3,
                 fibmodel_poly_order=3, use_default_fibmodel=False, 
                 fibmodel_nbins=15, make_fib_ind_plots=False, 
                 check_fibermodel=False, fsize=8., sigma=2.5, power=2.5,
                 fiber_group=8, col_group=48, mask=None, wave_nbins=21, 
                 wave_order=3, default_fib=0, init_lims=None, 
                 interactive=False, check_wave=False,filt_size_ind=21, 
                 filt_size_agg=51, filt_size_final=51, filt_size_sky=51):
        ''' 
        Initialize class
        ----------------
        :param filename:
            Filename for the raw data (which is just a single amplifier)
        :param path:
            The output path where to save the amplifier class and Fiber class
        :param refit:
            If set to yes, then nearly all attributes are refit.
            Specifically:
                1) Trace
                2) Fibermodel
                3) Wavelength Solution
                4) Fiber to Fiber Normalization for the amplifier in question
        :param calpath:
            This parameter sets the file path that points to the compressed
            Fiber Class files of previously reduced calibration frames.  
            This path is used to load trace, wavelength solution, fibermodel, 
            and fiber_to_fiber properties.
        :param skypath:
            This parameter sets the file path that points to the compressed
            Fiber Class files of previously reduced sky frames.  This path
            is used to load "sky_spectrum" properties.
        :param debug:
            This is used to print out run times for different stages.  In many
            ways this parameter is ever changing for various debugging 
            purposes.
        :param darkpath:
            This path is used to point to the library containing dark current
            images for each amplifier.
        :param biaspath:
            This path is used to point to the library containing bias images
            for each amplifier.
        :param virusconfig:
            This path points to a more general library of PixelFlats, IFUcen
            Files, Fplane files, and solar/twighlight spectra.
        :param dark_mult:
            This value is used to multiply the dark image before subtraction.
        :param bias_mult:
            This value is used to multiply the bias image before subtraction.
        :param use_pixelflat:
            A boolean for correcting the pixel to pixel variations.
        :param specname:
            Base filename for loading the solar/twighlight spectrum used
            in wavelength calibration.
        :param fdist:
            This parameter sets the distance allowable to match fiber
            centroid measurments for neighboring columns.  Ideal values are
            generous enough for glitches and gaps and yet small enough to not 
            jump to neighboring fibers.  
        :param check_trace:
            Plot trace images at the end.  Good for post checks.
        :param calculate_shift:
            Calculate the shift between calibration files and new trace 
            measurement.
        :param trace_poly_order:
            Not in use now, but remains for possible future use.
        :param fibmodel_poly_order:
            This parameter is not used anymore but remains for flexibility if
            it is useful again.
        :param use_default_fibmodel:
            If true, no fit is performed and instead the default fiber model
            is used.  The default profile is a gauss-hermite exponential and is
            defined by "bins", "fsize", "sigma", and "power". 
        :param fibmodel_nbins:
            The number of bins (+2) used to define the fiber profile.
        :param make_fib_ind_plots:
            If True, plots for each image stamp that is fit.  Produces 
            number of fibers times col_group images.  Typically 112x48 = 5376.
            This also takes awhile to save all of these images.
        :param check_fibermodel:
            If True, a summary plot is produced of the fibermodel fit for the 
            top/bottom/middle of the amplifier in fibers and columns 
            (3x3 plot).
        :param fsize:
            The pixel size from the center of the trace over which the 
            fibermodel is defined. In other words: [trace-fsize:trace+fsize].
        :param sigma:
            The sigma parameter for the gauss-hermite exponential for the 
            initial/default fibermodel.  This model is used to defining the
            binning of the empirical fibermodel through the CDF of the 
            second derivative.
        :param power:
            The power parameter for the gauss-hermite exponential for the 
            initial/default fibermodel.
        :param fiber_group:
            Total number of fibers used to constrain the profile at one time.
        :param col_group:
            Total number of columns used to constrain the profile at one time.
        :param mask:
            Used for masking pixels and avoids them in the spectral extraction.
        :param wave_nbins:
            Number of bins used across the wavelength range to fit
            linear solutions.
        :param wave_order:
            Polynomial order for fitting the wavelength solution.
        :param default_fib:
            The first fiber to fit for the wavelength solution.
        :param init_lims:
            List or tuple of two values for the starting and finishing 
            wavelengths for the amplifier.  General values are fine.  In
            other words, [3500,5500] for VIRUS.
        :param interactive:
            The first fiber wavelength solution fit is interactive.
        :param check_wave:
            Plots the wavelength solution in a 3x3 grid covering the 
            top/middle/bottom fibers and columns. 
        :param filt_size_ind:
            Biweight filter size for individual spectra
        :param filt_size_agg:
            Biweight filter size for the master spectrum
        :param filt_size_final:
            Biweight filter size for fiber to fiber
        :param filt_size_sky:
            Biweigth filter size for the master sky
        
        :init header:
            The fits header of the raw frame.
        :init basename:
            The basename for the filename input parameter.
        :init N:
            The number of rows (y-dir)
        :init D:
            The number of cols (x-dir) minus the overscan region
        :init overscan_value:
            Defaults to none but is evaluated as the biweight_location
            of the values in the overscan region
        :init gain:
            (e-/ADU) read from the header
        :init rdnoise:
            (e-) read from the header
        :init amp:
            The amplifier name you are looking at, e.g., "RU" or "LL"
        :init trimsec:
            The CCD section defining the DATASEC.
        :init biassec:
            The CCD section defining the BIASSEC. (overscan region)
        :init fibers:
            Initially an empty list to be filled with Fiber Class Objects,
            one for each fiber in the amplifier.
        :init image:
            The trimmed image that will be manipulated and cleaned during the
            reduction process.
        :init type:
            Observation type.
        :init specid:
            The spectrograph ID for this amplifier.
        :init ifuid:
            The IFU ID for the fiber bundle.
        :init ifuslot:
            The IFU slot position in the Fplane array.
        '''
        # Organization options
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
        self.refit = refit
        self.calpath = calpath
        self.skypath = skypath
        self.darkpath = darkpath
        self.biaspath = biaspath
        self.specname = specname
        self.debug = debug
        
        # Image manipulation options
        self.dark_mult = dark_mult
        self.bias_mult = bias_mult
        self.use_pixelflat = use_pixelflat
        
        # Trace options
        self.fdist = fdist
        self.check_trace = check_trace
        self.trace_poly_order = trace_poly_order
        self.calculate_shift = calculate_shift
        
        # Fibermodel options
        self.fibmodel_poly_order = fibmodel_poly_order
        self.use_default_fibmodel = use_default_fibmodel
        self.fibmodel_nbins = fibmodel_nbins
        self.make_fib_ind_plots = make_fib_ind_plots
        self.check_fibermodel = check_fibermodel
        self.fsize = fsize
        self.sigma = sigma
        self.power = power
        self.fiber_group = fiber_group
        self.col_group = col_group
        
        # Masking options (Fiberextract related)
        self.mask = mask
        
        # Wavelength Solution options
        self.wave_nbins = wave_nbins
        self.wave_order = wave_order
        self.default_fib = default_fib
        self.init_lims = init_lims
        self.interactive = interactive
        self.check_wave = check_wave        
        
        # Smoothing options for individual and master spectra
        self.filt_size_ind = filt_size_ind
        self.filt_size_agg = filt_size_agg
        self.filt_size_final = filt_size_final
        self.filt_size_sky = filt_size_sky
        
        # Initialized variables
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
        self.fibers = []
        self.image = None
        self.type = F[0].header['IMAGETYP'].replace(' ', '')
        self.specid = '%03d' %F[0].header['SPECID']
        self.ifuid = F[0].header['IFUID'].replace(' ', '')
        self.ifuslot ='%03d' %F[0].header['IFUSLOT']

        
    def check_overscan(self, image, recalculate=False):
        '''
        Evaluate the overscan_value using a biweight average of the BIASSEC
        region.  Only calculate the value if one does not exist or 
        recalculate is set to True.
        '''
        if (self.overscan_value is None) or recalculate:
            self.overscan_value = biweight_location(image[
                                              self.biassec[2]:self.biassec[3],
                                              self.biassec[0]:self.biassec[1]])
   
   
    def save(self):
        '''
        Save the entire amplifier include the list of fibers.  
        This property is not used often as "amp*.pkl" is large and typically
        the fibers can be loaded and the other amplifier properties quickly
        recalculated.
        '''
        fn = op.join(self.path, 'amp_%s.pkl' % self.basename)
        if not op.exists(self.path):
            os.mkdir(self.path)
        with open(fn, 'wb') as f:
           pickle.dump(self, f)
           

    def load_fibers(self):
        '''
        Load fibers in self.path.  Redefine the path for the fiber to self.path
        in case it was copied over. After loaded, evaluate the fibermodel as 
        well as the wavelength solution if available.  
        '''
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
                        self.fibers[-1].path = self.path
                        if self.fibers[-1].wave_polyvals is not None:
                            self.fibers[-1].eval_wave_poly()
    

    def convert_binning(self, fiber, prop):
        '''
        Occassionally the binning of property being loaded from calibration
        or otherwise is not binned the same as the data.  This function 
        converts between different binning if necessary.  The properties
        affected by this are "fibmodel_x", "trace", and "fiber_to_fiber".  
        Others are immune due to their evaluation methods.
        :param fiber:
            A fiber class object from which you are loading "prop".
        :param prop:
            A property to load from "fiber".
        '''
        if prop == 'fibmodel_x':
            values = getattr(fiber, prop)
            return (1.* values) / fiber.D * self.D
        elif prop in ['trace', 'fiber_to_fiber']:
            values = getattr(fiber, prop)
            return np.interp((1.*np.arange(self.D))/self.D, 
                             (1.*np.arange(fiber.D))/fiber.D, values)
            
        else:
            return 1.* getattr(fiber, prop)       
                             
                             
    def load_cal_property(self, prop, pathkind='calpath'):
        '''
        Load a specific property from a calibration object.  This can be
        "trace" from a twighlight frame or "sky_spectrum" from another
        science frame.  
        :param prop:
            A property to load from a set of fibers.
        :param pathkind:
            Used to defined the path from the larger class.  E.g., "calpath"
            or "skypath".
        '''
        path = getattr(self, pathkind)
        fn = op.join(path,'fiber_*_%s_%s_%s_%s.pkl' %(self.specid, 
                                                                  self.ifuslot,
                                                                  self.ifuid,
                                                                  self.amp))
        if isinstance(prop, basestring):
            prop = [prop]
        files = sorted(glob.glob(fn))
        for i, fiber_fn in enumerate(files):
            fcheck = op.join(path,'fiber_%03d_%s_%s_%s_%s.pkl' %(i+1, 
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
                try:
                    setattr(F, pro, self.convert_binning(F1,pro))
                except AttributeError:
                    if self.debug:
                        print("Cannot load attribute %s from %s" %(pro, 
                                                                   fiber_fn))
            if append_flag:
                self.fibers.append(F)
    
    
    def load_all_cal(self):
        '''
        Load all calibration properties to current fiber set and evaluate
        the fibermodel and wavelength solution.
        '''
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
        '''
        Save the fibers to self.path using the fiber class "save" function.
        '''
        for fiber in self.fibers:
            fiber.save(self.specid, self.ifuslot, self.ifuid, self.amp)
       
       
    def get_image(self):
        '''
        This many purpose function loads the image, finds the overscan value,
        subtracts the dark image and bias image, multplies the gain, 
        and divides the pixelflat.  It also orients the image using the
        "orient" function above.
        '''
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
        if self.use_pixelflat:
            pixelflat = np.array(fits.open(op.join(self.virusconfig, 
                                                   'PixelFlats','20161223',
                                                   'pixelflat_cam%s_%s.fits' 
                                            %(self.specid, self.amp)))[0].data,
                                  dtype=float)
        if self.dark_mult>0.0:
            image[:] = image - self.dark_mult * darkimage
        if self.bias_mult>0.0:
            image[:] = image - self.bias_mult * biasimage
        image[:] = image * self.gain
        if self.use_pixelflat:
            image[:] = np.where(pixelflat != 0, image / pixelflat, 0.0)
        self.image = self.orient(image)
        
    
    def find_shift(self):
        '''
        Find the shift in the trace compared to the calibration fibers in
        self.calpath.  If they are the same then the shift will be 0.  
        This has not been tested yet.
        '''
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
                        
       
    def get_trace(self):
        '''
        This function gets the trace for this amplifier.  It checks functional
        dependencies first: get_image().  If self.type is 'twi' or self.refit 
        is True then the trace is calculated, otherwise the trace is loaded 
        from calpath.
        '''                      
          
        if self.image is None:
            self.get_image()
        if self.type == 'twi' or self.refit:
            allfibers, xc = get_trace_from_image(self.image, interp_window=2.5,
                                                 debug=False)
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
                self.fibers[i].trace_poly_order = self.trace_poly_order
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
                                                  - yvals[floc]) < self.fdist):
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
                                                  - yvals[floc]) < self.fdist):
                        self.fibers[i].trace_x[c] = c
                        self.fibers[i].trace_y[c] = yvals[floc]
            for fib, fiber in enumerate(self.fibers):
                if np.sum(fiber.trace_x != fiber.flag)==len(xc):
                    fiber.fit_trace_poly()
                    fiber.eval_trace_poly()
            for fib, fiber in enumerate(self.fibers):
                if np.sum(fiber.trace_x != fiber.flag)!=len(xc):
                    sel = np.where(fiber.trace_x != fiber.flag)[0]
                    setdiff = np.setdiff1d(xc, sel)
                    fiber.fit_trace_poly()
                    k=1
                    done=False
                    if (fib+k)<=(len(self.fibers)-1):
                        if np.all(self.fibers[fib+k].trace > 0):
                            dif = np.interp(setdiff, fiber.trace_x[sel], 
                                            fiber.trace_y[sel] - 
                                            self.fibers[fib+k].trace[sel])
                            fiber.trace_x[setdiff] = setdiff
                            fiber.trace_y[setdiff] = (dif + 
                                             self.fibers[fib+k].trace[setdiff])
                            fiber.eval_trace_poly()
                            done=True
                    if (fib-k) >= 0 and not done:
                        if np.all(self.fibers[fib-k].trace > 0):
                            dif = np.interp(setdiff, fiber.trace_x[sel], 
                                            fiber.trace_y[sel] - 
                                            self.fibers[fib-k].trace[sel])
                            fiber.trace_x[setdiff] = setdiff
                            fiber.trace_y[setdiff] = (dif + 
                                             self.fibers[fib-k].trace[setdiff])
                            fiber.eval_trace_poly()
                            done=True
                    if not done:
                        fiber.trace_x[setdiff] = setdiff
                        fiber.trace_y[setdiff] = np.polyval(
                                                          fiber.trace_polyvals, 
                                                              setdiff / self.D)
                        fiber.eval_trace_poly()

                        
            if self.calculate_shift:
                self.net_trace_shift = self.find_shift()
                if self.net_trace_shift is not None:
                    for fiber in self.fibers:
                        fiber.trace_polyvals[-1] += self.net_trace_shift
                        fiber.eval_trace_poly()
                        fiber.trace[:] = self.trace + self.net_trace_shift
            
        else:
            self.load_cal_property('trace_polyvals')
            self.load_cal_property('trace')
            for fiber in self.fibers:
                fiber.eval_trace_poly()
                
        if self.check_trace:
            outfile = op.join(self.path,'trace_%s.png' %self.basename)
            check_fiber_trace(self.image, self.fibers, outfile)
            
                
    def get_fibermodel(self):
        '''
        This function gets the fibermodel for this amplifier.  It checks 
        functional dependencies first: get_image() and get_trace().  
        If self.type is 'twi' or self.refit is True then the fibermodel is 
        calculated, otherwise the fibermodel is loaded and evaluated
        from calpath.
        
        The fiber model is currently a linear interpolation between amplitudes
        of different bins over the fiber size.  Critically, "bins", "fsize",
        "sigma", and "power" set the initial x-values for each bin where x
        is the distance in pixels from the fiber trace.
        
        The fiber model is fit using image stamps containing a number of 
        fibers, "fiber_group", surrounding the target fiber 
        and number of columns, "col_group".  Over this
        stamp, a single fiber model is fit to all relevant fibers.  This is
        done for each fiber in a moving window sense, but done only once in
        the column direction in an gridded sense.
        
        '''
        if self.image is None:
            self.get_image()
        if not self.fibers:
            self.get_trace()
        if self.type == 'twi' or self.refit:
            sol, xcol, binx = fit_fibermodel_nonparametric(self.image, 
                                                                   self.fibers,
                                                                   debug=False,
                                         use_default=self.use_default_fibmodel,
                                                  plot=self.make_fib_ind_plots,
                                                           outfolder=self.path,
                                                  fiber_group=self.fiber_group,
                                                      col_group=self.col_group,
                                                      bins=self.fibmodel_nbins,
                                                              fsize=self.fsize,
                                                              sigma=self.sigma,
                                                              power=self.power)
            nfibs, ncols, nbins = sol.shape
            for i, fiber in enumerate(self.fibers):
                fiber.fibmodel_poly_order = self.fibmodel_poly_order
                fiber.fibmodel_x = xcol
                fiber.fibmodel_y = sol[i,:,:]
                fiber.binx = binx
                fiber.fit_fibmodel_poly()
                fiber.eval_fibmodel_poly()
        else:
            self.load_cal_property(['fibmodel_polyvals','binx'])
            for fiber in self.fibers:
                fiber.eval_fibmodel_poly()
                
        if self.check_fibermodel:
            if self.fibers[0].spectrum is None:
                norm = get_norm_nonparametric(self.image, self.fibers, 
                                              debug=False)
                for i, fiber in enumerate(self.fibers):
                    fiber.spectrum = norm[i,:]
            outfile = op.join(self.path,'fibmodel_%s.png' %self.basename)
            check_fiber_profile(self.image, self.fibers, outfile)


    def fiberextract(self):
        '''
        This function gets the spectrum for each fiber.  It checks 
        functional dependencies first: get_image(), get_trace(), and
        get_fibermodel(). 

        '''
        if self.image is None:
            self.get_image()
        if not self.fibers:
            self.get_trace()
        if self.fibers[0].fibmodel is None:
            self.get_fibermodel()
        else:
            for fiber in self.fibers:
                fiber.eval_fibmodel_poly()
        norm = get_norm_nonparametric(self.image, self.fibers, 
                                      debug=False, mask=self.mask)
        for i, fiber in enumerate(self.fibers):
            fiber.spectrum = norm[i,:]
    
    
    def get_wavelength_solution(self):
        '''
        This function gets the wavelength solution for each fiber.  It checks 
        functional dependencies first: get_image(), get_trace(),
        get_fibermodel(), and fiberextract().
        
        The wavelength solution is done for an initial fiber first.  In bins
        a linear solution is fit using a chi^2 minimization comparing a
        default solar/twi spectrum compared to this amplifiers twi spectrum.
        Both are normalized before chi^2 minimization.  A polynomial is then
        fit to the collection of linear solutions.  After the first fiber
        neighboring fibers start with previous solution as starting points
        for the binned linear fits.

        '''                                    
        solar_spec = np.loadtxt(op.join(self.virusconfig,
                                        'solar_spec/%s_temp.txt' 
                                        %self.specname))
        if self.image is None:
            self.get_image()
        if not self.fibers:
            self.get_trace()
        if self.fibers[0].fibmodel is None:
            self.get_fibermodel()
        else:
            for fiber in self.fibers:
                fiber.eval_fibmodel_poly() 
        if self.type == 'twi' or self.refit:
            if self.fibers[0].spectrum is None:
                norm = get_norm_nonparametric(self.image, self.fibers, 
                                              debug=False)
                for i, fiber in enumerate(self.fibers):
                    fiber.spectrum = norm[i,:]
            if self.init_lims is None:
                print("Please provide initial wavelength endpoint guess")
                sys.exit(1)
            for k in xrange(2):
                content=False
                cnting = 0
                while content==False:
                    fiber = self.fibers[self.default_fib]
                    fw, fwp = calculate_wavelength_chi2(np.arange(self.D), 
                                                        fiber.spectrum, 
                                                        solar_spec, 
                                                      init_lims=self.init_lims, 
                                                        debug=False, 
                                                  interactive=self.interactive,
                                                        nbins=self.wave_nbins)
                    fiber.wavelength = fw*1.
                    fiber.wave_polyvals = fwp*1.
                    # Boundary check:
                    if (np.abs(fiber.wavelength.min()-self.init_lims[0])>100. 
                     or np.abs(fiber.wavelength.max()-self.init_lims[1])>100.):
                        self.default_fib = np.random.choice(112) 
                        cnting+=1
                    else:
                        content=True
                    if cnting>9:
                        print("Ran through 9 loops but no solution found.")
                        sys.exit(1)
                fc = np.arange(len(self.fibers))
                if self.default_fib==0:
                    fibs1 = []
                else:
                    fibs1 = fc[(self.default_fib-1)::-1]
                if self.default_fib==(len(self.fibers)-1):
                    fibs2 = []
                else:
                    fibs2 = fc[(self.default_fib+1)::1]
                for fib in fibs1:
                    fiber = self.fibers[fib]
                    fw, fwp = calculate_wavelength_chi2(np.arange(self.D), 
                                                        fiber.spectrum, 
                                                        solar_spec, 
                                                      init_lims=self.init_lims, 
                                                        debug=False, 
                                                       interactive=False,
                                     init_sol=self.fibers[fib+1].wave_polyvals,
                                                        nbins=self.wave_nbins)
                    fiber.wavelength = fw*1.
                    fiber.wave_polyvals = fwp*1.
                for fib in fibs2:
                    fiber = self.fibers[fib]
                    fw, fwp = calculate_wavelength_chi2(np.arange(self.D), 
                                                        fiber.spectrum, 
                                                        solar_spec, 
                                                      init_lims=self.init_lims, 
                                                        debug=False, 
                                                       interactive=False,
                                     init_sol=self.fibers[fib-1].wave_polyvals,
                                                        nbins=self.wave_nbins)
                    fiber.wavelength = fw*1.
                    fiber.wave_polyvals = fwp*1.
                if k==0:
                    self.get_master_sky(norm=True)
                    solar_spec = np.zeros((len(self.masterwave),2))
                    solar_spec[:,0] = self.masterwave
                    solar_spec[:,1] = self.mastersmooth
        else:
            self.load_cal_property(['wave_polyvals'])
            for fiber in self.fibers:
                fiber.eval_wave_poly()
        if self.check_wave:
            if self.fibers[0].spectrum is None:
                norm = get_norm_nonparametric(self.image, self.fibers, 
                                                  debug=False)
                for i, fiber in enumerate(self.fibers):
                    fiber.spectrum = norm[i,:]
            outfile = op.join(self.path,'wavesolution_%s.png' %self.basename)
            check_wavelength_fit(self.fibers, solar_spec, outfile)
                
                    
    def get_fiber_to_fiber(self):
        '''
        This function gets the fiber to fiber normalization for this amplifier. 
        It checks functional dependencies first: get_image(), get_trace(), 
        get_fibermodel(), fiberextract(), and get_wavelength_solution().
        '''
        if self.image is None:
            if self.debug:
                print("Building image for %s" %self.basename)
            self.get_image()
        if not self.fibers:
            if self.debug:
                print("Getting trace for %s" %self.basename)
            self.get_trace()
        if self.fibers[0].fibmodel is None:
            if self.debug:
                print("Getting fibermodel for %s" %self.basename)
            self.get_fibermodel()
        else:
            for fiber in self.fibers:
                fiber.eval_fibmodel_poly()
        if self.type == 'twi' or self.refit:
            if self.fibers[0].spectrum is None:
                if self.debug:
                    print("Fiberextracting %s" %self.basename)
                self.fiberextract()
            if self.fibers[0].wave_polyvals is None:
                if self.debug:
                    print("Getting Wavelength solution for %s" %self.basename)
                self.get_wavelength_solution()
            else:
                for fiber in self.fibers:
                    fiber.eval_wave_poly()
            if self.debug:
                print("Getting Fiber to Fiber for %s" %self.basename)
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
            self.averagespec = biweight_filter(smoothspec, self.filt_size_agg)
            for fib, fiber in enumerate(self.fibers):
                fiber.fiber_to_fiber = biweight_filter(masterspec[fib] 
                                      / np.interp(fiber.wavelength, masterwave, 
                                                  self.averagespec), 
                                                       self.filt_size_final)

        else:
            self.load_cal_property(['fiber_to_fiber'])   
            
        
    def sky_subtraction(self):
        '''
        This function gets the master sky spectrum and 
        evaluates the sky_spectrum for each fiber. It then builds a sky image
        and a sky-subtracted image.  It checks functional dependencies first: 
        get_image(), get_trace(), get_fibermodel(), fiberextract(), 
        get_wavelength_solution(), and fiber_to_fiber().        
        '''
        if self.image is None:
            self.get_image()
        if not self.fibers:
            self.get_trace()
        if self.fibers[0].fibmodel is None:
            self.get_fibermodel()
        else:
            for fiber in self.fibers:
                fiber.eval_fibmodel_poly()
        if self.fibers[0].fiber_to_fiber is None:
            self.load_cal_property(['fiber_to_fiber'])
        if self.fibers[0].spectrum is None:
            self.fiberextract()
        if self.fibers[0].wave_polyvals is None:
            self.get_wavelength_solution()
        else:
            for fiber in self.fibers:
                fiber.eval_wave_poly()
        if self.fibers[0].fiber_to_fiber is None:
            self.get_fiber_to_fiber()
        if self.skypath is not None:
            self.load_cal_property('sky_spectrum', pathkind='skypath')
            if self.fibers[0].sky_spectrum is None:
                print("Loading sky spectrum from %s did not work." 
                      %self.skypath)
                print("Setting skypath to None for this amplifier.")
                self.skypath = None
        if self.skypath is None:
            self.get_master_sky(sky=True)
            for fib, fiber in enumerate(self.fibers):
                fiber.sky_spectrum = (fiber.fiber_to_fiber 
                                 * np.interp(fiber.wavelength, self.masterwave, 
                                             self.mastersky))
         
        self.skyframe = get_model_image(self.image, self.fibers, 
                                        'sky_spectrum', debug=False)
        self.clean_image = self.image - self.skyframe
        
    def clean_cosmics(self):
         cc = cosmics.cosmicsimage(self.clean_image, gain=1.0, 
                                   readnoise=self.rdnoise, 
                                   sigclip=25.0, sigfrac=0.001, objlim=0.001,
                                   satlevel=-1.0)
         cc.run(maxiter=4)
         c = np.where(cc.mask == True)
         self.mask = np.zeros(self.image.shape)
         for x, y in zip(c[0], c[1]):
             self.mask[x][y] = -1.0 
             
    def get_master_sky(self, sky=False, norm=False):
        masterwave = []
        if sky:
            masterspec = []
        if norm:
            mastersmooth=[]
        for fib, fiber in enumerate(self.fibers):
            if norm:
                y = biweight_filter(fiber.spectrum, self.filt_size_ind)
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
            self.mastersky = biweight_filter(masterspec, self.filt_size_sky)
        if norm:
            mastersmooth = np.hstack(mastersmooth)
            mastersmooth[:] = mastersmooth[ind]
            self.mastersmooth = biweight_filter(mastersmooth, 
                                                self.filt_size_sky)     