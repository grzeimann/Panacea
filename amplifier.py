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
from utils import biweight_location, biweight_filter, biweight_midvariance
from utils import biweight_bin, is_outlier
from astropy.io import fits
import os.path as op
import numpy as np
import sys
import re
import glob
from fiber_utils import get_trace_from_image, get_indices, get_interp_list
from fiber_utils import check_fiber_trace, measure_background
from fiber_utils import calculate_wavelength_chi2, get_model_image
from fiber_utils import check_fiber_profile, check_wavelength_fit
from fiber_utils import fast_nonparametric_fibermodel, new_fast_norm
from fiber_utils import calculate_significance, get_wavelength_offsets
from fiber_utils import bspline_x0
from fiber import Fiber
from scipy.signal import medfilt
from scipy.interpolate import splrep, splev
from scipy.ndimage.filters import percentile_filter
import cosmics
from datetime import datetime
import logging
import matplotlib.pyplot as plt

__all__ = ["Amplifier"]

class Amplifier:
    def __init__(self, filename, path, name=None, refit=False, calpath=None, 
                 skypath=None, verbose=True, darkpath=None, biaspath=None,
                 pixflatpath=None,
                 virusconfig=None, dark_mult=0., bias_mult=0., 
                 use_pixelflat=True, specname=None, fdist=2., fdist_ref=4., 
                 check_trace=True, adjust_trace=False, trace_poly_order=3,
                 fibmodel_poly_order=3, use_default_fibmodel=False, 
                 fibmodel_nbins=31, make_fib_ind_plots=False, 
                 check_fibermodel=False, fsize=8., sigma=2.5, power=2.5,
                 fiber_group=8, col_group=48, mask=None, wave_nbins=21, 
                 wave_order=3, default_fib=0, init_lims=None, collapse_lims=None,
                 interactive=False, check_wave=False, filt_size_ind=51, 
                 filt_size_agg=51, filt_size_final=151, filt_size_sky=151,
                 col_frac = 0.47, use_trace_ref=False, fiber_date=None,
                 cont_smooth=25, make_residual=True, do_cont_sub=True,
                 make_skyframe=True, wave_res=1.9, trace_step=4,
                 fibmodel_slope=0.000, fibmodel_intercept=0.000,
                 fibmodel_breakpoint=5.5, fibmodel_step=4,
                 fibmodel_interpkind='linear', cosmic_iterations=1,
                 sky_scale=1.0, make_model_image=False, init_sol=None,
                 wavestepsize=1., nknots=51, bspline_binsize=200.,
                 bspline_waveres=1.0, sky_iterations=3,
                 sky_sigthresh=2.5, adjust_ftf=False, trace_y_window=3.,
                 trace_repeat_length=2):
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
        :param fdist_ref:
            This parameter sets the distance allowable to match fiber
            centroid measurments for reference fibers.
        :param check_trace:
            Plot trace images at the end.  Good for post checks.
        :param adjust_trace:
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
        :param col_frac:
            Column that is this fraction through the image.
            In other words, col = int(num_cols * col_frac)
        :param use_trace_ref:
            This is a boolean variable that sets whether the fibers found
            are compared with reference fibers and define dead fibers if
            there are missing fibers in the search.
        :param fiber_date:
            The date on which the fiber reference is located.
            
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
        self.verbose = verbose
        self.setup_logging()

        # Organization options
        try:
            F = fits.open(filename)
        except:
            self.log.error("%s could not be opened with Amplifier class" %filename)
            return None
        self.header = F[0].header
        try:
            self.header.remove('bzero')
        except:
            self.log.warning('No BZERO to remove')
        self.name = name
        self.filename = op.abspath(F.filename())
        self.basename = op.basename(F.filename())[:-5]
        self.virusconfig = virusconfig
        self.path = path
        if not op.exists(self.path):
            mkpath(self.path)
        self.refit = refit
        self.calpath = calpath
        self.skypath = skypath
        self.darkpath = darkpath
        self.biaspath = biaspath
        self.pixflatpath = pixflatpath
        self.specname = specname
        
        # Image manipulation options
        self.dark_mult = dark_mult
        self.bias_mult = bias_mult
        self.use_pixelflat = use_pixelflat
        
        # Trace options
        self.fdist = fdist
        self.check_trace = check_trace
        self.trace_poly_order = trace_poly_order
        self.adjust_trace = adjust_trace
        self.col_frac = col_frac
        self.use_trace_ref = use_trace_ref
        self.fdist_ref = fdist_ref
        self.fiber_date = fiber_date
        self.trace_step = trace_step
        self.trace_y_window = trace_y_window
        self.trace_repeat_length = trace_repeat_length
        
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
        self.fibmodel_slope = fibmodel_slope
        self.fibmodel_intercept = fibmodel_intercept
        self.fibmodel_breakpoint = fibmodel_breakpoint
        self.fibmodel_step = fibmodel_step
        self.fibmodel_interpkind = fibmodel_interpkind
        
        # Masking options (Fiberextract related)
        self.mask = mask
        self.cosmic_iterations = cosmic_iterations
        
        # Wavelength Solution options
        self.wave_nbins = wave_nbins
        self.wave_order = wave_order
        self.default_fib = default_fib
        self.init_lims = init_lims
        self.collapse_lims = collapse_lims
        self.interactive = interactive
        self.check_wave = check_wave        
        self.wave_res = wave_res    
        self.init_sol = init_sol
        
        # Smoothing options for individual and master spectra
        self.filt_size_ind = filt_size_ind
        self.filt_size_agg = filt_size_agg
        self.filt_size_final = filt_size_final
        self.filt_size_sky = filt_size_sky
        self.sky_scale = sky_scale
        self.wavestepsize = wavestepsize
        # Continuum subtraction
        self.cont_smooth = cont_smooth
        
        # Sky / Fiber to Fiber
        self.nknots = nknots
        self.bspline_binsize = bspline_binsize
        self.bspline_waveres = bspline_waveres
        self.sky_iterations = sky_iterations
        self.sky_sigthresh = sky_sigthresh
        self.adjust_ftf = adjust_ftf
        self.back = 0.0
        
        # Image Options
        self.make_residual = make_residual
        self.do_cont_sub = do_cont_sub
        self.make_skyframe = make_skyframe
        self.make_model_image = make_model_image
        
        # Initialized variables
        self.N, self.D = F[0].data.shape
        if self.D == 1064:
            self.D -= 32
        if self.D == 2128:
            self.D -= 64
        if self.N == 1074:
            self.N -= 50
            self.D -= 45
        if self.N == 2048:
            self.D -= 44
        self.trimmed = False
        self.overscan_value = None
        try:
            self.gain = F[0].header['GAIN']
        except:
            self.gain = 1.
        if self.gain == 0.0:
            self.gain = 1.
        new_keywords = ['TRAJCRA', 'TRAJCDEC', 'PARANGLE',
                        'RHO_STRT']
        att_list = ['ra', 'dec', 'pa', 'rho']
        for newk, att in zip(new_keywords, att_list):
            try:
                setattr(self, att, F[0].header[newk])
            except:
                setattr(self, att, None)
        try:
            self.rdnoise = F[0].header['RDNOISE']
        except:
            self.rdnoise = 3.
        self.amp = (F[0].header['CCDPOS'].replace(' ', '') 
                    + F[0].header['CCDHALF'].replace(' ', ''))
        trim = re.split('[\[ \] \: \,]', F[0].header['TRIMSEC'])[1:-1]
        self.trimsec = [int(t)-((i+1)%2) for i,t in enumerate(trim)]        
        bias = re.split('[\[ \] \: \,]', F[0].header['BIASSEC'])[1:-1]
        self.biassec = [int(t)-((i+1)%2) for i,t in enumerate(bias)]        
        self.fibers = []
        self.type = F[0].header['IMAGETYP'].replace(' ', '')
        self.specid = '%03d' %F[0].header['SPECID']
        self.ifuid = '%03d' %(int(F[0].header['IFUID'].replace(' ', '')))
        self.ifuslot ='%03d' %F[0].header['IFUSLOT']
        datetemp = re.split('[-,T]',F[0].header['DATE-OBS'])
        self.date = datetime(int(datetemp[0]), int(datetemp[1]), 
                             int(datetemp[2]))
        self.image = np.array(F[0].data, dtype=float)
        if self.rdnoise<0.1:
            self.rdnoise = 2.5
        self.image_prepped = False
        if self.gain>0:
            self.error = self.rdnoise / self.gain * np.ones((self.N,self.D), 
                                                            dtype=float)
        else:
            self.error = np.zeros((self.N, self.D), dtype=float)
        
        self.exptime = F[0].header['EXPTIME']
        self.ifupos = None
        try:
            self.ampname = F[0].header['AMPNAME']
        except:
            self.ampname = None
        if self.ampname is None:
            try:
                self.ampname = F[0].header['AMPLIFIE']
            except:
                self.ampname = None            
            
    def setup_logging(self):
        '''Set up a logger for shuffle with a name ``amplifier``.
    
        Use a StreamHandler to write to stdout and set the level to DEBUG if
        verbose is set from the command line
        '''
        self.log = logging.getLogger('panacea')
        if not len(self.log.handlers):
            fmt = '[%(levelname)s - %(asctime)s] %(message)s'
            if not self.verbose:
                level = logging.WARNING
            else:
                level = logging.INFO
           
            fmt = logging.Formatter(fmt)
        
            handler = logging.StreamHandler()
            handler.setFormatter(fmt)
            handler.setLevel(level)
        
            self.log = logging.getLogger('panacea')
            self.log.setLevel(logging.DEBUG)
            self.log.addHandler(handler)

        
    def write_to_fits(self, hdu, outname):
        '''
        Writing fits file to outname
        '''
        mkpath(op.dirname(outname))
        try:
            hdu.writeto(outname, overwrite=True)
        except TypeError:
            hdu.writeto(outname, clobber=True)

    def write_header(self, hdu):
        hdu.header = self.header
        hdu.header['RA'] = self.ra
        hdu.header['DEC'] = self.dec
        hdu.header['PA'] = self.pa
        hdu.header['RHO'] = self.rho
        hdu.header['EXPTIME'] = self.exptime
        hdu.header['OSCANMN'] = self.overscan_value
        hdu.header['BACKLVL'] = self.back
        hdu.header['OSCANSTD'] = self.overscan_noise * self.gain
        hdu.header['RAWPATH'] = self.path
        hdu.header['RAWFN'] = self.filename
        hdu.header['FSIZE'] = self.fsize
        hdu.header['IFUSLOT'] = self.ifuslot
        hdu.header['SPECID'] = self.specid
        hdu.header['IFUSID'] = self.ifuid
        hdu.header['AMP'] = self.amp
        if self.biaspath is not None:
            hdu.header['BIAPATH'] = self.biaspath
        return hdu
            
    def save(self, image_list=[], spec_list=[]):
        '''
        Save the entire amplifier include the list of fibers.  
        This property is not used often as "amp*.pkl" is large and typically
        the fibers can be loaded and the other amplifier properties quickly
        recalculated.
        '''
        self.log.info('Saving images/properties to %s' %self.path)
        fn = op.join(self.path, 'multi_%s_%s_%s_%s.fits' %(self.specid, 
                                                           self.ifuslot,
                                                           self.ifuid,
                                                           self.amp))
        fits_list = []
        for i,image in enumerate(image_list):
            if i==0:
                fits_list.append(fits.PrimaryHDU(np.array(getattr(self, image), dtype='float32')))
            else:
                fits_list.append(fits.ImageHDU(np.array(getattr(self, image), dtype='float32')))

            fits_list[-1].header['EXTNAME'] = image
            
        for i, spec in enumerate(spec_list):
            try:
                s = np.array([getattr(fiber, spec) for fiber in self.fibers], dtype='float32')
                fits_list.append(fits.ImageHDU(s))
                fits_list[-1].header['EXTNAME'] = spec
            except AttributeError:
                self.log.warning('Attribute %s does not exist to save' %spec)
        if fits_list:
            fits_list[0] = self.write_header(fits_list[0])
            hdu = fits.HDUList(fits_list)
            self.write_to_fits(hdu, fn)
           
    def load(self, path='path', image_list=[], spec_list=[]):
        '''
        Save the entire amplifier include the list of fibers.  
        This property is not used often as "amp*.pkl" is large and typically
        the fibers can be loaded and the other amplifier properties quickly
        recalculated.
        '''
        self.log.info('Loading images/properties from %s' %getattr(self, path))
        fn = op.join(getattr(self, path), 'multi_%s_%s_%s_%s.fits' %(self.specid, 
                                                           self.ifuslot,
                                                           self.ifuid,
                                                           self.amp))
        try:
            F = fits.open(fn)
        except IOError:
            self.log.error('Failed to open %s' %fn)
            return None
        
        for i,image in enumerate(image_list):
            try:
                setattr(self, image, F[image].data)
            except KeyError:
                self.log.error('Failed to open extension %s for %s' %(image, fn))
                
        for i, spec in enumerate(spec_list):
            try:
                a = F[spec].data.shape[0]
            except KeyError:
                self.log.error('Failed to open extension %s for %s' %(spec, fn))
                return None
            for j in np.arange(a):
                try:
                    f = self.fibers[j]
                except IndexError:    
                    f = Fiber(self.D, j+1, self.path, self.filename)                
                    self.fibers.append(f)
                if spec=='spectrum' and path=='calpath':
                    setattr(self.fibers[j], 'twi_'+spec, F[spec].data[j])
                else:
                    setattr(self.fibers[j], spec, F[spec].data[j])
            if spec=='trace':
                get_indices(self.image, self.fibers, self.fsize)
            if spec=='dead':
                self.good_fibers = [fiber for fiber in self.fibers 
                                      if not fiber.dead]
                self.dead_fibers = [fiber for fiber in self.fibers 
                                      if fiber.dead]
            


    def convert_binning(self, values, b):
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
        
        return np.interp((1.*np.arange(self.D))/self.D, 
                         (1.*np.arange(b))/b, values)
            

    def save_fibmodel(self):
        '''
        Save the fibers to fits file with two extensions
        The first is 3-d for trace, wavelength and fiber_to_fiber
        The second is which fibers are good and not dead
        '''
        try: 
            self.fibers[0].fibmodel
        except:
            self.log.warning('Trying to save fibermodel but none exist.')
            return None
        ylims = np.linspace(-1.*self.fsize, self.fsize, self.fibmodel_nbins)
        fibmodel = np.zeros((len(self.fibers), self.fibmodel_nbins, self.D))
        for i,fiber in enumerate(self.fibers):
            fibmodel[i,:,:] = fiber.fibmodel
        s = fits.PrimaryHDU(np.array(fibmodel,dtype='float32'))
        t = fits.ImageHDU(np.array(ylims,dtype='float32'))
        hdu = fits.HDUList([s,t])
        fn = op.join(self.path, 'fibermodel_%s_%s_%s_%s.fits' %(self.specid, 
                                                               self.ifuslot,
                                                               self.ifuid,
                                                               self.amp))
        self.write_to_fits(hdu, fn)
            
    def load_fibmodel(self, path='path'):
        '''
        Load fibers in self.path. Redefine the path for the fiber to self.path
        in case it was copied over. After loaded, evaluate the fibermodel as 
        well as the wavelength solution if available.  
        '''
        self.log.info('Loading fiber models from %s' %getattr(self, path))
        fn = op.join(getattr(self, path), 'fibermodel_%s_%s_%s_%s.fits' %(self.specid, 
                                                            self.ifuslot,
                                                            self.ifuid,
                                                            self.amp))
        try:
            F = fits.open(fn)
        except IOError:
            self.log.warning('Error opening fibers from %s' %fn)
            self.log.warning('%s does not exist.' %fn)
            return None
        a,b,c = F[0].data.shape
        ygrid, xgrid = np.indices(self.image.shape)
        for i in xrange(a):
            try:
                self.fibers[i]
            except IndexError:    
                self.log.warning('Fiber %i does not have the necessary trace info.'  %i)
                return None
            self.fibers[i].yoff = self.fibers[i].yind - self.fibers[i].trace[self.fibers[i].xind]
            interp_list = get_interp_list(self.fsize, self.fibmodel_nbins, 
                                          self.fibmodel_interpkind)
            self.fibers[i].core = np.zeros((len(self.fibers[i].xind),))
            for j,iv in enumerate(interp_list):
                self.fibers[i].core += (iv(self.fibers[i].yoff)
                                        *F[0].data[i, j, self.fibers[i].xind])

      
    def orient_image(self):
        '''
        Orient the images from blue to red (left to right)
        Fibers are oriented to match configuration files
        '''
        self.log.info('Orienting Image for %s' %self.basename)
        if self.amp == "LU":
            self.image[:] = self.image[::-1,::-1]
            self.error[:] = self.error[::-1,::-1]
            if self.mask is not None:
                self.mask[:] = self.mask[::-1,::-1]
        if self.amp == "RL":
            self.image[:] = self.image[::-1,::-1]
            self.error[:] = self.error[::-1,::-1]
            if self.mask is not None:
                self.mask[:] = self.mask[::-1,::-1]
        if self.ampname is not None:
            if self.ampname == 'LR' or self.ampname == 'UL':
                self.image[:] = self.image[:,::-1]
                self.error[:] = self.error[:,::-1]
                if self.mask is not None:
                    self.mask[:] = self.mask[:,::-1]
        
    def subtract_overscan(self):
        '''
        Evaluate the overscan_value using a biweight average of the BIASSEC
        region.  Only calculate the value if one does not exist or 
        recalculate is set to True.
        '''
        if self.overscan_value is None:
            ly = self.biassec[2]
            hy = self.biassec[3]
            lx = self.biassec[0]+1
            hx = self.biassec[1]
            self.overscan_value = biweight_location(self.image[ly:hy, lx:hx])
            self.overscan_col = biweight_location(self.image[ly:hy, lx:hx],
                                                  axis=(1,))
            self.overscan_noise = biweight_midvariance(self.image[ly:hy,lx:hx])
            self.image[:] = self.image - self.overscan_col[:,np.newaxis]
            self.log.info('Subtracting overscan value %0.3f from %s' 
                           %(self.overscan_value, self.basename))



    def trim_image(self):
        '''
        Trim the image to just the datasec/trimsec.
        '''
        self.log.info('Trimming image %s' %self.basename)
        if not self.trimmed:
            self.image = self.image[self.trimsec[2]:self.trimsec[3], 
                                       self.trimsec[0]:self.trimsec[1]]
            self.trimmed = True
      
    def get_closest_date(self, path, name):
        fn = sorted(glob.glob(op.join(path, '*', name)))
        dates = [op.basename(op.dirname(f)) for f in fn]
        timediff = np.zeros((len(dates),))
        for i, date in enumerate(dates):
            d = datetime(int(date[:4]), int(date[4:6]),
                         int(date[6:]))
            timediff[i] = (self.date - d).days
        sel = np.where(timediff>=0)[0]
        ind = sel[np.argmin(timediff[sel])]
        return op.join(path, dates[ind], name)
     
    def subtract_dark(self):
        if self.dark_mult>0.0:
            self.log.info('Subtracting dark with multiplication %0.2f for %s' 
                      %(self.dark_mult, self.basename))
            darkpath = self.get_closest_date(self.darkpath,
                                             'masterdark_%s_%s.fits' %
                                             (self.specid, self.amp))
            darkimage = np.array(fits.open(darkpath)[0].data, 
                                 dtype=float)
            self.image[:] = self.image - self.dark_mult * darkimage
            # error in ADU calculation
            self.error[:] = np.sqrt(self.error**2 
                                    + 1./self.gain*self.dark_mult*darkimage)
            
            
    def subtract_bias(self):
        if self.bias_mult>0.0:
            self.log.info('Subtracting bias frame from folder: %s' % self.biaspath)
            self.log.info('Subtracting bias with multiplication %0.2f for %s' 
                      %(self.bias_mult, self.basename))
            biaspath = self.get_closest_date(self.biaspath,
                                             'masterbias_%s_%s.fits' %
                                             (self.specid, self.amp))
            biasimage = np.array(fits.open(biaspath)[0].data, 
                                 dtype=float)
            self.image[:] = self.image - self.bias_mult * biasimage
            #self.error[:] = np.sqrt(self.error**2 + self.gain*self.bias_mult*biasimage)
            
            
    def multiply_gain(self):
        self.log.info('Multiplying gain for %s' % self.basename)
        self.image *= self.gain
        self.error *= self.gain
        
        
    def calculate_photonnoise(self):
        self.log.info('Calculating photon noise for %s' % self.basename)
        self.error[:] = np.sqrt(self.error**2 
                                + np.where(self.image>1e-8,self.image, 0.))
        
        
    def divide_pixelflat(self):
        if self.use_pixelflat:
            self.log.info('Dividing pixelflat for %s' % self.basename)
            pixelpath = self.get_closest_date(self.pixflatpath,
                                              'pixelflat_cam%s_%s.fits' %
                                              (self.specid, self.amp))
            pixelflat = np.array(fits.open(pixelpath)[0].data,
                                 dtype=float)
            self.image[:] = np.where(pixelflat > 1e-1, self.image / pixelflat, 
                                     0.0)
            self.error[pixelflat<=0.0] = -1.0
            self.mask = np.where(pixelflat > 1e-1, 0.0, -1.0)
            sel = np.isfinite(self.image)
            if (~sel).sum()>0:
                self.image[~sel] = 0.0
                self.error[~sel] = -1.0

    def subtract_background(self):
        if not self.image_prepped:
            self.prepare_image()
        if not self.fibers:
            self.get_trace()
        self.log.info('Subtracting background using "blank" pixels for %s' 
                      %self.basename)
        if self.fibers[0].xind is None:
            get_indices(self.image, self.fibers, self.fsize)
        self.back = measure_background(self.image, self.fibers)
        self.log.info('Background level is: %0.2f e-' % self.back) 

        self.image[:] = self.image - self.back
        
    def prepare_image(self):
        '''
        This many purpose function loads the image, finds the overscan value,
        subtracts the dark image and bias image, multplies the gain, 
        and divides the pixelflat.  It also orients the image using the
        "orient" function above.
        '''
        self.subtract_overscan()
        self.trim_image()
        self.subtract_bias()
        self.subtract_dark()
        self.multiply_gain()
        self.calculate_photonnoise()        
        self.divide_pixelflat()
        self.orient_image()
        self.image_prepped = True
        
    
    def find_shift(self):
        '''
        Find the shift in the trace compared to the calibration fibers in
        self.calpath.  If they are the same then the shift will be 0.  
        This has not been tested yet.
        '''
        self.log.info('Finding trace shift from fibers in %s' 
                      %getattr(self, 'calpath'))
        if getattr(self, 'calpath') is None:
            return None
        fn = op.join(getattr(self, 'calpath'), 'multi_%s_%s_%s_%s.fits' %(self.specid, 
                                                                self.ifuslot,
                                                                self.ifuid,
                                                                self.amp))
        try:
            F = fits.open(fn)
        except IOError:
            self.log.warning('Error opening fibers from %s' %fn)
            self.log.warning('%s does not exist.' %fn)
            return None
        
        try:
            trace_cal = F['trace'].data[:,:]
        except KeyError:
            self.log.warning('Error opening trace from %s' %fn)
            return None       
        if trace_cal.shape[0] != len(self.fibers):
            self.log.error('Mismatch in number of fibers from cal and sci.')
            self.log.error('Cal Fibers: %i, Sci Fibers: %i' %(trace_cal.shape[0],
                                                              len(self.fibers)))
            self.log.error('BEYOND THIS POINT, THE AMP REDUCTIONS WILL FAIL.')
            
        shift = []
        for i,fiber in enumerate(self.fibers):
            col = 3.*self.D/5.
            width = 20
            low = int(col-width)
            high = int(col+width+1)
            shift.append(biweight_location(fiber.trace[low:high] 
                              - trace_cal[i,low:high]))
            fiber.trace = 1.*trace_cal[i,:]
            
        self.shift = shift
        smooth_shift = biweight_filter(self.shift, 7, 1)
        self.log.info("Shift for %s is %0.3f" %(self.amp, biweight_location(np.array(shift))))
        for i,fiber in enumerate(self.fibers):
            fiber.trace = fiber.trace + smooth_shift[i]
            self.log.info("Shift for fiber %i is %0.3f pixels" %(i+1, smooth_shift[i]))

        
        return biweight_location(np.array(shift))


    def get_attribute_from_neighbor_fiber(self, ind1, ind2, prop):
        '''
        This small function is used to get an attribute from the fiber below
        it.
        '''
        if isinstance(prop, basestring):
            prop = [prop]
        for i, pro in enumerate(prop):
            setattr(self.fibers[ind1], pro, 1.*getattr(self.fibers[ind2], pro))

    
    def fill_in_dead_fibers(self, props, return_ind=False):
        '''
        This function fills in dead fibers by looking at neighbors which are
        not dead and filling in their info.
        '''
        if return_ind:
            keep_ind=[]
        for fiber in self.dead_fibers:
            i = fiber.fibnum - 1
            j = 1
            no_neighbor = True
            while no_neighbor:
                if i-j >=0:
                    if not self.fibers[i-j].dead:
                        self.get_attribute_from_neighbor_fiber(i, i-j, props)
                        no_neighbor=False
                        if return_ind:
                            keep_ind.append(i-j)
                if i+j<len(self.fibers) and no_neighbor:
                    if not self.fibers[i+j].dead:
                        self.get_attribute_from_neighbor_fiber(i, i+j, props)
                        no_neighbor=False
                        if return_ind:
                            keep_ind.append(i+j)
                j+=1
                if i-j<0 and ((i+j)>=len(self.fibers)):
                    print(i,j,len(self.dead_fibers))
                    sys.exit(1)
        if return_ind:
            return keep_ind
            


    def get_trace(self):
        '''
        This function gets the trace for this amplifier.  It checks functional
        dependencies first: prepare_image().  If self.type is 'twi' or self.refit 
        is True then the trace is calculated, otherwise the trace is loaded 
        from calpath.
        '''                      
          
        if not self.image_prepped:
            self.prepare_image()
        
        if self.type == 'twi' or self.refit:
            self.log.info('Measuring the trace from %s' %self.basename)
            allfibers, xc = get_trace_from_image(self.image, interp_window=2.5,
                                                 debug=False, mx_cut=0.1,
                                                 x_step=self.trace_step,
                                                 y_window=self.trace_y_window,
                                                 repeat_length=self.trace_repeat_length)
            brcol = np.argmin(np.abs(xc-self.D*self.col_frac))
            if self.use_trace_ref:
                if self.fiber_date is None:
                    fn = glob.glob(op.join(self.virusconfig, 
                                           'Fiber_Locations','*')+'/')
                    dates = [op.basename(op.dirname(f)) for f in fn]
                    timediff = np.zeros((len(dates),))
                    for i,date in enumerate(dates):
                        d = datetime(int(date[:4]), int(date[4:6]),
                                     int(date[6:]))
                        timediff[i] = np.abs((self.date - d).days)
                    timesel = np.argsort(timediff)
                    closest_date = dates[np.argmin(timediff)]
                    loaded = False
                    k = 0
                    while not loaded:
                        closest_date = dates[timesel[k]]
                        try:
                            ref_file = np.loadtxt(op.join(self.virusconfig, 
                                                  'Fiber_Locations',
                                                  closest_date, 
                                                  'fiber_loc_%s_%s_%s_%s.txt'
                                                  %(self.specid, self.ifuslot, 
                                                    self.ifuid, self.amp)))
                            self.log.info('Loading trace reference from %s' %closest_date)
                            loaded = True
                        except:
                            k = k+1
                            if k >= len(timesel):
                                loaded=False
                                self.log.error('Could not find any trace refence file for: %s %s %s %s' 
                                               %(self.specid, self.ifuslot, 
                                                    self.ifuid, self.amp))
                                return None
                else:    
                    ref_file = np.loadtxt(op.join(self.virusconfig, 
                                                  'Fiber_Locations',
                                                  self.fiber_date, 
                                                  'fiber_loc_%s_%s_%s_%s.txt'
                                                  %(self.specid, self.ifuslot, 
                                                    self.ifuid, self.amp)))
                # sort by y values just in case the missing fiber is at the end
                ref_file = ref_file[ref_file[:,0].argsort(),:]
                standardcol=ref_file[:,0]
            else:
                standardcol = allfibers[brcol]
            cols1 = xc[brcol::-1]
            cols2 = xc[(brcol+1)::1]
            # Initialize fiber traces
            for i in xrange(len(standardcol)):
                try: 
                    F = self.fibers[i]
                except IndexError:    
                    F = Fiber(self.D, i+1, self.path, self.filename)
                    self.fibers.append(F)
                self.fibers[i].trace_poly_order = self.trace_poly_order
                self.fibers[i].init_trace_info()
                fdist = standardcol[i]-allfibers[brcol]
                floc = np.argmin(np.abs(fdist))
                if np.abs(fdist[floc])<self.fdist_ref:    
                    self.fibers[i].trace_x[int(xc[brcol])] = xc[brcol]
                    self.fibers[i].trace_y[int(xc[brcol])] = allfibers[brcol][floc]
                else:
                    self.fibers[i].dead=True
            self.good_fibers = [fiber for fiber in self.fibers 
                                      if not fiber.dead]
            self.dead_fibers = [fiber for fiber in self.fibers 
                                      if fiber.dead]  
            for c in cols1:
                loc = np.where(xc==c)[0][0]
                for i,fiber in enumerate(self.good_fibers):
                    yvals = allfibers[int(loc)]
                    if not yvals.size:
                        continue
                    xloc = np.argmin(np.abs(fiber.trace_x - c))
                    floc = np.argmin(np.abs(fiber.trace_y[xloc] 
                                            - yvals))
                    if (np.abs(fiber.trace_y[xloc]-yvals[floc]) < self.fdist):
                        fiber.trace_x[c] = c
                        fiber.trace_y[c] = yvals[floc]
            for c in cols2:
                loc = np.where(xc==c)[0][0]
                for i,fiber in enumerate(self.good_fibers):
                    yvals = allfibers[int(loc)]
                    if not yvals.size:
                        continue
                    xloc = np.argmin(np.abs(fiber.trace_x - c))
                    floc = np.argmin(np.abs(fiber.trace_y[xloc] 
                                            - yvals))
                    if (np.abs(fiber.trace_y[xloc]-yvals[floc]) < self.fdist):
                        fiber.trace_x[c] = c
                        fiber.trace_y[c] = yvals[floc]
            # Evaluate good fibers with trace defined everywhere measured
            for fib, fiber in enumerate(self.good_fibers):
                fiber.fit_trace_poly()
                fiber.eval_trace_poly()

            ind = self.fill_in_dead_fibers(['trace_x', 'trace_y', 'trace'], 
                                           return_ind=True)
            for i, fiber in enumerate(self.dead_fibers):
                fiber.trace_y += (standardcol[fiber.fibnum-1]
                                  -standardcol[ind[i]])
                fiber.trace += (standardcol[fiber.fibnum-1]
                                -standardcol[ind[i]])
            fn = op.join(self.path, self.basename + '_trace.txt')
            A = np.zeros((len(self.fibers),4))
            A[:,0] = [fiber.fibnum for fiber in self.fibers]
            A[:,1] = [fiber.trace[int(self.D/6.)] for fiber in self.fibers]
            A[:,2] = [fiber.trace[int(self.D/2.)] for fiber in self.fibers]
            A[:,3] = [fiber.trace[int(5.*self.D/6.)] for fiber in self.fibers]
            if self.check_trace:
                np.savetxt(fn, A)
            if self.adjust_trace and self.refit:
                self.net_trace_shift = self.find_shift()
            self.log.info('Trace measured from %s' %self.basename)
            get_indices(self.image, self.fibers, self.fsize)
        else:
            self.load(path='calpath', spec_list=['trace','dead'])       
        
        if self.check_trace:
            outfile = op.join(self.path,'trace_%s.png' %self.basename)
            check_fiber_trace(self.image, self.fibers, outfile)
            
                
    def get_fibermodel(self):
        '''
        This function gets the fibermodel for this amplifier.  It checks 
        functional dependencies first: prepare_image() and get_trace().  
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
        if not self.image_prepped:
            self.prepare_image()
        if not self.fibers:
            self.get_trace()
        if self.type == 'twi' or self.refit:
            self.log.info('Measuring the fiber model from %s' %self.basename)
            fast_nonparametric_fibermodel(self.image, self.good_fibers, 
                                          self.fsize, self.fibmodel_nbins, 
                                          self.sigma, self.power, 
                                          c1=self.fibmodel_slope, 
                                          c2=self.fibmodel_intercept,
                                          break_fix=self.fibmodel_breakpoint,
                                          fib_group=self.fiber_group, 
                                          col_group=self.col_group,
                                          cols=np.arange(0, self.D, 
                                                         self.fibmodel_step),
                                          kind=self.fibmodel_interpkind,
                                          do_fit=True)

            self.fill_in_dead_fibers(['core', 'fibmodel'])
            for fib, fiber in enumerate(self.dead_fibers):
                fiber.spectrum = np.zeros((self.D,))
            get_indices(self.image, self.dead_fibers, self.fsize)
                        
        else:
            self.load_fibmodel(path='calpath')
                
        if self.check_fibermodel:
            outfile = op.join(self.path,'fibmodel_%s.png' %self.basename)
            check_fiber_profile(self.image, self.good_fibers, outfile, self.fsize)
        
        if self.make_model_image:
            for fib,fiber in enumerate(self.fibers):
                fiber.fib_mult = np.ones((self.D,))
            self.fibmodel_image = get_model_image(self.image, self.fibers, 
                                            'fib_mult', debug=False)

    def fiberextract(self, cols=None):
        '''
        This function gets the spectrum for each fiber.  It checks 
        functional dependencies first: prepare_image(), get_trace(), and
        get_fibermodel(). 

        '''
        if not self.image_prepped:
            self.prepare_image()
        if not self.fibers:
            self.get_trace()
        if self.fibers[0].core is None:
            self.get_fibermodel()
                
        self.log.info('Calculating the spectrum for each fiber for %s' %self.basename)
        new_fast_norm(self.image, self.fibers, cols=cols, mask=self.mask)        
        
    
    def get_wavelength_solution(self):
        '''
        This function gets the wavelength solution for each fiber.  It checks 
        functional dependencies first: prepare_image(), get_trace(),
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
        if not self.image_prepped:
            self.prepare_image()
        if not self.fibers:
            self.get_trace()
        if self.fibers[0].core is None:
            self.get_fibermodel()
        if self.type == 'twi' or self.refit:
            if self.fibers[0].spectrum is None:
                self.fiberextract()
            if self.init_lims is None:
                self.log.warning("Please provide initial wavelength endpoint guess")
                return None
            self.log.info('Calculating the wavelength solution for %s' %self.basename)
            for k in xrange(2):
                content=False
                cnting = 0
                while content==False:
                    fiber = self.good_fibers[self.default_fib]
                    fw, fwp = calculate_wavelength_chi2(np.arange(self.D), 
                                                        solar_spec,
                                                        self.good_fibers,
                                                        self.default_fib,
                                                        self.fiber_group,
                                                      init_lims=self.init_lims, 
                                                        debug=False, 
                                                  interactive=self.interactive,
                                                        nbins=self.wave_nbins,
                                                        res=self.wave_res,
                                                        init_sol=self.init_sol)
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
                        self.log.warning("Ran through 9 loops but no solution found.")
                        return None
                fc = np.arange(len(self.good_fibers))
                if self.default_fib==0:
                    fibs1 = []
                else:
                    fibs1 = fc[(self.default_fib-1)::-1]
                if self.default_fib==(len(self.good_fibers)-1):
                    fibs2 = []
                else:
                    fibs2 = fc[(self.default_fib+1)::1]
                for fibn in fibs1:
                    fiber = self.good_fibers[fibn]
                    fw, fwp = calculate_wavelength_chi2(np.arange(self.D), 
                                                        solar_spec,
                                                        self.good_fibers,
                                                        fibn,
                                                        self.fiber_group,
                                                      init_lims=self.init_lims, 
                                                        debug=False, 
                                                       interactive=False,
                                     init_sol=self.fibers[fibn+1].wave_polyvals,
                                                        nbins=self.wave_nbins,
                                                        res=self.wave_res)
                    fiber.wavelength = fw*1.
                    fiber.wave_polyvals = fwp*1.
                for fibn in fibs2:
                    fiber = self.good_fibers[fibn]
                    fw, fwp = calculate_wavelength_chi2(np.arange(self.D), 
                                                        solar_spec,
                                                        self.good_fibers,
                                                        fibn,
                                                        self.fiber_group,
                                                      init_lims=self.init_lims, 
                                                        debug=False, 
                                                       interactive=False,
                                     init_sol=self.fibers[fibn-1].wave_polyvals,
                                                        nbins=self.wave_nbins,
                                                        res=self.wave_res)
                    fiber.wavelength = fw*1.
                    fiber.wave_polyvals = fwp*1.
                if k==0:
                    self.get_master_sky(norm=True)
                    solar_spec = np.zeros((len(self.masterwave),2))
                    solar_spec[:,0] = self.masterwave
                    solar_spec[:,1] = self.mastersmooth
            #for k in np.arange(2):        
            #    offset = get_wavelength_offsets(self.good_fibers, 25)
            #    for off, fiber in zip(offset, self.good_fibers):
            #        self.log.info('Fiber wavelength offset is: %0.2f' % off)
            #        fiber.wavelength = fiber.wavelength - off
            for k in np.arange(self.good_fibers[0].D):
                x = [fiber.trace[k] for fiber in self.good_fibers]
                y = [fiber.wavelength[k] for fiber in self.good_fibers]
                p0 = np.polyfit(x,y,3)
                z = np.polyval(p0,x)
                for i,fiber in enumerate(self.good_fibers):
                    fiber.wavelength[k] = z[i]
            self.fill_in_dead_fibers(['wavelength', 'wave_polyvals'])        
            
        else:
            self.load(path='calpath', spec_list=['wavelength'])       
            
        if self.check_wave:
            if self.fibers[0].spectrum is None:
                self.fiberextract()
            outfile = op.join(self.path,'wavesolution_%s.png' %self.basename)
            check_wavelength_fit(self.good_fibers, solar_spec, outfile)
                
                    
    def get_fiber_to_fiber(self):
        '''
        This function gets the fiber to fiber normalization for this amplifier. 
        It checks functional dependencies first: prepare_image(), get_trace(), 
        get_fibermodel(), fiberextract(), and get_wavelength_solution().
        '''
        if not self.image_prepped:
            self.prepare_image()
        if not self.fibers:
            self.get_trace()
        if self.fibers[0].core is None:
            self.get_fibermodel()
        if self.type == 'twi' or self.refit:
            if self.fibers[0].spectrum is None:
                self.fiberextract()
            if self.fibers[0].wavelength is None:
                self.get_wavelength_solution()
            self.log.info('Getting Fiber to Fiber for %s' %self.basename)
            
            masterwave = []
            masterspec = []
            k=0
            inds = []
            for fib, fiber in enumerate(self.good_fibers):
                masterwave.append(fiber.wavelength)
                masterspec.append(fiber.spectrum)
                inds.append(k)
                k += self.D
            masterwave = np.hstack(masterwave)
            masterspec = np.hstack(masterspec)
            B, c = bspline_x0(masterwave, nknots=self.nknots)

            smooth = []
            for fiber,ind in zip(self.good_fibers,inds):
                m = medfilt(fiber.spectrum,9)
                sel = np.where(is_outlier(m-fiber.spectrum)==0)[0]
                sol = np.linalg.lstsq(c[int(ind):int(ind+self.D),:][sel,:],
                                      fiber.spectrum[sel])[0]
                smooth.append(np.dot(c[int(ind):int(ind+self.D),:],sol))
 
            wv = np.arange(masterwave.min(),masterwave.max()+self.wavestepsize,
                           self.wavestepsize)
            asmooth = np.hstack(smooth)
            avgsmooth = biweight_bin(wv,masterwave,asmooth)
            for fiber, sm in zip(self.good_fibers, smooth):
                sm_tot = np.interp(fiber.wavelength,wv,avgsmooth)
                fiber.fiber_to_fiber = sm / sm_tot   
            self.get_bspline_sky()
            k=0
            #B, c = bspline_x0(masterwave, nknots=21)
            for fiber,ind in zip(self.good_fibers,inds):
                spec = np.interp(fiber.wavelength ,self.masterwave, 
                                 self.mastersky)
                m = medfilt(fiber.spectrum/spec,25)
                sel = np.where(is_outlier(m-fiber.spectrum/spec)==0)[0]
                sel = np.intersect1d(sel,np.arange(5,self.D-5))
                sel2 = np.setdiff1d(np.arange(self.D),sel)
                sol=np.linalg.lstsq(c[int(ind):int(ind+self.D),:][sel,:],
                                    (fiber.spectrum/spec)[sel])[0]
                fiber.fiber_to_fiber = np.dot(c[int(ind):int(ind+self.D),:], 
                                              sol)
#                plt.figure(figsize=(8,6))
#                plt.scatter(fiber.wavelength,fiber.spectrum/spec,alpha=0.2)
#                plt.scatter(fiber.wavelength[sel2],
#                            fiber.spectrum[sel2]/spec[sel2],marker='x',
#                            color='k')
#                plt.plot(fiber.wavelength,m,'r--')
#                plt.plot(fiber.wavelength,fiber.fiber_to_fiber,'k-')
#                plt.axis([self.init_lims[0],self.init_lims[1],0.6,1.4])
#                plt.savefig(op.join(self.path,'test_%s_%i.png'%(self.amp,k)))
#                plt.close()
                k+=1
            for fib,fiber in enumerate(self.dead_fibers):
                fiber.fiber_to_fiber = np.zeros(fiber.spectrum.shape)
                
            
        else:
            self.load(path='calpath', spec_list=['fiber_to_fiber'])
            if self.adjust_ftf:
                self.log.info('Adjusting Fiber to Fiber from file')
                try:
                    Fcor = np.loadtxt(op.join(self.virusconfig, 'FtFcor', 
                                         '%s%s.sdat' %(self.ifuslot,self.amp)))
                    for fiber,cor in zip(self.fibers,Fcor):
                        fiber.fiber_to_fiber = fiber.fiber_to_fiber + cor[1]
                        self.log.info('Adjusting fiber %i, by %0.3f' 
                                      %(fiber.fibnum,cor[1]))    
                except:
                    self.log.warning('No fiber to fiber file found.')
                
  
            
        
    def sky_subtraction(self):
        '''
        This function gets the master sky spectrum and 
        evaluates the sky_spectrum for each fiber. It then builds a sky image
        and a sky-subtracted image.  It checks functional dependencies first: 
        prepare_image(), get_trace(), get_fibermodel(), fiberextract(), 
        get_wavelength_solution(), and fiber_to_fiber().        
        '''
        if not self.image_prepped:
            self.prepare_image()
        if not self.fibers:
            self.get_trace()
        if self.fibers[0].core is None:
            self.get_fibermodel()
        if self.fibers[0].fiber_to_fiber is None:
            self.get_fiber_to_fiber()
        if self.fibers[0].spectrum is None:
            self.fiberextract()
        if self.fibers[0].wavelength is None:
            self.get_wavelength_solution()
        if self.skypath is not None:
            self.load(path='skypath', spec_list=['sky_spectrum'])       
            if self.fibers[0].sky_spectrum is None:
                self.log.warning("Loading sky spectrum from %s did not work." 
                      %self.skypath)
                self.log.info("Setting skypath to None for this amplifier.")
                self.skypath = None
            else: 
                for fiber in self.fibers:
                    fiber.sky_spectrum *= self.sky_scale
                self.build_mastersky()
        self.log.info('Subtracting sky for %s' %self.basename)
        if self.skypath is None:
            self.alt_sky_model()                
            for fib, fiber in enumerate(self.fibers):
                fiber.sky_spectrum = (fiber.fiber_to_fiber 
                                 * np.interp(fiber.wavelength, self.masterwave, 
                                             self.mastersky))
        for fib, fiber in enumerate(self.fibers):    
            fiber.sky_subtracted = fiber.spectrum - fiber.sky_spectrum 
            if hasattr(fiber, 'twi_spectrum'):
                
                fiber.corrected_spectrum = np.where(fiber.twi_spectrum>0.0, 
                                                fiber.spectrum 
                                                / biweight_filter(fiber.twi_spectrum,self.filt_size_final),
                                                0.0)
            
                fiber.corrected_sky_subtracted = np.where(fiber.twi_spectrum>0.0, 
                                                fiber.sky_subtracted 
                                                / biweight_filter(fiber.twi_spectrum,self.filt_size_final),
                                                0.0)  
            else:
                fiber.corrected_spectrum = np.where(fiber.fiber_to_fiber>0.0, 
                                                fiber.spectrum 
                                                / fiber.fiber_to_fiber,
                                                0.0)
            
                fiber.corrected_sky_subtracted = np.where(fiber.fiber_to_fiber>0.0, 
                                                fiber.sky_subtracted 
                                                / fiber.fiber_to_fiber,
                                                0.0) 
        if self.make_skyframe: 
            self.skyframe = get_model_image(self.image, self.fibers, 
                                            'sky_spectrum', debug=False)
            self.flat_image = get_model_image(self.image, self.fibers, 
                                            'fiber_to_fiber', debug=False)
            self.clean_image = self.image - self.skyframe
        if self.do_cont_sub:
            for fib, fiber in enumerate(self.fibers):
                fiber.continuum = biweight_filter(fiber.spectrum 
                                                  -fiber.sky_spectrum, 
                                                  self.cont_smooth, 
                                                  ignore_central=7)
            self.cont_frame = get_model_image(self.image, self.fibers, 
                                            'continuum', debug=False)
            self.continuum_sub = self.image - self.skyframe - self.cont_frame
        if self.make_residual:
            self.model = get_model_image(self.image, self.fibers, 'spectrum', 
                                         debug=False)
            self.residual = self.image - self.model
        #self.plot_fib_sky()
        
    def clean_cosmics(self):
        '''
        We use a direct copy of Malte Tewes and Pieter Van Dokkum's cosmics.py
        which is a python-interface of the LA cosmics algorithm to remove
        cosmics from sky-subtracted frames.
        
        '''
        if self.cosmic_iterations>0:
            self.log.info('Cleaning cosmics for %s' %self.basename)
            try:
                cc = cosmics.cosmicsimage(self.clean_image, gain=1.0, 
                                      readnoise=self.rdnoise, 
                                      sigclip=25.0, sigfrac=0.001, objlim=0.001,
                                      satlevel=-1.0)
            except:
                self.log.info('Running clean cosmics on regular image since '
                              'sky subtraction has not been performed.')
                cc = cosmics.cosmicsimage(self.image, gain=1.0, 
                                      readnoise=self.rdnoise, 
                                      sigclip=25.0, sigfrac=0.001, objlim=0.001,
                                      satlevel=-1.0)
            cc.run(maxiter=self.cosmic_iterations)
            c = np.where(cc.mask == True)
            if self.mask is None:
                self.mask = np.zeros(self.image.shape)
            for x, y in zip(c[0], c[1]):
                self.mask[x][y] = -1.0
            self.error[c[0],c[1]] = -1.0
             
             
    def get_master_sky(self, sky=False, norm=False):
        '''
        This builds a master sky spectrum from the spectra of all fibers
        using a biweight average of a wavelength ordered master array.
        
        '''
        
        masterwave = []
        if sky:
            masterspec = []
        if norm:
            mastersmooth=[]
        for fib, fiber in enumerate(self.good_fibers):
            if norm:
                y = biweight_filter(fiber.spectrum, self.filt_size_ind)
                mastersmooth.append(y/fiber.spectrum)
                masterwave.append(fiber.wavelength)
            if sky:
                flag = np.ones(fiber.spectrum.shape,dtype=bool)
                flag[fiber.fiber_to_fiber<1e-8] = False
                flag[fiber.spectrum == 0.0] = False
                masterspec.append(fiber.spectrum[flag]
                                  /fiber.fiber_to_fiber[flag])
                masterwave.append(fiber.wavelength[flag])

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
                                                
    def get_significance_map(self, attr='clean_image'):
        '''
        Rudimentary significance map algorithm
        '''
        self.log.info('Calculating significance map for %s' %self.basename)
        self.sig, self.sigwave = calculate_significance(self.fibers, 
                                                        getattr(self,attr), 
                                                        self.error)
    def make_error_analysis(self, attr='clean_image'):
        wave_step = 4.
        wavelength = np.array([fiber.wavelength for fiber in self.fibers])
        ftf = np.array([fiber.fiber_to_fiber for fiber in self.fibers])
        trace = np.array([np.round(fiber.trace) for fiber in self.fibers], dtype=int)
                
        
        # Calculate Noise
        wvmin = wavelength.min()
        wvmax = wavelength.max()
        wave_array = np.arange(wvmin,wvmax+2*wave_step, wave_step)
        wave_array_fit = wave_array[:-1]+wave_step/2.
        xind=[]
        yind=[]
        lx = 0
        ly = len(wave_array)-1
        for i in np.arange(ly):
            x,y = np.where((wavelength>=wave_array[i]) 
                            * (wavelength<wave_array[i+1]))
            lx = np.max([lx,len(x)])
            xind.append(x)
            yind.append(y)
        A = np.ones((ly,lx))*-999.
        B = np.ones((ly,lx))*-999.
        for i in np.arange(ly):
            A[i,:len(xind[i])] = np.where(ftf[xind[i],yind[i]]>1e-8, 
                                          getattr(self,attr)[trace[xind[i],
                                                                 yind[i]],
                                                           yind[i]]
                                          / ftf[xind[i],yind[i]],
                                          0.0)
            B[i,:len(xind[i])] = self.error[trace[xind[i],yind[i]],yind[i]]
        C = np.ma.array(A, mask=(A==-999.).astype(np.int))    
        D = np.ma.array(B, mask=(B==-999.).astype(np.int))    
        E = biweight_midvariance(C, axis=(1,))
        F = biweight_location(D, axis=(1,))
        self.error_analysis = np.zeros((3,len(wave_array_fit)))
        self.error_analysis[0,:] = wave_array_fit
        self.error_analysis[1,:] = F 
        self.error_analysis[2,:] = E 
        
    def get_binned_sky(self):
        masterwave,masterspec = [],[]
        for fiber in self.good_fibers:
            masterwave.append(fiber.wavelength)
            masterspec.append(fiber.spectrum/fiber.fiber_to_fiber)
        masterwave = np.hstack(masterwave)
        masterspec = np.hstack(masterspec)

        self.masterwave = np.arange(masterwave.min(),masterwave.max()
                                               +self.wavestepsize,
                               self.wavestepsize)
                               
        self.mastersky = biweight_bin(self.masterwave,masterwave,
                                      masterspec)                                                         

    def get_bspline_sky(self):
        self.log.info('Calculating master sky using bspline approximation')
        masterwave,masterspec,masterfib = [],[],[]
        for i,fiber in enumerate(self.good_fibers):
            masterwave.append(fiber.wavelength)
            masterspec.append(fiber.spectrum/fiber.fiber_to_fiber)
            masterfib.append(i*np.ones(fiber.wavelength.shape, dtype=int))
        masterwave = np.hstack(masterwave)
        masterspec = np.hstack(masterspec)
        masterfib = np.hstack(masterfib)
        ind = np.argsort(masterwave)
        masterwave = masterwave[ind]
        masterspec = masterspec[ind]
        masterfib = masterfib[ind]
        self.masterwave = masterwave
        bins = np.arange(masterwave.min(), 
                         masterwave.max()+self.bspline_binsize, 
                         self.bspline_binsize)
        self.mastersky = np.zeros(masterwave.shape)
        for i in np.arange(len(bins)-1):
            sel = np.where((masterwave>=bins[i])*(masterwave<bins[i+1]))[0]
            B,c = bspline_x0(masterwave[sel],
                             nknots=int(self.bspline_binsize
                                      /self.bspline_waveres))
            m = percentile_filter(masterspec[sel], 50, 15)
            d = masterspec[sel]-m
            bv = biweight_midvariance(d)
            nothigh = np.abs(d) < (self.sky_sigthresh*bv)
            sel1 = np.where(nothigh)[0]

            sol = np.linalg.lstsq(c[sel1,:],masterspec[sel][sel1])[0]
            for j in np.arange(self.sky_iterations):
                v = np.dot(c,sol)
                d = masterspec[sel] - v
                bv = biweight_midvariance(d)
                sel1 = np.abs(d) < (self.sky_sigthresh * bv)
                sol = np.linalg.lstsq(c[sel1,:],masterspec[sel][sel1])[0]
            self.mastersky[sel] = np.dot(c,sol) 
        
     
    def alt_sky_model(self):
        fac = 10
        mn = 1e9
        mx = 0.
        for fiber in self.good_fibers:
            mn = np.min([mn, fiber.wavelength.min()])
            mx = np.max([mx, fiber.wavelength.max()])
        xs = np.linspace(0, 1, self.D*fac)
        A = np.zeros((len(xs), len(self.good_fibers)))
        for i, fiber in enumerate(self.good_fibers):
            y = fiber.spectrum / fiber.fiber_to_fiber
            xp = np.interp(fiber.wavelength, np.linspace(mn, mx, self.D*fac),
                           xs, left=0.0, right=0.0)
            tck = splrep(xp, y)
            A[:, i] = splev(xs, tck)
        ys = biweight_location(A, axis=(1,))
        self.masterwave = np.linspace(mn, mx, self.D*fac)
        B, c = bspline_x0(self.masterwave, nknots=self.D)
        sol = np.linalg.lstsq(c, ys)[0]
        self.mastersky = np.dot(c, sol)
        
    def build_mastersky(self):
        masterwave = []
        mastersky = []
        for fiber in self.good_fibers:
            masterwave.append(fiber.wavelength)
            mastersky.append(fiber.sky_spectrum/fiber.fiber_to_fiber)
        masterwave = np.hstack(masterwave)
        mastersky = np.hstack(mastersky)
        ind = np.argsort(masterwave)
        self.masterwave = masterwave[ind]
        self.mastersky = mastersky[ind]
        
    def plot_fib_sky(self):
        s = np.arange(self.masterwave.min(),
                      self.masterwave.max()+self.bspline_binsize, 
                      self.bspline_binsize)
        nfibs = len(self.fibers)
        colors = plt.get_cmap('RdBu')(np.arange(nfibs)/(nfibs-1.))

        plt.figure(figsize=(8,6))
        ax1 = plt.subplot(2,1,1)
        ax2 = plt.subplot(2,1,2)
        ax1.set_position([0.15,0.35,0.7,0.5])
        ax2.set_position([0.15,0.15,0.7,0.2])
        
        for fib,fiber in enumerate(self.fibers):
            v = np.interp(fiber.wavelength,self.masterwave,self.mastersky)
            ax1.scatter(fiber.wavelength, 
                        fiber.spectrum/fiber.fiber_to_fiber,
                        s=5, color=colors[fib])
            ax2.scatter(fiber.wavelength, 
                        (fiber.spectrum/fiber.fiber_to_fiber - v) / v,
                        s=10, color=colors[fib],alpha=0.3)
        ax1.plot(self.masterwave,self.mastersky,'k-')
        ax2.plot(self.init_lims[0],self.init_lims[1],[0,0],'k--')
        for v in s[:-1]:                
            ax1.set_xlim([v,v+self.bspline_binsize])
            ax2.set_xlim([v,v+self.bspline_binsize])
            xl = np.searchsorted(self.masterwave,v)
            xh = np.searchsorted(self.masterwave,v+270)
            mn = self.mastersky[xl:xh].min()
            mx = self.mastersky[xl:xh].max()
            ax1.set_ylim([mn-(mx-mn)*0.1, mn+(mx-mn)*1.1])
            ax2.set_ylim([-0.05,0.05])
            plt.savefig(op.join(self.path,'testfib_%s_%0.1f.png'%(self.amp,v)))
        plt.close()