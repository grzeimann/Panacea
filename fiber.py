#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Fiber Class
-----------
To be used in conjuction with IFU reduction code, Panacea


"""

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

import numpy as np
import cPickle as pickle
import os.path as op
import os
from utils import biweight_filter

__all__ = ["Fiber"]

class Fiber:
    def __init__(self, D, fibnum, path, filename, trace_poly_order=3, 
                 fibmodel_poly_order=3,wave_poly_order=3):
        ''' 
        Initialize class
        ----------------
        :param N:
            The number of rows in the spectrograph detector image.
        :param D:
            The number of columns in the spectrograph detector  image.
        :param trace_poly_order:
            The order of the polynomial used to fit the trace.   
        :init N:
            Same as described above.             
        :init D:
            Same as described above.
        :init trace_poly_order:
            Same as described above.
        :init flag:
            Flag value for undefined values.
        :init trace_x:
            Columns values for trace.
        :init trace_y:
            Row values for the trace.
        :init trace: 
            Polynomial fitted values for the trace.
        :init polyvals:
            Polynomial coefficients for trace fit.
        :init wavelength:
            Wavelength for each column along the trace.
        :init throughput:
            Fiber throughput.
        '''
        self.D = D
        self.fibnum = fibnum
        self.filename = filename
        self.basename = op.basename(filename)[:-5]
        self.path = path
        self.trace_poly_order = trace_poly_order
        self.fibmodel_poly_order = fibmodel_poly_order
        self.wave_poly_order = wave_poly_order
        self.flag = -999
        self.init_trace_info()
        self.wave = None
        self.throughput = None
        self.fiber_to_fiber = None
        self.amp_to_amp = None
        self.RA = None
        self.Dec = None
        self.rot = None
        self.theta = None
        self.gain = None
        self.dark_mult = None
        self.bias_mult = None
        self.file = None
        self.ifuslot = None
        self.specid = None
        self.ifuid = None
        self.object = None
        self.datetime = None
        self.fibmodel = None
        self.fibmodel_x = None
        self.fibmodel_y = None
        self.binx = None
        self.fibmodel_polyvals = None
        self.trace_polyvals = None
        self.spectrum = None
        self.wave_polyvals = None
        
    def init_trace_info(self):
        self.trace_x = self.flag * np.ones((self.D,),dtype = np.int)
        self.trace_y = np.zeros((self.D,))
        self.trace = np.zeros((self.D,))
        self.trace_polyvals = np.zeros((self.trace_poly_order+1,))
        
    def fit_trace_poly(self):
        sel = self.trace_x != self.flag
        self.trace_polyvals = np.polyfit(self.trace_x[sel] / self.D, 
                                         self.trace_y[sel], 
                                         self.trace_poly_order)
                                   
    def fit_fibmodel_poly(self):
        self.fibmodel_polyvals = np.polyfit(self.fibmodel_x / self.D, 
                                            self.fibmodel_y, 
                                            self.fibmodel_poly_order)
                                   
    def eval_trace_poly(self, use_poly=False, smoothing_length=15):
        sel = self.trace_x != self.flag
        if use_poly:
            self.trace = np.polyval(self.trace_polyvals, 
                                    1. * np.arange(self.D) / self.D)
        else:
            self.trace = np.zeros((self.D,))
            init_x = np.where(sel)[0][0]
            fin_x = np.where(sel)[0][-1]
            self.trace[init_x:fin_x] = np.interp(np.arange(init_x,fin_x), 
                                                self.trace_x[sel],
                                                self.trace_y[sel])
            self.trace[sel] = biweight_filter(self.trace_y[sel], 
                                              smoothing_length)
            ix = int(init_x+smoothing_length/2+1)
            fx = int(init_x+smoothing_length/2+1 + smoothing_length*2)
            p1 = np.polyfit(np.arange(ix,fx), self.trace[ix:fx], 1)
            self.trace[:ix] = np.polyval(p1, np.arange(ix))
            ix = int(fin_x-smoothing_length/2-1 - smoothing_length*2)
            fx = int(fin_x-smoothing_length/2) 
            pf = np.polyfit(np.arange(ix,fx), self.trace[ix:fx], 1)
            self.trace[fx:self.D] = np.polyval(pf, np.arange(fx,self.D))
            
    
    def eval_fibmodel_poly(self, use_poly=False):
        porder, nbins = self.fibmodel_polyvals.shape
        self.fibmodel = np.zeros((self.D, nbins))
        for i in xrange(nbins):
            if use_poly:
                self.fibmodel[:,i] = np.polyval(self.fibmodel_polyvals[:,i],
                                                1.* np.arange(self.D) / self.D)
            else:
                self.fibmodel[:,i] = np.interp(np.arange(self.D), 
                                               self.fibmodel_x, 
                                               self.fibmodel_y[:,i])
    def eval_wave_poly(self):
        self.wavelength = np.polyval(self.wave_polyvals, 
                                     1.* np.arange(self.D) / self.D)        
        
    def save(self, specid, ifuslot, ifuid, amp):
        self.fibmodel = None
        self.trace_x = None
        self.trace_y = None
        #self.spectrum = None
        self.wavelength = None
        self.fn = op.join(self.path, 'fiber_%03d_%s_%s_%s_%s.pkl' 
                                % (self.fibnum, specid, ifuslot, ifuid, amp))
        if not op.exists(self.path):
            os.mkdir(self.path)
        with open(self.fn, 'wb') as f:
           pickle.dump(self, f)