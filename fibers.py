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

__all__ = ["Fibers"]

class Fibers:
    def __init__(self, N, D, trace_poly_order=3):
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
        :init weight:
            Fiber Model Weight over the CCD image.
        '''
        self.N = N
        self.D = D
        self.trace_poly_order = trace_poly_order
        self.flag = -99999
        self.trace_x = self.flag * np.ones((D,),dtype = np.int)
        self.trace_y = np.zeros((D,))
        self.trace = np.zeros((D,))
        self.polyvals = np.zeros((trace_poly_order,))
        self.wave = np.zeros((D,))
        self.throughput = np.zeros((D,))
        self.RA = None
        self.Dec = None
        self.rot = None
        self.theta = None
        self.gain = None
        self.dark_mult = None
        self.bias_mult = None
        self.file = None
        
    def fitpoly(self):
        sel = self.x != self.flag
        self.polyvals = np.polyfit(self.x[sel] / self.D, self.y[sel], 
                                   self.trace_poly_order)
        
    def evalpoly(self):
        self.trace = np.polyval(self.polyvals, np.arange(len(self.x)) / self.D)