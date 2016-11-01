# -*- coding: utf-8 -*-
"""
Fiber Utilities
-----------
To be used in conjuction with IFU reduction code, Panacea


"""

from utils import biweight_location
import numpy as np
import time

def get_trace_from_image(image, y_window=3, x_window=5, repeat_length=2,
                         order=3, mx_cut=0.1, debug=False):
    '''
    :param image:
        Image [FITS] from which to fit the trace.
    :param y_window:
        Within [-y_window, +y_window], we evaluate the index of the maximum
        value.  If that index is repeated :repeat_length: times then it is 
        logged as a potential maximum, unless it is already logged in which
        case it is skipped as to not repeat.
    :param x_window:
        We average the columns over [-x_window,+x_window] to avoid pixel 
        defects and other related issues.
    :param repeat_length:
        As discussed in :param y_window:, repeat_length is used to identify
        maxima via being the maximum point within [-y_window,+y_window] for a
        repeated length of :repeat_length:.
    :param order:
        The order of the polynomial used to smooth the trace via fitting.
    :param mx_cut:
        The cut used to define maxima from noise.  Must be greater than 
        mx_cut * <avg height of all peaks>.
    '''
    allfibers=[]
    a, b = image.shape
    x = np.arange(a)
    if debug:
        t1 = time.time()
    for i in xrange(x_window, b-1-x_window):
        y = biweight_location(image[:,(i-x_window):(i+x_window+1)], axis=(1,))
        mxval = np.zeros((a-2-y_window*2+1,))
        k=0
        peaks = []
        peaks_refined = []
        peaks_height = []
        for j in xrange(y_window,a - 1 - y_window):
            mxval[k] = np.argmax(y[(j-y_window):(j+y_window+1)])+j-y_window
            if k>repeat_length:
                if mxval[k] in mxval[(k-repeat_length):j]:
                    if mxval[k] not in peaks:
                        peaks.append(int(mxval[k]))
                        lwa = ((y[(i-y_window):(i+y_window+1)]
                                *x[(i-y_window):(i+y_window+1)]).sum() 
                               /(y[(i-y_window):(i+y_window+1)]).sum())
                        peaks_refined.append(lwa)
                        peaks_height.append(mxval[k])
            k+=1
        mh = max(peaks_height)
        peaks_final = [peaks_refined[j] for j,h in enumerate(peaks_height) 
                                        if h > (mx_cut * mh)]
    allfibers.append(peaks_final) 
    if debug:
        t2 = time.time()
        print("Time Taken for Trace: %0.2f" %(t2-t1))
    