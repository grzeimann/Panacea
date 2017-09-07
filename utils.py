# -*- coding: utf-8 -*-
"""
Created on Wed Oct  5 18:12:13 2016

@author: gregz

This code relies on original software from:
Copyright (c) 2011-2016, Astropy Developers    
Copyright (c) 2012, Free Software Foundation    

"""

import numpy as np
import sys

def median_absolute_deviation(a, axis=None):
    """
    Copyright (c) 2011-2016, Astropy Developers    
    
    Compute the median absolute deviation.

    Returns the median absolute deviation (MAD) of the array elements.
    The MAD is defined as ``median(abs(a - median(a)))``.

    Parameters
    ----------
    a : array-like
        Input array or object that can be converted to an array.
    axis : tuple, optional
        Axis along which the medians are computed. The default (axis=None)
        is to compute the median along a flattened version of the array.

    Returns
    -------
    median_absolute_deviation : ndarray
        A new array holding the result. If the input contains
        integers, or floats of smaller precision than 64, then the output
        data-type is float64.  Otherwise, the output data-type is the same
        as that of the input.

    Examples
    --------

    This will generate random variates from a Gaussian distribution and return
    the median absolute deviation for that distribution::

        >>> from utils import median_absolute_deviation
        >>> from numpy.random import randn
        >>> randvar = randn(10000)
        >>> mad = median_absolute_deviation(randvar)

    See Also
    --------
    numpy.median
    
    Note
    --------
    Copy of the astropy function with the "axis" argument added appropriately.

    """

    # Check if the array has a mask and if so use np.ma.median
    # See https://github.com/numpy/numpy/issues/7330 why using np.ma.median
    # for normal arrays should not be done (summary: np.ma.median always
    # returns an masked array even if the result should be scalar). (#4658)
    if isinstance(a, np.ma.MaskedArray):
        func = np.ma.median
    else:
        func = np.nanmedian

    a = np.asanyarray(a)

    a_median = func(a, axis=axis)

    # re-broadcast the output median array to subtract it
    if axis is not None:
        for i in axis:
            a_median = np.expand_dims(a_median, axis=i)

    # calculated the median average deviation
    return func(np.abs(a - a_median), axis=axis)
    
def biweight_bin(xv, x, y):
    '''
    Compute the biweight location with a moving window of size "order"

    '''
    diff_array = np.hstack([np.diff(xv),np.diff(xv)[-1]]) / 2.
    mxcol = 0
    sel_list = []
    for i,v in enumerate(xv):
        sel_list.append(np.where( (x>(xv[i]-diff_array[i])) 
                                 *(x<(xv[i]+diff_array[i])))[0]) 
        mxcol = np.max([mxcol, len(sel_list[-1])])
    
    A = np.ones((len(xv),mxcol))*-999.
    for i,v in enumerate(xv):
        A[i,:len(sel_list[i])] = y[sel_list[i]]
    C = np.ma.array(A, mask=(A == -999.).astype(np.int))
    return biweight_location(C, axis=(1,))


def biweight_filter2d(a, Order, Ignore_central=(3,3), c=6.0, M=None, func=None):
    '''
    Compute the biweight location with a moving window of size "order"

    '''
    if not isinstance(Order, tuple):
        print("The order should be an tuple")
        sys.exit(1)
    if not isinstance(Ignore_central, tuple):
        print("The ignore_central value should be an tuple")
        sys.exit(1)   
    for order in Order:
        if order%2==0:
            order+=1
    for ignore_central in Ignore_central:
        if ignore_central%2==0:
            ignore_central+=1
    if ((Order[0]-3 <= Ignore_central[0] and Order[0]>1) 
        or (Order[1]-3 <= Ignore_central[1] and Order[1]>1)):
        print("The max order-3 should be larger than max ignore_central.")
        sys.exit(1)
    if func is None:
        func = biweight_location
        
    a = np.array(a, copy=False)
    if a.ndim != 2:
        print("Input array/list should be 2-dimensional")
        sys.exit()

    yc = np.arange(Ignore_central[0]/2+1,Order[0]/2+1)
    xc = np.arange(Ignore_central[1]/2+1,Order[1]/2+1)
    ly = np.max([len(yc)*2,1])
    lx = np.max([len(xc)*2,1])
    size = lx * ly
    A = np.ones(a.shape + (size,))*-999.
    k=0
    if Order[0]>1:
        if Order[1]>1:
            for i in yc:
                for j in xc:
                    A[:-i,:-j,k] = a[i:,j:]
                    k+=1
            for i in yc:
                for j in xc:
                    A[i:,j:,k] = a[:-i,:-j]
                    k+=1
            for i in yc:
                for j in xc:
                    A[:-i,j:,k] = a[i:,:-j]
                    k+=1
            for i in yc:
                for j in xc:
                    A[i:,:-j,k] = a[:-i,j:]
                    k+=1
        else:
            for i in yc:
                A[:-i,:,k] = a[i:,:]
                k+=1
            for i in yc:
                A[i:,:,k] = a[:-i,:]
                k+=1
    else:
        for j in xc:
            A[:,:-j,k] = a[:,j:]
            k+=1
        for j in xc:
            A[:,j:,k] = a[:,:-j]
            k+=1    
    C = np.ma.array(A, mask=(A == -999.).astype(np.int))
    return func(C, axis=(2,))

def biweight_filter(a, order, ignore_central=3, c=6.0, M=None, func=None):
  
    if order%2==0:
        order+=1
    if ignore_central%2==0:
        ignore_central+=1
    if func is None:
        func = biweight_location
    if ignore_central+3 >= order:
        print('The order, %i, should be +3 higher than ignore_central, %i.' 
              %(order, ignore_central))
        sys.exit(1)
    a = np.array(a, copy=False)
    if ignore_central > 0:
        yc = np.arange(ignore_central/2+1,order/2+1)
    else:
        yc = np.arange(0,order/2+1)
    ly = np.max([len(yc)*2,1])
    A = np.ones(a.shape + (ly,))*-999.
    k=0
    if order>1:
        for i in yc:
            A[:-i,k] = a[i:]
            k+=1
        for i in yc:
            A[i:,k] = a[:-i]
            k+=1  
    C = np.ma.array(A, mask=(A == -999.).astype(np.int))
    return func(C, axis=(1,))


#def biweight_filter(a, order, ignore_central=3, c=6.0, M=None, func=None):
#    '''
#    Compute the biweight location with a moving window of size "order"
#
#    '''
#    if not isinstance(order, int):
#        print("The order should be an integer")
#        sys.exit(1)
#    if not isinstance(ignore_central, int):
#        print("The ignore_central value should be an integer")
#        sys.exit(1)        
#    if order%2==0:
#        order+=1
#    if ignore_central%2==0:
#        ignore_central+=1
#    if order-3 <= ignore_central:
#        print("The order-3 should be larger than ignore_central.")
#        sys.exit(1)
#    if func is None:
#        func = biweight_location
#    a = np.array(a, copy=False)
#    if a.ndim != 1:
#        print("Input array/list should be 1-dimensional")
#        sys.exit()
#    ignore = [order/2]
#    for i in xrange(ignore_central/2):
#        ignore.append(order/2 - i - 1)
#        ignore.append(order/2 + i + 1)
#    half_order = order / 2    
#    order_array = np.delete(np.arange(order), ignore)
#    A = np.zeros((len(a)-order+1, len(order_array)))
#    k=0
#    for i in order_array:
#        if (order-i-1) == 0:
#            A[:,k] = a[i:]
#        else:
#            A[:,k] = a[i:-(order-i-1)]
#        k+=1
#    
#    Ab = func(A, axis=(1,))
#    A1 = np.zeros((half_order,))
#    A2 = np.zeros((half_order,))
#    for i in xrange(half_order):
#        ignore_l = [i]
#        ignore_h = [half_order]
#        for j in xrange(ignore_central/2):
#            if (i - j - 1) >=0:
#                ignore_l.append(i - j - 1)
#            ignore_l.append(i + j + 1)
#            if i > j:
#                ignore_h.append(half_order+j+1)
#            ignore_h.append(half_order-j-1)
#        A1[i] = func(np.delete(a[:(half_order+i+1)], ignore_l))
#        A2[-(i+1)] = func(np.delete(a[-(half_order+i+1):], ignore_h))
#    return np.hstack([A1,Ab,A2])
    
    

def biweight_location(a, c=6.0, M=None, axis=None, eps=1e-8):
    """
    Copyright (c) 2011-2016, Astropy Developers        
    
    Compute the biweight location for an array.

    Returns the biweight location for the array elements.
    The biweight is a robust statistic for determining the central
    location of a distribution.

    The biweight location is given by the following equation

    .. math::

        C_{bl}= M+\\frac{\Sigma_{\|u_i\|<1} (x_i-M)(1-u_i^2)^2}
        {\Sigma_{\|u_i\|<1} (1-u_i^2)^2}

    where M is the sample mean or if run iterative the initial guess,
    and u_i is given by

    .. math::

      u_{i} = \\frac{(x_i-M)}{cMAD}

    where MAD is the median absolute deviation.

    For more details, see Beers, Flynn, and Gebhardt, 1990, AJ, 100, 32B

    Parameters
    ----------
    a : array-like
        Input array or object that can be converted to an array.
    c : float, optional
        Tuning constant for the biweight estimator.  Default value is 6.0.
    M : float, optional
        Initial guess for the biweight location.
    axis : tuple, optional
        tuple of the integer axis values ot calculate over.  Should be sorted.
    Returns
    -------
    biweight_location : float
        Returns the biweight location for the array elements.

    Examples
    --------

    This will generate random variates from a Gaussian distribution and return
    the biweight location of the distribution::

    >>> from utils import biweight_location
    >>> from numpy.random import randn
    >>> randvar = randn(10000)
    >>> cbl = biweight_location(randvar)

    See Also
    --------
    median_absolute_deviation, biweight_midvariance
    
    Note
    --------
    Copy of the astropy function with the "axis" argument added appropriately.
    """


    if M is None:
        if isinstance(a, np.ma.MaskedArray):
            func = np.ma.median
        else:
            a = np.array(a, copy=False)
            func = np.median
        M = func(a, axis=axis)
    else:
        a = np.array(a, copy=False)
                
    N = M*1.      
    # set up the difference
    if axis is not None:
        for i in axis:
            N = np.expand_dims(N, axis=i)
            
    d = a - N
    
    # set up the weighting
    if axis is not None:
        MAD = median_absolute_deviation(a, axis=axis)
        for i in axis:
            MAD = np.expand_dims(MAD, axis=i)
    else:
        MAD = median_absolute_deviation(a)

    u = np.where(MAD < eps, 0., d / c / MAD)
    
    # now remove the outlier points
    if isinstance(a, np.ma.MaskedArray):
        mask = (np.abs(u) < 1).astype(np.int) * (1-a.mask.astype(np.int))
    else:
        mask = (np.abs(u) < 1).astype(np.int)
    u = (1 - u ** 2) ** 2
    return M + (d * u * mask).sum(axis=axis) / (u * mask).sum(axis=axis)
    
    
def biweight_midvariance(a, c=15.0, M=None, axis=None, eps=1e-8, niter=1):
    """
    Copyright (c) 2011-2016, Astropy Developers    
    
    Compute the biweight midvariance for an array.

    Returns the biweight midvariance for the array elements.
    The biweight midvariance is a robust statistic for determining
    the midvariance (i.e. the standard deviation) of a distribution.

    The biweight location is given by the following equation

    .. math::

      C_{bl}= (n')^{1/2} \\frac{[\Sigma_{|u_i|<1} (x_i-M)^2(1-u_i^2)^4]^{0.5}}
      {|\Sigma_{|u_i|<1} (1-u_i^2)(1-5u_i^2)|}

    where :math:`u_i` is given by

    .. math::

      u_{i} = \\frac{(x_i-M)}{cMAD}

    where MAD is the median absolute deviation.

    :math:`n'` is the number of data for which :math:`|u_i| < 1` holds, while the
    summations are over all i up to n:

    .. math::

        n' = \Sigma_{|u_i|<1}^n 1

    This is slightly different than given in the reference below, but
    results in a value closer to the true midvariance.

    The midvariance parameter c is typically 9.0.

    For more details, see Beers, Flynn, and Gebhardt, 1990, AJ, 100, 32B

    Parameters
    ----------
    a : array-like
        Input array or object that can be converted to an array.
    c : float
        Tuning constant for the biweight estimator.  Default value is 9.0.
    M : float, optional
        Initial guess for the biweight location.
    axis : tuple, optional
        tuple of the integer axis values ot calculate over.  Should be sorted.

    Returns
    -------
    biweight_midvariance : float
        Returns the biweight midvariance for the array elements.

    Examples
    --------

    This will generate random variates from a Gaussian distribution and return
    the biweight midvariance of the distribution::

    >>> from utils import biweight_midvariance
    >>> from numpy.random import randn
    >>> randvar = randn(10000)
    >>> scl = biweight_midvariance(randvar)

    See Also
    --------
    median_absolute_deviation, biweight_location
    
    Note
    --------
    Copy of the astropy function with the "axis" argument added appropriately.
    """
    for k in np.arange(niter):
        if isinstance(a, np.ma.MaskedArray):
            func = np.ma.median
        else:
            a = np.array(a, copy=False)
            func = np.median
        if M is None or k>0:
            M = func(a, axis=axis)
        N = M*1.      
        # set up the difference
        if axis is not None:
            for i in axis:
                N = np.expand_dims(N, axis=i)
    
        # set up the difference
        d = np.asarray(a - N)
    
        # set up the weighting
        if axis is not None:
            MAD = median_absolute_deviation(a, axis=axis)
            for i in axis:
                MAD = np.expand_dims(MAD, axis=i)
        else:
            MAD = median_absolute_deviation(a)
        # set up the weighting
        u = np.where(MAD < eps, 0., d / c / MAD)
    
        # now remove the outlier points
        if isinstance(a, np.ma.MaskedArray):
            mask = (np.abs(u) < 1).astype(np.int) * (1-a.mask.astype(np.int))
            a.mask = 1 - mask
        else:
            mask = (np.abs(u) < 1).astype(np.int)
    u = u ** 2
    n = mask.sum(axis=axis)
    return n ** 0.5 * (mask * d**2 * (1 - u) ** 4).sum(axis=axis) ** 0.5\
        / np.abs((mask * (1 - u) * (1 - 5 * u)).sum(axis=axis))


def is_outlier(points, thresh=3.5):
    """
    Copyright (c) 2012, Free Software Foundation   
    
    
    Returns a boolean array with True if points are outliers and False 
    otherwise.

    Parameters:
    -----------
        points : An numobservations by numdimensions array of observations
        thresh : The modified z-score to use as a threshold. Observations with
            a modified z-score (based on the median absolute deviation) greater
            than this value will be classified as outliers.

    Returns:
    --------
        mask : A numobservations-length boolean array.

    References:
    ----------
        Boris Iglewicz and David Hoaglin (1993), "Volume 16: How to Detect and
        Handle Outliers", The ASQC Basic References in Quality Control:
        Statistical Techniques, Edward F. Mykytka, Ph.D., Editor. 
    """
    if len(points.shape) == 1:
        points = points[:,None]
    median = np.median(points, axis=0)
    diff = np.sum((points - median)**2, axis=-1)
    diff = np.sqrt(diff)
    med_abs_deviation = np.median(diff)

    modified_z_score = 0.6745 * diff / med_abs_deviation

    return modified_z_score > thresh



   
def matrixCheby2D_7(x, y):
    if isinstance(x, (tuple, list)):
        x = np.asarray(x)
    if isinstance(y, (tuple, list)):
        y = np.asarray(y)

    T2x = 2. * x**2 - 1.
    T3x = 4. * x**3 - 3. * x
    T4x = 8. * x**4 - 8. * x**2 + 1.
    T5x = 16. * x**5 - 20. * x**3 + 5. * x
    T6x = 32. * x**6 - 48. * x**4 + 18. * x**2 - 1.
    T7x = 64. * x**7 - 112. * x**5 + 56. * x**3 - 7. * x
    T2y = 2. * y**2 - 1.
    T3y = 4. * y**3 - 3. * y
    T4y = 8. * y**4 - 8. * y**2 + 1.
    T5y = 16. * y**5 - 20. * y**3 + 5. * y
    T6y = 32. * y**6 - 48. * y**4 + 18. * y**2 - 1
    T7y = 64. * y**7 - 112. * y**5 + 56. * y**3 - 7 * y
    
    return np.vstack((T7x, T6x, T5x, T4x, T3x, T2x, x, T7y, T6y, T5y, 
                      T4y, T3y, T2y, y, T6x*y, x*T6y, T5x*T2y, T2x*T5y,
                      T4x*T3y, T3x*T4y, T5x*y, x*T5y, T4x*T2y, T2x*T4y, 
                      T3x*T3y, T4x*y, x*T4y, T3x*T2y, T2x*T3y, T3x*y, 
                      x*T3y, T2x*T2y, T2x*y, x*T2y, x*y, np.ones(x.shape))).swapaxes(0,1)