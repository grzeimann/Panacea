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
        func = np.median

    a = np.asanyarray(a)

    a_median = func(a, axis=axis)

    # re-broadcast the output median array to subtract it
    if axis is not None:
        for i in axis:
            a_median = np.expand_dims(a_median, axis=i)

    # calculated the median average deviation
    return func(np.abs(a - a_median), axis=axis)
    

def biweight_filter(a, order, c=6.0, M=None):
    '''
    Compute the biweight location with a moving window of size "order"

    '''
    if not isinstance(order, int):
        print("The order should be an integer")
        sys.exit(1)
    if order%2==0:
        order+=1
    a = np.array(a, copy=False)
    if a.ndim != 1:
        print("Input array/list should be 1-dimensional")
        sys.exit()
    
    half_order = order / 2    
    A = np.zeros((len(a)-order+1,order))
    for i in xrange(order):
        if (order-i-1) == 0:
            A[:,i] = a[i:]
        else:
            A[:,i] = a[i:-(order-i-1)]
    
    Ab = biweight_location(A, axis=(1,))
    A1 = np.zeros((half_order,))
    A2 = np.zeros((half_order,))
    for i in xrange(half_order):
        A1[i] = biweight_location(a[:(half_order+i+1)])
        A2[-(i+1)] = biweight_location(a[-(half_order+i+1):])
    return np.hstack([A1,Ab,A2])
    
    

def biweight_location(a, c=6.0, M=None, axis=None):
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

    a = np.array(a, copy=False)

    if M is None:
        if axis is None:
            M = np.median(a)
        else:   
            M = np.median(a, axis=axis)
            
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
        
    u = d / c / MAD
    
    # now remove the outlier points
    mask = (np.abs(u) < 1).astype(np.int)
    u = (1 - u ** 2) ** 2
    return M + (d * u * mask).sum(axis=axis) / (u * mask).sum(axis=axis)
    
    
def biweight_midvariance(a, c=9.0, M=None, axis=None):
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

    a = np.array(a, copy=False)

    if M is None:
        if axis is None:
            M = np.median(a)
        else:   
            M = np.median(a, axis=axis)

    N = M*1.      
    # set up the difference
    if axis is not None:
        for i in axis:
            N = np.expand_dims(N, axis=i)

    # set up the difference
    d = a - N

    # set up the weighting
    if axis is not None:
        MAD = median_absolute_deviation(a, axis=axis)
        for i in axis:
            MAD = np.expand_dims(MAD, axis=i)
    else:
        MAD = median_absolute_deviation(a)
    # set up the weighting
    u = d / c / MAD

    # now remove the outlier points
    mask = (np.abs(u) < 1).astype(np.int)

    u = u ** 2
    n = mask.sum(axis=axis)
    return n ** 0.5 * ((d*mask)**2 * (1 - u * mask) ** 4).sum(axis=axis) ** 0.5\
        / np.abs(((1 - u * mask) * (1 - 5 * u * mask)).sum(axis=axis))


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