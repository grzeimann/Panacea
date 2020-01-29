"""
Created on Wed Oct  5 18:12:13 2016

@author: gregz


"""

import numpy as np

    
def biweight(a, nan_treatment=True, c=9.0, axis=None, calc_std=False):
    """
    Returns the biweight location for the array elements along
    an axis if given.
    
    For more details, see Beers, Flynn, and Gebhardt, 1990, AJ, 100, 32B

    Parameters
    ----------
    a: ndarray
        Input array or object that can be converted to an array
    nan_treatment: boolean
        If true, nan's are ignored in the computation.
    c: float
        Tuning constant for the biweight estimator
    axis : integer or tuple
        Axis along which the biweight_location values are computed. 
        The default (axis=None) is to compute along a flattened version
        of the array.

    Returns
    -------
    biweight_location: float or array
        Returns the biweight location for the array elements.
    biweight_std: float or array
        Returns the biweight std for the array elements.
    """
    if axis is not None:
        if isinstance(axis, int):
            axis = (axis,)
        else:
            axis = tuple(sorted(axis))
            
    a = np.asanyarray(a)
    
    if nan_treatment:
        func = np.nanmedian
        sumfunc = np.nansum
    else: 
        func = np.median
        sumfunc = np.sum

    a_median = func(a, axis=axis)

    # re-broadcast the output median array to subtract it
    if axis is not None:
        for i in axis:
            a_median_exp = np.expand_dims(a_median, axis=i)
    else:
        a_median_exp = a_median
    
    # calculated the median average deviation
    d = a - a_median_exp
    MAD = func(np.abs(d), axis=axis)
    if axis is not None:
        for i in axis:
            MAD = np.expand_dims(MAD, axis=i)
    
    u = d / c / MAD
     
    mask = ((np.abs(u) < 1.) * np.isfinite(u)).astype(np.int)
    u1 = (1 - u**2)**2
    BL = (a_median + sumfunc(d * u1 * mask, axis=axis) /
          sumfunc(u1 * mask, axis))
    if type(BL) is np.ndarray:
        badpix = np.isnan(BL)
        BL[badpix] = a_median[badpix]
    else:
        if np.isnan(BL):
            BL = a_median
    if calc_std:
         u = u ** 2
         n = sumfunc(mask, axis=axis)
         BM = (n**0.5 * sumfunc(mask * d**2 * (1. - u)**4, axis=axis)**0.5 /
               np.abs(sumfunc(mask * (1. - u) * (1. - 5. * u), axis=axis)))
         return BL, BM
    else:
        return BL