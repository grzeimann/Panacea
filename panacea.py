#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
IFU Reduction Code 
------------------
Built for the VIRUS instrument as well as LRS2 on HET


"""

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

__all__ = ["Panacea"]

import numpy as np
from .utils import biweight_location, biweight_midvariance


class Panacea(object):
    """
    An abstract sampler object that implements various helper functions
    :param dim:
        The number of dimensions in the parameter space.
    :param lnpostfn:
        A function that takes a vector in the parameter space as input and
        returns the natural logarithm of the posterior probability for that
        position.
    :param args: (optional)
        A list of extra positional arguments for ``lnpostfn``. ``lnpostfn``
        will be called with the sequence ``lnpostfn(p, *args, **kwargs)``.
    :param kwargs: (optional)
        A list of extra keyword arguments for ``lnpostfn``. ``lnpostfn``
        will be called with the sequence ``lnpostfn(p, *args, **kwargs)``.
    """