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
    A reduction object 
    :param dim:
        The number of dimensions in the parameter space.

    """
    def __init__(self, args=[], kwargs={}):
        self.args = args
        self.kwargs = kwargs
        
    def model(self):
        model = bias + dark * self.time + fiber_weight * spectrum
        