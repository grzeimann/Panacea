# -*- coding: utf-8 -*-
"""
Created on Mon Feb 12 20:59:53 2018

@author: gregz
"""

from astropy.modeling import Fittable2DModel, Parameter
from astropy.units import UnitsError
import numpy as np
from collections import OrderedDict


class AsymMoffat2D(Fittable2DModel):
    """
    Two dimensional Moffat model with asymmetry.

    Parameters
    ----------
    amplitude : float
        Amplitude of the model.
    x_0 : float
        x position of the maximum of the Moffat model.
    y_0 : float
        y position of the maximum of the Moffat model.
    gamma : float
        Core width of the Moffat model.
    alpha : float
        Power index of the Moffat model.
    ratio : float
        Ratio of y / x axis

    See Also
    --------
    Gaussian2D, Box2D, Moffat2D from astropy.modeling

    Notes
    -----
    Model formula:

    .. math::

        f(x, y) = A \\left(1 + \\frac{\\left(x - x_{0}\\right)^{2} +
        \\left(y - y_{0}\\right)^{2}}{\\gamma^{2}}\\right)^{- \\alpha}
    """

    amplitude = Parameter(default=1)
    x_0 = Parameter(default=0)
    y_0 = Parameter(default=0)
    gamma = Parameter(default=1)
    alpha = Parameter(default=1)
    ratio = Parameter(default=1)

    @property
    def fwhm(self):
        """
        Moffat full width at half maximum.
        Derivation of the formula is available in
        """
        return 2.0 * self.gamma * np.sqrt(2.0 ** (1.0 / self.alpha) - 1.0)

    @staticmethod
    def evaluate(x, y, amplitude, x_0, y_0, gamma, alpha, ratio):
        """ Two dimensional Moffat model function """
        rr_gg = ((ratio * (x - x_0)) ** 2 + (y - y_0) ** 2) / gamma ** 2
        return amplitude * (1 + rr_gg) ** (-alpha)

    @staticmethod
    def fit_deriv(x, y, amplitude, x_0, y_0, gamma, alpha, ratio):
        """ Two dimensional Moffat model derivative """

        rr_gg = ((ratio * (x - x_0)) ** 2 + (y - y_0) ** 2) / gamma ** 2
        d_A = (1 + rr_gg) ** (-alpha)
        d_x_0 = (-amplitude * alpha * d_A * ratio**2 * (-2 * x + 2 * x_0) /
                 (gamma ** 2 * (1 + rr_gg)))
        d_y_0 = (-amplitude * alpha * d_A * (-2 * y + 2 * y_0) /
                 (gamma ** 2 * (1 + rr_gg)))
        d_alpha = -amplitude * d_A * np.log(1 + rr_gg)
        d_gamma = 2 * amplitude * alpha * d_A * (rr_gg / (gamma * (1 + rr_gg)))
        d_ratio = (-amplitude * alpha * d_A * 2 * ratio * (x - x_0)**2 /
                   (gamma ** 2 * (1 + rr_gg)))
        return [d_A, d_x_0, d_y_0, d_gamma, d_alpha, d_ratio]

    @property
    def input_units(self):
        if self.x_0.unit is None:
            return None
        else:
            return {'x': self.x_0.unit,
                    'y': self.y_0.unit}

    def _parameter_units_for_data_units(self, inputs_unit, outputs_unit):
        # Note that here we need to make sure that x and y are in the same
        # units otherwise this can lead to issues since rotation is not well
        # defined.
        if inputs_unit['x'] != inputs_unit['y']:
            raise UnitsError("Units of 'x' and 'y' inputs should match")
        return OrderedDict([('x_0', inputs_unit['x']),
                            ('y_0', inputs_unit['x']),
                            ('gamma', inputs_unit['x']),
                            ('amplitude', outputs_unit['z'])])
