#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Shot Class
-----------
To be used in conjuction with IFU reduction code, Panacea


"""

import numpy as np
from pyhetdex.het.ifu_centers import IFUCenter
from astropy.modeling.models import Moffat2D, Gaussian2D
from photutils import CircularAperture, aperture_photometry
from astropy.coordinates import SkyCoord
from astropy.coordinates import match_coordinates_sky
from astropy.nddata import Cutout2D
from astropy.wcs import WCS
from astropy.wcs.utils import skycoord_to_pixel
from astropy import units as u
from pyhetdex.het.fplane import FPlane
from pyhetdex.coordinates.tangent_projection_astropy import TangentPlane as TP

