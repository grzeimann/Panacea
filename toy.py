#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 30 13:16:47 2019

@author: gregz
"""

import sys
import numpy as np
import astropy.units as u

from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.table import Table
from hetdex_api.extract import Extract 
from hetdex_api.survey import Survey


survey = Survey('hdr1')
t = Table(survey.hdfile.root.Survey[:])

filename = '/work/03730/gregz/maverick/MUSE_WIDE_sources_for_hetdex.fits'

fitsfile = fits.open(filename)
bintable = fitsfile[1].data

ID = bintable['ID']

coords = SkyCoord(bintable['RA']*u.deg, bintable['DEC']*u.deg)
max_sep = 11.0 * u.arcminute

matched_sources = {}
shots_of_interest = []
for i, coord in enumerate(survey.coords):
    dist = coords.separation(coord)
    sep_constraint = dist < max_sep
    name = '%sv%03d' % (t['date'][i], t['obsid'][i])
    idx = np.where(sep_constraint)[0]
    matched_sources[name] = idx
    if len(idx) > 0:
        shots_of_interest.append(name)
print(shots_of_interest)