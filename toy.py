#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 30 13:16:47 2019

@author: gregz
"""

import astropy.units as u
import numpy as np
import sys
from input_utils import setup_logging
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.table import Table
from hetdex_api.extract import Extract 
from hetdex_api.survey import Survey

log = setup_logging('toy')
log.info('Loading Survey')
survey = Survey('hdr1')
t = Table(survey.hdfile.root.Survey[:])

log.info('Loading External File')
filename = '/work/03730/gregz/maverick/MUSE_WIDE_sources_for_hetdex.fits'

fitsfile = fits.open(filename)
bintable = fitsfile[1].data

ID = bintable['ID']

s = np.random.rand(10000)*360.
t = np.random.rand(10000)*180. - 90.
coords = SkyCoord(bintable['RA']*u.deg, bintable['DEC']*u.deg)
coords = SkyCoord(s*u.deg, t*u.deg)

max_sep = 11.0 * u.arcminute

log.info('Finding shots of interest')

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
log.info('Done.')
print(shots_of_interest)