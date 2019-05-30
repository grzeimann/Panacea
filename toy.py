#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 30 13:16:47 2019

@author: gregz
"""

import astropy.units as u
import numpy as np
import pickle

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

coords = SkyCoord(bintable['RA']*u.deg, bintable['DEC']*u.deg)

max_sep = 11.0 * u.arcminute

log.info('Finding shots of interest')

matched_sources = {}
shots_of_interest = []
E = Extract()
# Build aperture PSF for aperture extraction
fixed_aperture = 4.

# Using box size of 10.5 (length of box side) and pixel scale of 0.25
# To see documentation use: help(E.tophat_psf)
aperture = E.tophat_psf(fixed_aperture, 10.5, 0.25)

Sources = {}
for i in ID:
    Sources[i] = []
for i, coord in enumerate(survey.coords):
    dist = coords.separation(coord)
    sep_constraint = dist < max_sep
    name = '%sv%03d' % (t['date'][i], t['obsid'][i])
    idx = np.where(sep_constraint)[0]
    matched_sources[name] = idx
    if len(idx) > 0:
        shots_of_interest.append(name)
log.info('Number of shots of interest: %i' % len(shots_of_interest))

for i, coord in enumerate(survey.coords):
    dist = coords.separation(coord)
    sep_constraint = dist < max_sep
    name = '%sv%03d' % (t['date'][i], t['obsid'][i])
    idx = np.where(sep_constraint)[0]
    matched_sources[name] = idx
    if len(idx) > 0:
        log.info('Working on shot: %s' % name)
        E.load_shot(name)
        for ind in idx:
            info_result = E.get_fiberinfo_for_coord(coords[ind], radius=7.)
            if info_result is not None:
                log.info('Extracting %i' % ID[ind])
                ifux, ifuy, xc, yc, ra, dec, data, error, mask = info_result
                weights = E.build_weights(xc, yc, ifux, ifuy, aperture)
                result = E.get_spectrum(data, error, mask, weights)
                spectrum_aper, spectrum_aper_error = [res for res in result]
                Sources[ID[ind]].append([spectrum_aper, spectrum_aper_error,
                                         weights.sum(axis=0)])
        E.fibers.close()
pickle.dump(Sources, open( "save.p", "wb" ))
log.info('Done.')