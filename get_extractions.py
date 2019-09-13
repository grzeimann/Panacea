#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 30 13:16:47 2019

@author: gregz
"""

import argparse as ap
import astropy.units as u
import numpy as np

from input_utils import setup_logging
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.table import Table, Column
from astropy.table import vstack
from astropy.time import Time
from datetime import datetime as dt
from hetdex_api.extract import Extract 
from hetdex_api.survey import Survey

def get_spectrum(data, error, mask, weights):
        '''
        Weighted spectral extraction
        
        Parameters
        ----------
        data: numpy 2d array (number of fibers by wavelength dimension)
            Flux calibrated spectra for fibers within a given radius of the 
            source.  The radius is defined in get_fiberinfo_for_coord().
        error: numpy 2d array
            Error of the flux calibrated spectra 
        mask: numpy 2d array (bool)
            Mask of good wavelength regions for each fiber
        weights: numpy 2d array (float)
            Normalized extraction model for each fiber
        
        Returns
        -------
        spectrum: numpy 1d array
            Flux calibrated extracted spectrum
        spectrum_error: numpy 1d array
            Error for the flux calibrated extracted spectrum
        '''
        spectrum = (np.sum(data * mask * weights, axis=0) /
                    np.sum(mask * weights**2, axis=0))
        spectrum_error = (np.sqrt(np.sum(error**2 * mask * weights, axis=0)) /
                          np.sqrt(np.sum(mask * weights**2, axis=0)))
        # Only use wavelengths with enough weight to avoid large noise spikes
        w = np.sum(mask * weights, axis=0)
        bad = w < 0.05
        for i in np.arange(1, 4):
            bad[i:] += bad[:-i]
            bad[:-i] += bad[i:]
        spectrum[bad] = np.nan
        spectrum_error[bad] = np.nan
        
        return spectrum, spectrum_error, w

parser = ap.ArgumentParser(add_help=True)

parser.add_argument("filename",
                    help='''Name of fits file with Bin Table''',
                    type=str)

parser.add_argument("outputname",
                    help='''Name of fits file output''',
                    type=str)

args = parser.parse_args(args=None)

log = setup_logging('toy')
log.info('Loading Survey')
survey = Survey('hdr1')
t = Table(survey.hdfile.root.Survey[:])
sel = (t['response_4540'] > 0.08) * (t['fwhm_moffat'] < 2.6)
t = t[sel]
survey.coords = survey.coords[sel]

log.info('Loading External File')

fitsfile = fits.open(args.filename)
bintable = fitsfile[1].data
table = Table(bintable)

ID = bintable['source_id']

coords = SkyCoord(bintable['ra']*u.deg, bintable['dec']*u.deg)

max_sep = 11.0 * u.arcminute

log.info('Finding shots of interest')

matched_sources = {}
shots_of_interest = []
E = Extract()
# Build aperture PSF for aperture extraction
fixed_aperture = 3.

# Using box size of 10.5 (length of box side) and pixel scale of 0.25
# To see documentation use: help(E.tophat_psf)
aperture = E.tophat_psf(fixed_aperture, 10.5, 0.25)

Sources, Spectra, Error, Weights = ([], [], [], [])

for i, coord in enumerate(survey.coords):
    dist = coords.separation(coord)
    sep_constraint = dist < max_sep
    name = '%sv%03d' % (t['date'][i], t['obsid'][i])
    idx = np.where(sep_constraint)[0]
    matched_sources[name] = idx
    if len(idx) > 0:
        shots_of_interest.append([coord, t['date'][i], i])
log.info('Number of shots of interest: %i' % len(shots_of_interest))

N = len(shots_of_interest)

graceful_exit = 0
table.add_column(Column(np.zeros((len(table),), dtype=int), name='obs_id'))
for j, _info in enumerate(shots_of_interest):
    coord = _info[0]
    date = str(_info[1])
    i = _info[2]
    fwhm = t['fwhm_moffat'][i]
    moffat = E.moffat_psf(fwhm, 10.5, 0.25)

    epoch = Time(dt(int(date[:4]), int(date[4:6]), int(date[6:8]))).byear
    try:
        deltaRA = ((epoch - 2015.5) * bintable['pmra'] / 1e3 / 3600. /
                   np.cos(bintable['dec'] * np.pi / 180.))
        deltaDE = (epoch - 2015.5) * bintable['pmdec'] / 1e3 / 3600.
        deltaRA[np.isnan(deltaRA)] = 0.0
        deltaDE[np.isnan(deltaDE)] = 0.0
        ncoords = SkyCoord((bintable['ra']+deltaRA)*u.deg,
                           (bintable['dec']+deltaDE)*u.deg)
    except:
        log.warning("Can't convert proper motion for epoch")
        ncoords = coords
    dist = ncoords.separation(coord)
    sep_constraint = dist < max_sep
    name = '%sv%03d' % (t['date'][i], t['obsid'][i])
    intname = int(name.replace('v', ''))
    idx = np.where(sep_constraint)[0]
    matched_sources[name] = idx
    if len(idx) > 0:
        log.info('Working on shot [%i / %i]: %s' % (j+1, N, name))
        E.load_shot(name)
        for ind in idx:
            info_result = E.get_fiberinfo_for_coord(ncoords[ind], radius=7.)
            if info_result is not None:
                log.info('Extracting %s' % str(ID[ind]))
                ifux, ifuy, xc, yc, ra, dec, data, error, mask = info_result
                weights = E.build_weights(xc, yc, ifux, ifuy, moffat)
                second_mask = np.sqrt((ifux-xc)**2 + (ifuy-yc)**2) < 3.
                result = get_spectrum(data, error,
                                      mask*second_mask[:, np.newaxis],
                                      weights)
                spectrum_aper, spectrum_aper_error, w = [res for res in result]
                dtable = Table(table[ind])
                dtable['obs_id'] = intname
                Sources.append(dtable)
                Spectra.append(spectrum_aper)
                Error.append(spectrum_aper_error)
                Weights.append(w)
        E.fibers.close()

F1 = fits.BinTableHDU(vstack(Sources))
F1.header['EXTNAME'] = 'catalog'
F2 = fits.ImageHDU(np.array(Spectra))
F2.header['EXTNAME'] = 'spectra'
F2.header['CRVAL1'] = 3470.
F2.header['CDELTA1'] = 2.
F3 = fits.ImageHDU(np.array(Error))
F3.header['EXTNAME'] = 'error'
F3.header['CRVAL1'] = 3470.
F3.header['CDELTA1'] = 2.
F4 = fits.ImageHDU(np.array(Weights))
F4.header['EXTNAME'] = 'weight'
F4.header['CRVAL1'] = 3470.
F4.header['CDELTA1'] = 2.
F = fits.HDUList([fits.PrimaryHDU(), F1, F2, F3, F4])
F.writeto(args.outputname, overwrite=True)
