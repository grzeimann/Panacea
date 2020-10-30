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
from astropy.modeling.models import Moffat2D
from astropy.table import Table, Column
from astropy.table import vstack
from astropy.time import Time
from datetime import datetime as dt
from hetdex_api.extract import Extract 
from hetdex_api.survey import Survey
from photutils.centroids import centroid_2dg
from scipy.interpolate import LinearNDInterpolator


def build_weights(self, xc, yc, ifux, ifuy, psf):
    '''
    Build weight matrix for spectral extraction
    
    Parameters
    ----------
    xc: float
        The ifu x-coordinate for the center of the collapse frame
    yc: float 
        The ifu y-coordinate for the center of the collapse frame
    xloc: numpy array
        The ifu x-coordinate for each fiber
    yloc: numpy array
        The ifu y-coordinate for each fiber
    psf: numpy 3d array
        zeroth dimension: psf image, xgrid, ygrid
    
    Returns
    -------
    weights: numpy 2d array (len of fibers by wavelength dimension)
        Weights for each fiber as function of wavelength for extraction
    '''
    S = np.zeros((len(ifux), 2))
    T = np.array([psf[1].ravel(), psf[2].ravel()]).swapaxes(0, 1)
    I = LinearNDInterpolator(T, psf[0].ravel(),
                             fill_value=0.0)
    weights = np.zeros((len(ifux), len(self.wave)))
    for i in np.arange(len(self.wave)):
        S[:, 0] = ifux - self.ADRx[i] - xc
        S[:, 1] = ifuy - self.ADRy[i] - yc
        weights[:, i] = I(S[:, 0], S[:, 1])
    self.log.info('Average weight is: %0.2f' % np.median(np.sum(weights, axis=0)))
    return weights

def get_fiberinfo_for_coord(self, coord, radius=8.0, ffsky=False):
    """ 
    Grab fibers within a radius and get relevant info
    
    Parameters
    ----------
    coord: SkyCoord Object
        a single SkyCoord object for a given ra and dec
    radius:
        radius to extract fibers in arcsec
    ffsky: bool
        Flag to choose local (ffsky=False) or full frame (ffsky=True)
        sky subtraction
    
    Returns
    -------
    ifux: numpy array (length of number of fibers)
        ifu x-coordinate accounting for dither_pattern
    ifuy: numpy array (length of number of fibers)
        ifu y-coordinate accounting for dither_pattern
    ra: numpy array (length of number of fibers)
        Right ascension of fibers
    dec: numpy array (length of number of fibers)
        Declination of fibers
    spec: numpy 2d array (number of fibers by wavelength dimension)
        Calibrated spectra for each fiber
    spece: numpy 2d array (number of fibers by wavelength dimension)
        Error for calibrated spectra
    mask: numpy 2d array (number of fibers by wavelength dimension)
        Mask of good values for each fiber and wavelength
    """
    
    fiber_lower_limit = 7
    
    idx = self.fibers.query_region_idx(coord, radius=radius)

    if len(idx) < fiber_lower_limit:
        return None

    ifux = self.fibers.table.read_coordinates(idx, "ifux")
    ifuy = self.fibers.table.read_coordinates(idx, "ifuy")
    ra = self.fibers.table.read_coordinates(idx, "ra")
    dec = self.fibers.table.read_coordinates(idx, "dec")
    mname = [x.decode("utf-8") 
             for x in self.fibers.table.read_coordinates(idx, "multiframe")]
    if ffsky:
        spec = self.fibers.table.read_coordinates(idx, "spec_fullsky_sub") / 2.0
    else:
        spec = self.fibers.table.read_coordinates(idx, "calfib") / 2.0

    spece = self.fibers.table.read_coordinates(idx, "calfibe") / 2.0
    ftf = self.fibers.table.read_coordinates(idx, "fiber_to_fiber")
    
    if self.survey == 'hdr1':
        mask = self.fibers.table.read_coordinates(idx, "Amp2Amp")
        mask = (mask > 1e-8) * (np.median(ftf, axis=1) > 0.5)[:, np.newaxis]
    else:
        mask = self.fibers.table.read_coordinates(idx, "calfibe")
        mask = (mask > 1e-8) * (np.median(ftf, axis=1) > 0.5)[:, np.newaxis]
    expn = np.array(
        self.fibers.table.read_coordinates(idx, "expnum"), dtype=int
    )
    
    ifux[:] = ifux + self.dither_pattern[expn - 1, 0]
    ifuy[:] = ifuy + self.dither_pattern[expn - 1, 1]
    xc, yc = self.convert_radec_to_ifux_ifuy(
        ifux, ifuy, ra, dec, coord.ra.deg, coord.dec.deg
    )
    return ifux, ifuy, xc, yc, ra, dec, spec, spece, mask, mname

def moffat_psf(seeing, boxsize, scale, alpha=3.5):
    '''
    Moffat PSF profile image
    
    Parameters
    ----------
    seeing: float
        FWHM of the Moffat profile
    boxsize: float
        Size of image on a side for Moffat profile
    scale: float
        Pixel scale for image
    alpha: float
        Power index in Moffat profile function
    
    Returns
    -------
    zarray: numpy 3d array
        An array with length 3 for the first axis: PSF image, xgrid, ygrid
    '''
    M = Moffat2D()
    M.alpha.value = alpha
    M.gamma.value = 0.5 * seeing / np.sqrt(2**(1./ M.alpha.value) - 1.)
    xl, xh = (0. - boxsize / 2., 0. + boxsize / 2. + scale)
    yl, yh = (0. - boxsize / 2., 0. + boxsize / 2. + scale)
    x, y = (np.arange(xl, xh, scale), np.arange(yl, yh, scale))
    xgrid, ygrid = np.meshgrid(x, y)
    Z = moffat_psf_integration(xgrid.ravel(), ygrid.ravel(),
                                    seeing, boxsize=boxsize+1.5,
                                    alpha=alpha)
    Z = np.reshape(Z, xgrid.shape)
    zarray = np.array([Z, xgrid, ygrid])
    return zarray

def moffat_psf_integration(xloc, yloc, seeing, boxsize=14.,
                           scale=0.1, alpha=3.5):
    '''
    Moffat PSF profile image
    
    Parameters
    ----------
    seeing: float
        FWHM of the Moffat profile
    boxsize: float
        Size of image on a side for Moffat profile
    scale: float
        Pixel scale for image
    alpha: float
        Power index in Moffat profile function
    
    Returns
    -------
    zarray: numpy 3d array
        An array with length 3 for the first axis: PSF image, xgrid, ygrid
    '''
    M = Moffat2D()
    M.alpha.value = alpha
    M.gamma.value = 0.5 * seeing / np.sqrt(2**(1./ M.alpha.value) - 1.)
    xl, xh = (0. - boxsize / 2., 0. + boxsize / 2. + scale)
    yl, yh = (0. - boxsize / 2., 0. + boxsize / 2. + scale)
    x, y = (np.arange(xl, xh, scale), np.arange(yl, yh, scale))
    xgrid, ygrid = np.meshgrid(x, y)
    Z = M(xgrid, ygrid)
    Z = Z / Z.sum()
    V = xloc * 0.
    cnt = 0
    for xl, yl in zip(xloc, yloc):
        d = np.sqrt((xgrid-xl)**2 + (ygrid-yl)**2)
        sel = d <= 0.75
        adj = np.pi * 0.75**2 / (sel.sum() * scale**2)
        V[cnt] = np.sum(Z[sel]) * adj
        cnt += 1
    return V

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

parser.add_argument("-s", "--survey", type=str,
		    help='''survey name; hdrX''',
		    default='hdr1')

parser.add_argument("-ra", "--RA", type=str,
		    help='''ra name''',
		    default='ra')

parser.add_argument("-dec", "--Dec", type=str,
		    help='''dec name''',
		    default='dec')

parser.add_argument("-r", "--recenter",
                    help='''Re-centroid source''',
                    action="count", default=0)

parser.add_argument("-iba", "--ignore_bad_amps",
                    help='''Still use bad amplifiers ignoring that they are bad''',
                    action="count", default=0)

parser.add_argument("-rv", "--recenter_var2",
                    help='''Re-centroid source variation''',
                    action="count", default=0)

args = parser.parse_args(args=None)

if args.survey == 'hdr1':
    fname = 'fwhm_moffat'
else:
    fname = 'fwhm_virus'
log = setup_logging('toy')
log.info('Loading Survey')
survey = Survey(args.survey)
t = Table(survey.hdfile.root.Survey[:])
sel = (t['response_4540'] > 0.01) * (t[fname] < 3.5)
t = t[sel]
survey.coords = survey.coords[sel]

log.info('Loading External File')

fitsfile = fits.open(args.filename)
bintable = fitsfile[1].data
table = Table(bintable)

# Bad amplifier information
badamprec = fits.open('/scratch/03946/hetdex/hdr2.1/survey/amp_flag.fits')[1].data

ID = bintable['source_id']

coords = SkyCoord(bintable[args.RA]*u.deg, bintable[args.Dec]*u.deg)

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
if args.recenter:
    Images = []

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
    fwhm = t[fname][i]
    moffat = moffat_psf(fwhm, 10.5, 0.25)

    epoch = Time(dt(int(date[:4]), int(date[4:6]), int(date[6:8]))).byear
    try:
        deltaRA = ((epoch - 2015.5) * bintable['pmra'] / 1e3 / 3600. /
                   np.cos(bintable[args.Dec] * np.pi / 180.))
        deltaDE = (epoch - 2015.5) * bintable['pmdec'] / 1e3 / 3600.
        deltaRA[np.isnan(deltaRA)] = 0.0
        deltaDE[np.isnan(deltaDE)] = 0.0
        ncoords = SkyCoord((bintable[args.RA]+deltaRA)*u.deg,
                           (bintable[args.Dec]+deltaDE)*u.deg)
    except:
        log.warning("Can't convert proper motion for epoch")
        ncoords = coords
    dist = ncoords.separation(coord)
    sep_constraint = dist < max_sep
    name = '%sv%03d' % (t['date'][i], t['obsid'][i])
    shotid = int('%s%03d' % (t['date'][i], t['obsid'][i]))
    intname = int(name.replace('v', ''))
    idx = np.where(sep_constraint)[0]
    matched_sources[name] = idx
    if len(idx) > 0:
        log.info('Working on shot [%i / %i]: %s' % (j+1, N, name))
        E.load_shot(name, survey=args.survey)
        ampinds = np.where(badamprec['shotid'] == shotid)[0]
        ampnames = badamprec['multiframe'][ampinds]
        ampflags = badamprec['flag'][ampinds]

        for ind in idx:
            info_result = get_fiberinfo_for_coord(E, ncoords[ind], radius=10.)
            if info_result is not None:
                log.info('Extracting %s' % str(ID[ind]))
                ifux, ifuy, xc, yc, ra, dec, data, error, mask, mname = info_result
                any_in_bad_amp = False
                for jiter in np.arange(len(mname)):
                    s = mname[jiter] == ampnames
                    if ampflags[s] < 1:
                        if not args.ignore_bad_amps:
                            mask[jiter] = False
                        any_in_bad_amp = True
                if any_in_bad_amp:
                    log.info('Some Fibers in Bad Amplifier')
                if args.recenter:
                    wl = 5001. * (1 + bintable['z'][ind])
                    wh = 5013. * (1 + bintable['z'][ind])
                    log.info('Collapsing around %0.2f-%0.2f' % (wl, wh))
                    try:
                        dra = (np.cos(np.deg2rad(np.median(dec))) * 3600. *
                               (ra - ncoords[ind].ra.deg))
                        ddec = (3600. * (dec - ncoords[ind].dec.deg))
                        zarray1 = E.make_collapsed_image(0., 0., dra, ddec, data, mask,
                                                      scale=0.25, seeing_fac=1.5, boxsize=14.,
                                                      wrange=[wl, wh], nchunks=1,
                                                      convolve_image=True,
                                                      interp_kind='linear')
                        if args.recenter_var2:
                            nx, ny = centroid_2dg(zarray1[0])
                            nxc = np.interp(nx, np.arange(zarray1[1].shape[1]),
                                            zarray1[1][0, :])
                            nyc = np.interp(ny, np.arange(zarray1[2].shape[0]),
                                            zarray1[2][:, 0])
                            for n in [nxc, nyc]:
                                if np.isnan(n):
                                    n = 0.0
                            nra = ncoords[ind].ra.deg + nxc / np.cos(np.deg2rad(np.median(dec))) / 3600.
                            ndec = ncoords[ind].dec.deg + nyc / 3600.
                            log.info('%s: Shift: %0.2f, %0.2f, New: %0.6f, %0.5f' %
                                     (str(ID[ind]), nxc, nyc, nra, ndec))
                            zarray = E.make_collapsed_image(xc, yc, ifux, ifuy, data, mask,
                                                          scale=0.25, seeing_fac=1.5, boxsize=10.,
                                                          wrange=[wl, wh], nchunks=1,
                                                          convolve_image=True,
                                                          interp_kind='linear')
                            nx, ny = centroid_2dg(zarray[0])
                            nxc = np.interp(nx, np.arange(zarray[1].shape[1]),
                                            zarray[1][0, :])
                            nyc = np.interp(ny, np.arange(zarray[2].shape[0]),
                                            zarray[2][:, 0])
                            for n in [nxc, nyc]:
                                if np.isnan(n):
                                    n = 0.0
                            log.info('Original: %0.2f, %0.2f, Change: %0.2f, %0.2f' %
                                     (xc, yc, nxc, nyc))
                            xc, yc = (nxc+xc, nyc+yc)
                    except:
                        log.warning('Image Collapse Failed')
                        N1 = int(14. / 0.25)
                        zarray1 = [np.zeros((N1, N1)), np.zeros((N1, N1)),
                                  np.zeros((N1, N1))]
                    Images.append(zarray1[0])
                weights = build_weights(E, xc, yc, ifux, ifuy, moffat)
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
if args.recenter:
    F5 = fits.ImageHDU(np.array(Images))
    F = fits.HDUList([fits.PrimaryHDU(), F1, F2, F3, F4, F5])
else:
    F = fits.HDUList([fits.PrimaryHDU(), F1, F2, F3, F4])
F.writeto(args.outputname, overwrite=True)
