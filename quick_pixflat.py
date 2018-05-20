# -*- coding: utf-8 -*-
"""
Created on Sun May 20 07:39:27 2018

@author: gregz
"""

import numpy as np
import os.path as op
import glob
from distutils.dir_util import mkpath
from input_utils import setup_basic_parser, setup_logging
from amplifier import Amplifier
from astropy.io import fits
from scipy.signal import savgol_filter
from photutils import Background2D, SExtractorBackground
try:
    from photutils import SigmaClip
except:
    from astropy.stats import SigmaClip
from photutils import detect_sources
from astropy.convolution import Gaussian2DKernel
from astropy.stats import gaussian_fwhm_to_sigma


kwargs = {'refit': True, 'use_pixelflat': False, 'dark_mult': 0.0,
          'virusconfig': '/Users/gregz/cure/virus_early/virus_config',
          'bias_mult': 0.0, 'specname': 'virus', 'use_trace_ref': False,
          'check_fibermodel': True, 'fibmodel_slope': 0.001, 'fsize': 6.,
          'fibmodel_intercept': 0.002, 'fibmodel_breakpoint': 4.}

parser = setup_basic_parser()
parser.add_argument("-fx", "--filter_x_size",
                    help='''List of filter sizes in x-direction for each loop''',
                    type=str, default='75, 63, 51, 39')
args = parser.parse_args(args=None)
args.filter_x_size = [int(x) for x in args.filter_x_size.split(',')]
args.log = setup_logging(logname='pixelflat')

AMPS = ['LL', 'LU', 'RL', 'RU']


def measure_image_background(image, mask):
    '''
    image
    '''
    sigma_clip = SigmaClip(sigma=3., iters=5)
    bkg_estimator = SExtractorBackground()
    bkg = Background2D(image, (100, 100), filter_size=(3, 3),
                       sigma_clip=sigma_clip, bkg_estimator=bkg_estimator,
                       mask=mask)
    return image-bkg.background, bkg


def detect_in_image(image_sub, bkg, thresh=3, fwhm=1.5, scale=1.0):
    '''
    Make detections in sky-subtracted image
    image_sub, bkg
    '''
    threshold = (thresh * bkg.background_rms)
    fwhm_i = fwhm / scale
    sigma_i = fwhm_i * gaussian_fwhm_to_sigma
    kernel_size = int(sigma_i*4.)
    kernel = Gaussian2DKernel(sigma_i, x_size=kernel_size, y_size=kernel_size)
    kernel.normalize()
    segm = detect_sources(np.array(image_sub), threshold, npixels=3,
                          filter_kernel=None)
    return segm


def fit_continuum(wv, sky, mask, fil_len=95, func=np.array):
    skym_s = 1. * sky
    skym_s[~mask] = np.interp(wv[~mask], wv[mask], skym_s[mask])
    sky_sm = savgol_filter(skym_s, fil_len, 3)
    return sky_sm


def get_filenames(args, date, obsid, amp, ifuslot='022', specname='virus'):
    name = op.join(args.rootdir, date, specname, specname + '%07d' % int(obsid),
                   'exp*', specname, '*%s%s*.fits' % (ifuslot, amp))
    return sorted(glob.glob(name))


def make_master_amp(args, date, obs, amp, ifuslot='022'):
    amp_list = []
    filenames = get_filenames(args, date, obs, amp, ifuslot)
    for fn in filenames:
        amp_list.append(Amplifier(fn, '', **kwargs))
        amp_list[-1].subtract_overscan()
        amp_list[-1].trim_image()
    avgimg = np.median([a.image for a in amp_list], axis=0)
    amp_list[-1].image = avgimg
    return amp_list[-1]


def main():
    for j, amp in enumerate(AMPS):
        ldls = make_master_amp(args, args.date, args.observation, amp,
                               ifuslot=args.ifuslot)
        ldls.image /= np.nanmedian(ldls.image)
        x = np.arange(ldls.image.shape[1])
        y = np.arange(ldls.image.shape[0])

        mask = np.ones(ldls.image.shape, dtype=bool)
        mask[:, 514:518] = False
        for s, v in enumerate(args.filter_x_size):
            args.log.info('Starting iteration %i for %s%s' % (s+1,
                                                              args.ifuslot,
                                                              amp))
            modelflat_0 = np.zeros(ldls.image.shape)
            modelflat_1 = np.ones(ldls.image.shape)
            for k in np.arange(ldls.image.shape[0]):
                modelflat_0[k, :] = fit_continuum(x, ldls.image[k, :],
                                                  mask[k, :], fil_len=v)
            temp = ldls.image / modelflat_0
            for k in np.arange(ldls.image.shape[0]):
                if mask[:, k].sum() > (len(x) / 2):
                    modelflat_1[:, k] = fit_continuum(y, temp[:, k],
                                                      mask[:, k], fil_len=101)
            pixelflat = ldls.image / modelflat_0 / modelflat_1
            pixelflat[~np.isfinite(pixelflat)] = 0.0
            bad = (((pixelflat) < 0.97) + ((pixelflat) > 1.03))
            mask[bad] = False
            image_backsub, bkg = measure_image_background(pixelflat, ~mask)
            segm = detect_in_image(image_backsub, bkg)
            segm1 = detect_in_image(-1.*image_backsub, bkg)
            bad = ((segm.array > 0) + (segm1.array > 0))
            mask[bad] = False
            args.log.info('Number of bad pixels: %i' % ((~mask).sum()))
        ldls.pixelflat = pixelflat
        ldls.temp = temp
        ldls.mask = mask
        hdu = fits.PrimaryHDU(np.array(ldls.pixelflat, dtype='float32'))
        mkpath(op.join('pixelflat', args.date))
        hdu.writeto('pixelflat/%s/pixelflat_cam%s_%s.fits' %
                    (args.date, ldls.specid, amp), overwrite=True)

if __name__ == '__main__':
    main()
