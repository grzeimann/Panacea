# -*- coding: utf-8 -*-
"""

@author: gregz
"""

import logging
import sys

from astropy import wcs
from numpy import cos, sin, deg2rad

try:
    from pyhetdex.het.fplane import FPlane
    pyhetdex_flag = True
except:
    print('Could not find fplane class in pyhetdex.')
    print('Please check your pyhetdex installation.')
    pyhexdex_flag = False


class Astrometry:
    '''
    Astrometry Class

    This class holds basic astrometric solutions for fits imaging and has a
    suite of functions designed to manage the conversion from on sky
    coordinates to on image coordinates (and vice versa) as well as take
    updates to the solution.  Other functionaility has a specific design
    for serving as the HET's fplane astrometry and building astrometric
    solutions for given ifuslot's based on the fplane astrometry.  This is
    useful for reconstructed VIRUS images.
    '''
    def __init__(self, ra0, dec0, pa, x0, y0, x_scale=-1., y_scale=1.,
                 sys_rot=1.07, fplane_file=None, kind='fplane'):
        self.setup_logging()
        self.ra0 = ra0
        self.dec0 = dec0
        self.dra = 0.
        self.ddec = 0.
        self.dx = 0.
        self.dy = 0.
        self.x0 = x0
        self.y0 = y0
        self.pa = pa
        self.x_scale = x_scale
        self.y_scale = y_scale
        self.sys_rot = sys_rot
        self.fplane_file = fplane_file
        if self.fplane_file is None:
            self.log.info('No fplane file given.')
            self.log.info('Some functions will be unavailable until '
                          'an fplane is given.')
            self.fplane = None
        else:
            if not pyhetdex_flag:
                self.log.warning('Cannot load fplane because pyhetdex '
                                 'is not available')
                self.fplane = None
            else:
                self.fplane = FPlane(self.fplane_file)
        self.kind = kind
        self.set_effective_rotation()

        # Building tangent plane projection with scale 1"
        self.tp = self.setup_TP(self.ra0, self.dec0, self.rot, self.x0,
                                self.y0)
        self.tp_ifuslot = None

    def setup_TP(self, ra0, dec0, rot, x0=0.0, y0=0.0, x_scale=None,
                 y_scale=None):
        ''' TP is tangent plane '''
        ARCSECPERDEG = 1.0/3600.0

        # make a WCS object with appropriate FITS headers
        tp = wcs.WCS(naxis=2)
        tp.wcs.crpix = [x0, y0]  # central "pixel"
        tp.wcs.crval = [ra0, dec0]  # tangent point
        tp.wcs.ctype = ['RA---TAN', 'DEC--TAN']
        # pixel scale in deg.
        if x_scale is None:
            x_scale = self.x_scale
        if y_scale is None:
            y_scale = self.y_scale
        tp.wcs.cdelt = [ARCSECPERDEG * x_scale,
                        ARCSECPERDEG * y_scale]

        # Deal with PA rotation by adding rotation matrix to header
        rrot = deg2rad(rot)
        # clockwise rotation matrix
        tp.wcs.pc = [[cos(rrot), sin(rrot)], [-1.0*sin(rrot), cos(rrot)]]
        return tp

    def set_polynomial_platescale(self):
        ''' This has not been tested '''
        self.tp.wcs.a_0_0 = 0.177311
        self.tp.wcs.a_1_0 = -8.29099e-06
        self.tp.wcs.a_2_0 = -2.37318e-05
        self.tp.wcs.b_0_0 = 0.177311
        self.tp.wcs.b_0_1 = -8.29099e-06
        self.tp.wcs.b_0_2 = -2.37318e-05
        self.tp.wcs.a_order = 2
        self.tp.wcs.b_order = 2

    def setup_logging(self):
        '''Set up a logger for analysis with a name ``astrometry``.

        Use a StreamHandler to write to stdout and set the level to DEBUG if
        verbose is set from the command line
        '''
        self.log = logging.getLogger('astrometry')
        if not len(self.log.handlers):
            fmt = '[%(levelname)s - %(asctime)s] %(message)s'
            level = logging.INFO

            fmt = logging.Formatter(fmt)

            handler = logging.StreamHandler()
            handler.setFormatter(fmt)
            handler.setLevel(level)

            self.log = logging.getLogger('astrometry')
            self.log.setLevel(logging.DEBUG)
            self.log.addHandler(handler)

    def set_effective_rotation(self):
        ''' The rotation for the acam is correct '''
        if self.kind == 'fplane':
            self.rot = 360. - (90. + self.pa + self.sys_rot)
        elif self.kind == 'acam':
            self.rot = self.pa + self.sys_rot + 90.
        elif self.kind == 'lrs2':
            self.rot = 180. - self.pa - self.sys_rot
        else:
            self.log.error('"kind" was not set to available options.')
            self.log.info('Available options are: %s and %s' % ('fplane',
                                                                'acam'))
            self.log.info('Next time please choose one of the options above.')
            self.log.info('Exiting due to error.')
            sys.exit(1)

    def update_projection(self):
        ''' Use this for a new projection with small adjustments '''
        self.set_effective_rotation()
        self.tp = self.setup_TP(self.ra0 + self.dra, self.dec0 + self.ddec,
                                self.rot, self.x0 + self.dx, self.y0 + self.dy)

    def get_ifuslot_ra_dec(self, ifuslot):
        ''' Fplane functionality required for this '''
        if self.fplane is None:
            return None
        ifu = self.fplane.by_ifuslot(ifuslot)
        # remember to flip x,y
        return self.tp.wcs_pix2world(ifu.y, ifu.x, 1)

    def get_ifupos_ra_dec(self, ifuslot, x, y):
        ''' Fplane functionality required for this '''
        if self.fplane is None or not pyhetdex_flag:
            return None
        ifu = self.fplane.by_ifuslot(ifuslot)
        # remember to flip x,y
        return self.tp.wcs_pix2world(ifu.y + x, ifu.x + y, 1)

    def get_ifuslot_projection(self, ifuslot, imscale, crx, cry):
        ''' Fplane functionality required for this '''
        if self.fplane is None or not pyhetdex_flag:
            return None
        ra, dec = self.get_ifuslot_ra_dec(ifuslot)
        self.tp_ifuslot = self.setup_TP(ra, dec, self.rot, crx, cry,
                                        x_scale=-imscale, y_scale=imscale)

    def convert_ifuslot_xy_to_new_xy(self, x, y, wcs):
        ''' Fplane functionality required for this '''
        if self.tp_ifuslot is None or not pyhetdex_flag:
            self.log.error('You have not setup the ifuslot projection yet.')
            self.log.error('To do so, call '
                           '"get_ifuslot_projection(ifuslot, imscale')
            return None
        ra, dec = self.tp_ifuslot.wcs.wcs_pix2world(x, y, 1)
        return wcs.wcs_world2pix(ra, dec, 1)
