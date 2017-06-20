# -*- coding: utf-8 -*-
"""
Created on Fri Jun  2 09:19:46 2017

@author: gregz
"""

import numpy as np
from pyhetdex.het.fplane import FPlane
from pyhetdex.coordinates.tangent_projection import TangentPlane as TP
import logging
import sys

class Astrometry:
    def __init__(self, ra, dec, pa, sys_rot=1.3, fplane_file=None):
        self.setup_logging()
        self.ra = ra
        self.dec = dec
        self.dra = 0.
        self.ddec = 0.
        self.pa = pa
        self.sys_rot = sys_rot
        self.fplane_file = fplane_file
        if self.fplane_file is None:
            self.log.error('No fplane file given.')
            self.log.error('Please provide the path to the fplane file.')
            sys.exit(1)
        else:
            self.fplane = FPlane(self.fplane_file)
        self.set_effective_rotation()
        
        # Building tangent plane projection with scale 1"
        self.tp = TP(self.ra, self.dec, self.rot)
        self.tp_ifuslot = None
 

    def set_polynomial_platescale(self):
        self.tp.wcs.a_0_0 = 0.177311
        self.tp.wcs.a_1_0 = -8.29099e-06
        self.tp.wcs.a_2_0 = -2.37318e-05
        self.tp.wcs.b_0_0 = 0.177311
        self.tp.wcs.b_0_1 = -8.29099e-06
        self.tp.wcs.b_0_2 = -2.37318e-05
        self.tp.wcs.a_order=2
        self.tp.wcs.b_order=2
        
        
    def setup_logging(self):
        '''Set up a logger for analysis with a name ``shot``.
    
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
        # Making rotation from the PA
        # TODO: re-derive this formula 
        self.rot = 360. - (90. + self.pa + self.sys_rot)

    def update_projection(self):
        self.set_effective_rotation()
        # Building tangent plane projection with scale 1"
        #dra = self.dra / 3600. / np.cos(np.deg2rad(self.dec))
        #ddec = self.ddec / 3600.
        self.tp = TP(self.ra+self.dra, self.dec+self.ddec, self.rot)
            
    def get_ifuslot_ra_dec(self, ifuslot):
        ifu = self.fplane.by_ifuslot(ifuslot)
        # remember to flip x,y 
        return self.tp.xy2raDec(ifu.y, ifu.x)

    def get_ifuspos_ra_dec(self, ifuslot, x, y):
        ifu = self.fplane.by_ifuslot(ifuslot)
        # remember to flip x,y 
        return self.tp.xy2raDec(ifu.y + x, ifu.x + y)
        
    def get_ifuslot_projection(self, ifuslot, imscale, crx, cry):
        ra, dec = self.get_ifuslot_ra_dec(ifuslot)
        self.tp_ifuslot = TP(ra, dec, self.rot, -imscale, imscale)
        self.tp_ifuslot.wcs.wcs.crpix = [crx, cry]

        
    def convert_ifuslot_xy_to_new_xy(self, x, y, wcs):
        if self.tp_ifuslot is None:
            self.log.error('You have not setup the ifuslot projection yet.')
            self.log.error('To do so, call '
                           '"get_ifuslot_projection(ifuslot, imscale')
            return None
        ra, dec = self.tp_ifuslot.wcs.wcs_pix2world(x, y, 1)
        return wcs.wcs_world2pix(ra, dec, 1)
            
