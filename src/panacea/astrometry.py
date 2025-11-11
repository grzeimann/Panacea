# -*- coding: utf-8 -*-
"""

@author: gregz
"""

import logging
import sys

from astropy import wcs
from numpy import cos, sin, deg2rad
import numpy as np


class Astrometry:
    """Astrometric transformations and HET/VIRUS focal-plane utilities.

    Provides small utilities to convert between sky coordinates and image
    coordinates using a simple TAN (gnomonic) projection with configurable
    pixel scales and rotation. When an fplane file is supplied, also supports
    mapping IFU slot positions and generating per-IFU WCS projections useful
    for reconstructed VIRUS images.
    """
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
            self.fplane = FPlane(self.fplane_file)
        self.kind = kind
        self.set_effective_rotation()

        # Building tangent plane projection with scale 1"
        self.tp = self.setup_TP(self.ra0, self.dec0, self.rot, self.x0,
                                self.y0)
        self.tp_ifuslot = None

    def setup_TP(self, ra0, dec0, rot, x0=0.0, y0=0.0, x_scale=None,
                 y_scale=None):
        """Create a simple TAN WCS for a given tangent point, scale, and rotation.

        Args:
            ra0: Right ascension of the tangent point in degrees.
            dec0: Declination of the tangent point in degrees.
            rot: Rotation angle in degrees (counterclockwise positive).
            x0: Reference pixel x-coordinate. Defaults to 0.0.
            y0: Reference pixel y-coordinate. Defaults to 0.0.
            x_scale: Pixel scale in arcsec/pixel along X. Defaults to self.x_scale.
            y_scale: Pixel scale in arcsec/pixel along Y. Defaults to self.y_scale.

        Returns:
            astropy.wcs.WCS: Configured WCS object.
        """
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
        """Attach SIP-like polynomial plate scale terms to the WCS.

        Notes:
            Experimental: retained from legacy code and not thoroughly tested.
            Sets a 2nd-order distortion model on the current ``self.tp.wcs``.
        """
        self.tp.wcs.a_0_0 = 0.177311
        self.tp.wcs.a_1_0 = -8.29099e-06
        self.tp.wcs.a_2_0 = -2.37318e-05
        self.tp.wcs.b_0_0 = 0.177311
        self.tp.wcs.b_0_1 = -8.29099e-06
        self.tp.wcs.b_0_2 = -2.37318e-05
        self.tp.wcs.a_order = 2
        self.tp.wcs.b_order = 2

    def setup_logging(self):
        """Create a module logger named 'astrometry' writing to stdout.

        The handler is installed only once per process. The stream handler is
        configured at INFO level with a simple timestamped format, while the
        logger itself is set to DEBUG so downstream code can adjust verbosity by
        changing the logger level.
        """
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
        """Derive effective rotation angle based on instrument kind.

        Sets ``self.rot`` according to the provided ``kind``:
        - 'fplane': 360 - (90 + PA + sys_rot)
        - 'acam':   PA + sys_rot + 90
        - 'lrs2':   180 - PA - sys_rot

        Raises:
            SystemExit: If ``kind`` is not one of {'fplane','acam','lrs2'}.
        """
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
        """Rebuild the TAN WCS after applying small RA/Dec and XY offsets.

        Uses current values of dra/ddec/dx/dy together with PA/system rotation
        to recompute the internal ``self.tp`` WCS.
        """
        self.set_effective_rotation()
        self.tp = self.setup_TP(self.ra0 + self.dra, self.dec0 + self.ddec,
                                self.rot, self.x0 + self.dx, self.y0 + self.dy)

    def get_ifuslot_ra_dec(self, ifuslot):
        """Return sky coordinates (RA, Dec) of the IFU slot center.

        Requires an fplane to be loaded at construction time.

        Args:
            ifuslot: IFU slot identifier string.

        Returns:
            tuple: (ra_deg, dec_deg) in degrees, or None if no fplane available.
        """
        if self.fplane is None:
            return None
        ifu = self.fplane.by_ifuslot(ifuslot)
        # remember to flip x,y
        return self.tp.wcs_pix2world(ifu.y, ifu.x, 1)

    def get_ifupos_ra_dec(self, ifuslot, x, y):
        """Return sky coordinates (RA, Dec) for a position within an IFU.

        Requires an fplane to be loaded at construction time.

        Args:
            ifuslot: IFU slot identifier string.
            x: X offset from IFU center in IFU pixel units.
            y: Y offset from IFU center in IFU pixel units.

        Returns:
            tuple: (ra_deg, dec_deg) in degrees, or None if no fplane available.
        """
        if self.fplane is None:
            return None
        ifu = self.fplane.by_ifuslot(ifuslot)
        # remember to flip x,y
        return self.tp.wcs_pix2world(ifu.y + x, ifu.x + y, 1)

    def get_ifuslot_projection(self, ifuslot, imscale, crx, cry):
        """Create a TAN WCS centered on an IFU slot with given image scale.

        Args:
            ifuslot: IFU slot identifier string.
            imscale: Pixel scale (arcsec/pixel); applied as x_scale=-imscale, y_scale=imscale.
            crx: Reference pixel x for the new WCS.
            cry: Reference pixel y for the new WCS.

        Returns:
            None: Sets ``self.tp_ifuslot`` to the newly constructed WCS, or None if
            no fplane is available.
        """
        if self.fplane is None:
            return None
        ra, dec = self.get_ifuslot_ra_dec(ifuslot)
        self.tp_ifuslot = self.setup_TP(ra, dec, self.rot, crx, cry,
                                        x_scale=-imscale, y_scale=imscale)

    def convert_ifuslot_xy_to_new_xy(self, x, y, wcs):
        """Transform pixel coordinates from IFU-centered WCS to another WCS.

        Args:
            x: X coordinate(s) in the IFU-centered WCS ``self.tp_ifuslot``.
            y: Y coordinate(s) in the IFU-centered WCS ``self.tp_ifuslot``.
            wcs: Target astropy.wcs.WCS to transform into.

        Returns:
            tuple: (x_new, y_new) in the target WCS pixel coordinates, or None if
            ``self.tp_ifuslot`` has not been set via ``get_ifuslot_projection``.
        """
        if self.tp_ifuslot is None:
            self.log.error('You have not setup the ifuslot projection yet.')
            self.log.error('To do so, call "get_ifuslot_projection(ifuslot, imscale)"')
            return None
        ra, dec = self.tp_ifuslot.wcs.wcs_pix2world(x, y, 1)
        return wcs.wcs_world2pix(ra, dec, 1)


class IFU(object):
    """
    Class representing an Integral Field Unit (IFU).

    Attributes:
        ifuslot: Unique identifier for the IFU slot.
        x: Float representation of the X-coordinate position.
        y: Float representation of the Y-coordinate position.
        specid: Integer identifier for the spectrograph.
        specslot: Integer of the specific slot on the spectrograph.
        ifuid: String identifier for the IFU.
        ifurot: Float value denoting the rotation value of the IFU.
        platescl: Float value for the plate scale.
        xid: Integer extracted from the first two characters of 'ifuslot'.
        yid: Integer derived from the third character of 'ifuslot'.
    """
    def __init__(self, ifuslot, x, y, specid, specslot,
                 ifuid, ifurot, platescl):
        self.ifuslot = ifuslot
        self.x = float(x)
        self.y = float(y)
        self.specid = int(specid)
        self.specslot = int(specslot)
        self.ifuid = str(ifuid)
        self.ifurot = float(ifurot)
        self.platescl = float(platescl)
        self.xid = int(self.ifuslot[0:2])
        self.yid = int(self.ifuslot[2])

class FPlane(object):
    """
    Represents a focal plane configuration containing multiple Integral Field Units (IFUs).

    This class is responsible for loading focal plane configuration data from a file, creating
    and managing IFU objects, and providing access to IFUs by their slot identifiers. It allows
    for customization, excluding specific IFUs or skipping empty IFUs during loading.

    Attributes:
        _IFU: The class used to create IFU instances.
        _ifus_by_slot: A dictionary mapping IFU slots to their respective IFU objects.
    """
    def __init__(self, fplane_file, ifu_class=IFU, empty_specid='00',
                 empty_ifuid='000', exclude_ifuslot=[], skip_empty=False):
        self._fplane_file = fplane_file
        self._IFU = ifu_class
        self._ifus_by_slot = {}
        self._load_fplane(fplane_file, empty_specid, empty_ifuid,
                          exclude_ifuslot, skip_empty)

    def by_ifuslot(self, ifuslot):
        """Return the IFU object for a given slot identifier.

        Args:
            ifuslot: IFU slot identifier string.

        Returns:
            IFU: The IFU instance corresponding to ``ifuslot``.
        """
        return self._ifus_by_slot[ifuslot]

    def _load_fplane(self, fname, empty_specid, empty_ifuid, exclude_ifuslot,
                     skip_empty):
        """Load focal-plane (fplane) configuration and populate IFU map.

        Reads whitespace-separated columns while skipping comments, applies
        optional filtering and placeholder substitution for empty spec/IFU ids,
        and registers each IFU via ``add_ifu``.

        Args:
            fname: Path to the fplane configuration file.
            empty_specid: String that marks an empty spectrograph ID.
            empty_ifuid: String that marks an empty IFU ID.
            exclude_ifuslot: Iterable of IFU slot identifiers to skip entirely.
            skip_empty: If True, rows with empty specid or ifuid are skipped
                rather than filled with placeholders.

        Raises:
            ValueError: If the file cannot be parsed.
        """
        missing = 1
        data = np.genfromtxt(fname, dtype=str, comments="#", autostrip=True)

        for row in data:
            params = [str(i) for i in row[:8]]  # enforce 8 fields as strings
            if not params:
                continue
            if params[0] in exclude_ifuslot:
                continue

            changed = False
            if len(params) > 3 and params[3] == empty_specid:
                params[3] = '-%02d' % missing
                changed = True
                if skip_empty:
                    continue
            if len(params) > 5 and params[5] == empty_ifuid:
                params[5] = 'N%02d' % missing
                changed = True
                if skip_empty:
                    continue
            if changed:
                missing += 1
            self.add_ifu(params)

    def add_ifu(self, ifu_params):
        """
        Adds an IFU (Integral Field Unit) object to the appropriate mappings.

        The IFU data, provided as a tuple of parameters, is used to initialize
        an IFU object, which is then stored in multiple mapping structures
        for efficient access based on its identifier, slot, or spectrograph ID.

        Args:
            ifu_params (tuple): A tuple containing the parameters required to create
                an IFU object. The tuple must include:
                ifuslot: The identifier of the IFU slot.
                x: The x-coordinate of the IFU.
                y: The y-coordinate of the IFU.
                specid: The spectrograph ID to which the IFU belongs.
                speclot: The spectrograph lot identifier.
                ifuid: The unique identifier of the IFU.
                ifurot: The rotation value of the IFU.
                platescl: The plate scale of the IFU.
        """
        ifuslot, x, y, specid, speclot, ifuid, ifurot, platescl = ifu_params
        _ifu = self._IFU(ifuslot, x, y, specid, speclot,
                         ifuid, ifurot, platescl)
        self._ifus_by_slot[_ifu.ifuslot] = _ifu
