""" Panacea - throughput.py


1) Measures differential atmospheric refraction

.. moduleauthor:: Greg Zeimann <gregz@astro.as.utexas.edu>

"""

import numpy as np
import os.path as op
from fiber_utils import bspline_x0
from utils import biweight_location, is_outlier
from astropy.io import fits
from dar import Dar
from scipy.interpolate import splrep, splev
from scipy.signal import savgol_filter
from input_utils import setup_logging
from astrometry import Astrometry

# from telluricabs import TelluricAbs
try:
    from pyhetdex.het.telescope import HetpupilModel
    hetpupil_installed = True
except ImportError:
    print('Cannot find HETpupilModel.  Please check pyhetdex installation.')
    print('For now, using default 50m**2 for mirror illumination')
    hetpupil_installed = False


class ReduceLRS2:
    ''' Wrapper for reduction routines with processed data, multi*.fits '''
    def __init__(self, base_filename, side, use_twi=True,
                 standard_folder='/Users/gregz/cure/virus_early/virus_config/'
                 'standards', fplane_file=None):
        '''
        This serves as a wrapper for astrometry, DAR, extraction, and
        throughput.

        Parameters:
        base_filename : str
            The full path of the file, "multi_*_*_*", from the initial
            panacea2.py reduction. For example:
            "/work/03946/hetdex/maverick/reductions/20180206/lrs2/lrs20000035/"
            "exp01/lrs2/multi_503_056_7001"
            The "_{amp}.fits" is left off.
        side : str
            Either 'BL', 'BR', 'RL', 'RR'.  The input is case insensitive.
            This will choose either uv, orange, red or farred depending.
        '''
        self.base_filename = base_filename
        self.side = side
        self.standard_folder = standard_folder
        self.use_twi = use_twi
        self.fplane_file = fplane_file
        self.ifuslot = self.basename(self.base_filename).split('_')[2]
        self.log = setup_logging('throughput')
        if self.side.lower() == 'bl':
            self.amps = ['LL', 'LU']
            self.wave_lims = [3640., 4610.]
        elif self.side.lower() == 'br':
            self.amps = ['RL', 'RU']
            self.wave_lims = [4680., 6960.]
        if self.side.lower() == 'rl':
            self.amps = ['LL', 'LU']
            self.wave_lims = [6450., 8400.]
        elif self.side.lower() == 'rr':
            self.amps = ['RL', 'RU']
            self.wave_lims = [8275., 10500.]
        self.read_in_files()
        self.get_astrometry()

    def get_astrometry(self):
        ''' Set the RA, dec for each Fiber '''
        self.astrom = Astrometry(self.header['TRAJCRA']*15.,
                                 self.header['TRAJCDEC'],
                                 self.header['PARANGLE'], 0., 0.,
                                 fplane_file=self.fplane_file)
        self.ra, self.dec = self.astrom.get_ifuspos_ra_dec(self.ifuslot,
                                                           self.ifux,
                                                           self.ifuy)

    def get_dar_model(self):
        '''
        Use the Dar class to retreive or measure the Dar, rectifyied the
        spectrum, and extract a source using a PSF model.
        '''
        self.dar = Dar(self.ifux, self.ifuy, self.spec, self.wave)
        self.dar.measure_dar()
        self.dar.psfextract()
        list_var = self.restrict_wavelengths(self.dar.rect_wave, self.dar.flux,
                                             self.dar.back, self.dar.rect_spec)
        [self.dar.rect_wave, self.dar.flux,
         self.dar.back, self.dar.rect_spec] = list_var

    def restrict_wavelengths(self, wave, spec, back, rect_spec):
        ''' Cowardly avoid wavelengths with rapid response changes '''
        sel = np.where((wave > self.wave_lims[0]) *
                       (wave < self.wave_lims[1]))[0]
        return wave[sel], spec[sel], back[sel], rect_spec[:, sel]

    def get_telluric_abs(self):
        ''' TESTING MODE STILL '''
        self.telabs = TelluricAbs(self.dar.rect_wave, self.clam, self.RH,
                                  self.T, self.P, self.ZD)
        self.telabs.fit_telluric_abs()

    def read_in_files(self):
        '''
        Read in the multi* fits files for each amp and build the x, y positions
        in the ifu as well as the spectrum (corrected for fiber to fiber)
        and wavelength.
        '''
        self.log.info('Reading in initial reductions from %s' %
                      self.base_filename)
        x, y, spec, wave, twi = ([], [], [], [], [])
        for amp in self.amps:
            fn = self.base_filename + ('_%s.fits' % amp)
            F = fits.open(fn)
            x.append(F['ifupos'].data[:, 0])
            y.append(F['ifupos'].data[:, 1])
            spec.append(F['spectrum'].data)
            twi.append(F['twi_spectrum'].data)
            wave.append(F['wavelength'].data)
        self.header = dict(F[0].header)
        self.object = F[0].header['OBJECT'][:-6]
        self.exptime = F[0].header['EXPTIME']
        self.RH = F[0].header['HUMIDITY']
        self.T = F[0].header['AMBTEMP']
        self.P = F[0].header['BAROMPRE']
        self.ZD = F[0].header['ZD']
        self.ifux = np.hstack(x)
        self.ifuy = np.hstack(y)
        self.spec = np.vstack(spec)
        self.wave = np.vstack(wave)
        twi = np.vstack(twi)
        if self.use_twi:
            spec = twi * 1.
        else:
            spec = self.spec * 1.
        self.ftf = self.get_fiber_to_fiber(self.wave, spec)
        for i in np.arange(self.spec.shape[0]):
            self.spec[i, :] /= self.ftf[i, :]

    def get_standard_spectrum_from_file(self):
        ''' Read standard spectrum for self.object and convert to f_lam '''
        self.log.info('Grabbing standard spectrum for %s' % self.object)
        filename = op.join(self.standard_folder,
                           'm' + self.object.lower() + '.dat.txt')
        wave, standardmag = np.loadtxt(filename, usecols=(0, 1), unpack=True)
        fnu = 10**(0.4 * (-48.6 - standardmag))
        self.standard_flam = fnu * 2.99792e18 / wave**2
        self.standard_wave = wave

    def get_mirror_illumination(self, fn=None):
        ''' Use Hetpupil from Cure to calculate mirror illumination (cm^2) '''
        self.log.info('Getting mirror illumination')
        if hetpupil_installed:
            self.log.info('Using HetpupilModel from pyhetdex')
            if fn is None:
                fn = self.base_filename + ('_%s.fits' % self.amps[0])
            mirror_illum = HetpupilModel([fn], normalize=False)
            self.area = mirror_illum.fill_factor[0] * 55. * 1e4
        else:
            self.log.info('Using default mirror illumination value')
            self.area = 50. * 1e4
        self.log.info('Mirror illumination: %0.2f m^2' % (self.area/1e4))

    def convert_units(self):
        ''' convert cnts/A to cnts/A/s/cm^2 '''
        self.log.info('Converting cnts/A to cnts/A/s/cm^2')
        if not hasattr(self, 'area'):
            self.get_mirror_illumination()
        self.clam = self.dar.flux / self.exptime / self.area

    def compare_spectrum_to_standard(self):
        ''' Bin measured clam spectrum and calculate response, R '''
        self.log.info('Getting response function')
        if not hasattr(self, 'area'):
            self.get_mirror_illumination()
        if not hasattr(self, 'clam'):
            self.convert_units()
        if not hasattr(self, 'standard_wave'):
            self.get_standard_spectrum_from_file()
        xl = np.searchsorted(self.standard_wave, self.dar.rect_wave.min(),
                             side='left')
        xh = np.searchsorted(self.standard_wave, self.dar.rect_wave.max(),
                             side='right')
        self.binned_clam = np.interp(self.standard_wave[xl:xh],
                                     self.dar.rect_wave, self.clam)
        self.R = self.standard_flam[xl:xh] / self.binned_clam
        self.R_wave = self.standard_wave[xl:xh]
        B, c = bspline_x0(self.R_wave, nknots=11)
        sol = np.linalg.lstsq(c, self.R)[0]
        self.smooth_R = np.dot(c, sol)

    def adjust_ftf_from_twi_to_sky(self):
        '''
        Use pre-calculated relation which converts fiber to fiber calculated
        from the twilight frame to that of the a night sky.  The differences
        are not yet understood
        '''
        self.log.info('Adjusting Fiber to Fiber')
        pass

    def get_fiber_to_fiber(self, wavelength, spec, fac=10, knots=15):
        '''
        Get fiber to fiber from twilight or sky (input spec)

        Parameters
        ----------
        wavelength : 2-d array
            Wavelength for each fiber
        spec : 2-d array
            Spectrum (twi or sky) for each fiber
        fac : int
            Resample to a finer grid by this factor
        knots : int
            Number of knots in B-spline for smoothing the fiber to fiber
            evaluations
        '''
        self.log.info('Getting Fiber to Fiber')
        mn, mx = (wavelength.min(), wavelength.max())
        N, D = spec.shape
        wv = np.linspace(mn, mx, fac*D)
        xs = np.linspace(0, 1, D*fac)
        A = np.zeros((len(xs), N))
        i = 0
        for sp, wave in zip(spec, wavelength):
            y = sp
            xp = np.interp(wave, np.linspace(mn, mx, D*fac), xs)
            tck = splrep(xp, y)
            A[:, i] = splev(xs, tck)
            i += 1
        ys = biweight_location(A, axis=(1,))
        Anew = A / ys[:, np.newaxis]
        B, c = bspline_x0(wv, nknots=knots)
        smooth = np.zeros(Anew.shape)
        for i in np.arange(Anew.shape[1]):
            xl = np.searchsorted(wv, wavelength[i].min(), side='right')
            xh = np.searchsorted(wv, wavelength[i].max(), side='left')
            m = savgol_filter(Anew[xl:xh, i], 155, 1)
            if np.median(np.abs(m - Anew[xl:xh, i])) == 0.:
                sel = np.arange(0, xh-xl, dtype=int)
            else:
                sel = np.where(is_outlier(m - Anew[xl:xh, i]) == 0)[0]
            sol = np.linalg.lstsq(c[xl:xh, :][sel, :], Anew[xl:xh, i][sel],
                                  rcond=None)[0]
            smooth[:, i] = np.dot(c, sol)
        ftf = np.zeros(spec.shape)
        for i in np.arange(smooth.shape[1]):
            ftf[i, :] = np.interp(wavelength[i, :], wv, smooth[:, i])
        return ftf
