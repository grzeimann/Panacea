""" Panacea - throughput.py


1) Measures differential atmospheric refraction

.. moduleauthor:: Greg Zeimann <gregz@astro.as.utexas.edu>

"""

import numpy as np
import os.path as op
from utils import biweight_bin
from fiber_utils import bspline_x0
from astropy.io import fits
from dar import Dar
from telluricabs import TelluricAbs
try:
    from pyhetdex.het.telescope import HetpupilModel
    hetpupil_installed = True
except ImportError:
    print('Cannot find HETpupilModel.  Please check pyhetdex installation.')
    print('For now, using default 50m**2 for mirror illumination')
    hetpupil_installed = False


class Throughput:
    ''' Throughput for standard stars with LRS2 at the HET '''
    def __init__(self, base_filename, side,
                 standard_folder='/Users/gregz/cure/virus_early/virus_config/'
                 'standards'):
        '''
        Aimed at calculating the relative throughput for a given standard star
        which can be applied to science frames.

        Parameters:
        base_filename : str
            The full path of the file, "multi_*_*_*", from the initial
            panacea2.py reduction. For example:
            "/work/03946/hetdex/maverick/reductions/20180206/lrs2/lrs20000035/"
            "exp01/lrs2/multi_503_056_7001"
            The "_{amp}.fits" is left off.
        side : str
            Either 'L' or 'R'.  The input is case insensitive.  This will
            choose either the uv vs. orange or red vs. farred depending
            on the base_filename which selects between LRS2-B and LRS2-R.
        '''
        self.base_filename = base_filename
        self.side = side
        self.standard_folder = standard_folder
        if self.side.lower() == 'l':
            self.amps = ['LL', 'LU']
            self.wave_lims = [3640., 4610.]
        elif self.side.lower() == 'r':
            self.amps = ['RL', 'RU']
            self.wave_lims = [4650., 7000.]
        self.read_in_files()

    def get_dar_model(self):
        self.dar = Dar(self.ifux, self.ifuy, self.spec, self.wave)
        self.dar.measure_dar()
        self.dar.psfextract()
        self.dar.rect_wave, self.dar.flux = self.restrict_wavelengths(
                                                            self.dar.rect_wave,
                                                            self.dar.flux)

    def restrict_wavelengths(self, wave, spec):
        sel = np.where((wave > self.wave_lims[0]) *
                       (wave < self.wave_lims[1]))[0]
        return wave[sel], spec[sel]

    def get_telluric_abs(self):
        self.telabs = TelluricAbs(self.dar.rect_wave, self.clam, self.RH,
                                  self.T, self.P, self.ZD)
        self.telabs.fit_telluric_abs()

    def read_in_files(self):
        '''
        Read in the multi* fits files for each amp and build the x, y positions
        in the ifu as well as the spectrum (corrected for fiber to fiber)
        and wavelength.
        '''
        x, y, spec, wave = ([], [], [], [])
        for amp in self.amps:
            fn = self.base_filename + ('_%s.fits' % amp)
            F = fits.open(fn)
            x.append(F['ifupos'].data[:, 0])
            y.append(F['ifupos'].data[:, 1])
            spec.append(F['spectrum'].data / F['fiber_to_fiber'].data)
            wave.append(F['wavelength'].data)
        self.object = F[0].header['OBJECT'].split('_')[0]
        self.exptime = F[0].header['EXPTIME']
        self.RH = F[0].header['HUMIDITY']
        self.T = F[0].header['AMBTEMP']
        self.P = F[0].header['BAROMPRE']
        self.ZD = F[0].header['ZD']
        self.ifux = np.hstack(x)
        self.ifuy = np.hstack(y)
        self.spec = np.vstack(spec)
        self.wave = np.vstack(wave)

    def get_standard_spectrum_from_file(self):
        ''' Read standard spectrum for self.object and convert to f_lam '''
        filename = op.join(self.standard_folder,
                           'm' + self.object.lower() + '.dat.txt')
        wave, standardmag = np.loadtxt(filename, usecols=(0, 1), unpack=True)
        fnu = 10**(0.4 * (-48.6 - standardmag))
        self.standard_flam = fnu * 2.99792e18 / wave**2
        self.standard_wave = wave

    def get_mirror_illumination(self, fn=None):
        ''' Use Hetpupil from Cure to calculate mirror illumination (cm^2) '''
        if hetpupil_installed:
            if fn is None:
                fn = self.base_filename + ('_%s.fits' % self.amps[0])
            mirror_illum = HetpupilModel([fn], normalize=False)
            self.area = mirror_illum.fill_factor[0] * 55. * 1e4
        else:
            self.area = 50. * 1e4

    def convert_units(self):
        ''' convert cnts/A to cnts/A/s/cm^2 '''
        if not hasattr(self, 'area'):
            self.get_mirror_illumination()
        self.clam = self.dar.flux / self.exptime / self.area

    def compare_spectrum_to_standard(self):
        ''' Bin measured clam spectrum and calculate response, R '''
        if not hasattr(self, 'area'):
            self.get_mirror_illumination()
        if not hasattr(self, 'clam'):
            self.convert_units()
        xl = np.searchsorted(self.standard_wave, self.dar.rect_wave.min(),
                             side='left')
        xh = np.searchsorted(self.standard_wave, self.dar.rect_wave.max(),
                             side='right')
        self.binned_clam = np.array(biweight_bin(self.standard_wave[xl:xh],
                                                 self.dar.rect_wave,
                                                 self.clam))
        self.R = self.standard_flam[xl:xh] / self.binned_clam
        self.R_wave = self.standard_wave[xl:xh]
        sel = np.where(np.isfinite(self.R) * (self.binned_clam >
                       (.4*np.nanmedian(self.binned_clam))))[0]
        self.R_wave = self.R_wave[sel]
        self.R = self.R[sel]
        B, c = bspline_x0(self.R_wave, nknots=25)
        sol = np.linalg.lstsq(c, self.R)[0]
        self.smooth_R = np.dot(c, sol)
