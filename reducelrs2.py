""" Panacea - reducelrs2.py


1) STEPS

.. moduleauthor:: Greg Zeimann <gregz@astro.as.utexas.edu>

"""

import numpy as np
import os.path as op
from fiber_utils import bspline_x0
from utils import biweight_location, is_outlier
from astropy.io import fits
from dar import Dar
from scipy.interpolate import splrep, splev, interp1d
from scipy.signal import savgol_filter
from input_utils import setup_logging
from astrometry import Astrometry
import astropy.coordinates as coord
import astropy.units as u

try:
    from telluricabs import TelluricAbs
except ImportError:
    print('Telluric Fitter not accessible, please check installation')

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
                 config_folder='/Users/gregz/cure/virus_early/virus_config/',
                 fplane_file=None, ftf_file=None):
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
        self.path = op.dirname(base_filename)
        self.multi_name = op.basename(base_filename)
        self.side = side
        self.config_folder = config_folder
        self.standard_folder = op.join(config_folder, 'standards')
        self.ftfcor_folder = op.join(config_folder, 'FtFcor')
        self.use_twi = use_twi
        self.fplane_file = fplane_file
        self.ftf_file = ftf_file
        self.goodfibers = None
        self.ifuslot = op.basename(self.base_filename).split('_')[2]
        self.log = setup_logging('throughput')
        if self.side.lower() == 'bl':
            self.instr_name = 'lrs2_uv'
            self.side_name = 'uv'
            self.amps = ['LL', 'LU']
            self.wave_lims = [3640., 4640.]
        elif self.side.lower() == 'br':
            self.instr_name = 'lrs2_orange'
            self.side_name = 'orange'
            self.amps = ['RU', 'RL']
            self.wave_lims = [4655., 6960.]
        if self.side.lower() == 'rl':
            self.instr_name = 'lrs2_red'
            self.side_name = 'red'
            self.amps = ['LL', 'LU']
            self.wave_lims = [6450., 8400.]
        elif self.side.lower() == 'rr':
            self.instr_name = 'lrs2_farred'
            self.side_name = 'farred'
            self.amps = ['RU', 'RL']
            self.wave_lims = [8275., 10500.]
        self.read_in_files()
        self.get_astrometry()
        self.dar = Dar(self.ifux, self.ifuy, self.spec, self.wave,
                       goodfibers=self.goodfibers)

    def get_extinction_mag(self, ra, dec, wave):
        '''
        Get galactic extinction in magnitudes for ra, dec

        Parameters
        ----------
        ra : float
            Right Ascension in degrees
        dec : float
            Declination in degrees
        wave : array
            Wavelengths for galactic extinction

        Outputs
        -------
        Amag : array
            Galactic Extinction in magnitudes for input ra, dec at
            wavelengths set by input wave
        '''        
        from astroquery.irsa_dust import IrsaDust
        coo = coord.SkyCoord(ra*u.deg, dec*u.deg, frame='fk5')
        table = IrsaDust.get_extinction_table(coo)
        Amag = np.interp(wave, table['LamEff'], table['A_SFD'])
        return Amag

    def get_astrometry(self):
        ''' Set the RA, dec for each Fiber '''
        try:
            self.astrom = Astrometry(self.header_dict['TRAJCRA']*15.,
                                     self.header_dict['TRAJCDEC'],
                                     self.header_dict['PARANGLE'], 0., 0.,
                                     fplane_file=self.fplane_file)
            self.ra, self.dec = self.astrom.get_ifuspos_ra_dec(self.ifuslot,
                                                               self.ifuy,
                                                               self.ifux)
        except:
            self.log.warning('Astrometry did not work because of -999999.')
            self.ra, self.dec = (self.ifux, self.ifuy)

    def get_dar_model(self):
        '''
        Use the Dar class to retreive or measure the Dar, rectifyied the
        spectrum, and extract a source using a PSF model.
        '''
        self.log.info('Measuring DAR model and rectifying spectra')
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
        x, y, spec, wave, twi, tr, im, er = ([], [], [], [], [], [], [], [])
        for amp in self.amps:
            fn = self.base_filename + ('_%s.fits' % amp)
            F = fits.open(fn)
            x.append(F['ifupos'].data[:, 0])
            y.append(F['ifupos'].data[:, 1])
            im.append(F['PRIMARY'].data*1.)
            er.append(F['error'].data*1.)
            spec.append(F['spectrum'].data)
            twi.append(F['twi_spectrum'].data)
            wave.append(F['wavelength'].data)
            if amp in ['LU', 'RL']:
                addtr = F[0].data.shape[0]
            else:
                addtr = 0.
            tr.append(F['trace'].data + addtr)
        self.header_dict = dict(F[0].header)
        self.header = F[0].header
        self.object = F[0].header['OBJECT'][:-6]
        self.exptime = F[0].header['EXPTIME']
        self.RH = F[0].header['HUMIDITY']
        self.T = F[0].header['AMBTEMP']
        self.P = F[0].header['BAROMPRE']
        self.ZD = F[0].header['ZD']
        self.ifux = np.array(np.hstack(x), dtype='float64')
        self.ifuy = np.array(np.hstack(y), dtype='float64')
        self.spec = np.array(np.vstack(spec), dtype='float64')
        self.oldspec = self.spec * 1.
        self.wave = np.array(np.vstack(wave), dtype='float64')
        self.trace = np.array(np.vstack(tr), dtype='float64')
        self.image_name = np.array(np.vstack(im), dtype='float64')
        self.error = np.array(np.vstack(er), dtype='float64')
        self.twi = np.array(np.vstack(twi), dtype='float64')
        self.define_good_fibers()

    def get_standard_spectrum_from_file(self):
        ''' Read standard spectrum for self.object and convert to f_lam '''
        self.log.info('Grabbing standard spectrum for %s' % self.object)
        filename = op.join(self.standard_folder,
                           'm' + self.object.lower() + '.dat.txt')
        wave, standardmag = np.loadtxt(filename, usecols=(0, 1), unpack=True)
        fnu = 10**(0.4 * (-48.6 - standardmag))
        self.standard_flam = fnu * 2.99792e18 / wave**2
        self.standard_wave = wave
        test1 = np.any(self.standard_wave > 7500)
        test2 = np.all(self.standard_wave < 10500)
        if (test1 * test2):
            sel = np.where(self.standard_wave > 7500)[0]
            p0 = np.polyfit(self.standard_wave[sel],
                            np.log10(self.standard_flam[sel]), 1)
            mdiff = np.median(np.diff(self.standard_wave))
            x = np.arange(self.standard_wave[-1] + mdiff, 10500 + mdiff, mdiff)
            yext = 10**np.polyval(p0, x)
            self.standard_wave = np.hstack([self.standard_wave, x])
            self.standard_flam = np.hstack([self.standard_flam, yext])

    def get_mirror_illumination(self, fn=None):
        ''' Use Hetpupil from Cure to calculate mirror illumination (cm^2) '''
        self.log.info('Getting mirror illumination')
        if hetpupil_installed:
            self.log.info('Using HetpupilModel from pyhetdex')
            if fn is None:
                fn = self.base_filename + ('_%s.fits' % self.amps[0])
            try:
                mirror_illum = HetpupilModel([fn], normalize=False)
                self.area = mirror_illum.fill_factor[0] * 55. * 1e4
            except:
                self.log.info('Using default mirror illumination value')
                self.area = 50. * 1e4
        else:
            self.log.info('Using default mirror illumination value')
            self.area = 50. * 1e4
        self.log.info('Mirror illumination: %0.2f m^2' % (self.area/1e4))

    def get_error(self, fac=0.7):
        self.dar.rect_spec_error = np.sqrt((3.*1.15)**2 +
                                           self.dar.rect_spec*fac)

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
                             side='left') - 1
        xh = np.searchsorted(self.standard_wave, self.dar.rect_wave.max(),
                             side='right') + 1
        indices = np.digitize(self.standard_wave[xl:xh], self.dar.rect_wave)
        bins = np.array_split(self.clam, indices)
        self.binned_clam = np.hstack([np.median(bin0) for bin0 in bins])[1:]
        self.R = self.standard_flam[xl:xh] / self.binned_clam
        interpolator = interp1d(self.standard_wave[xl:xh], self.R,
                                kind='cubic', fill_value='extrapolate')
        self.R = interpolator(self.dar.rect_wave)
        self.R = (np.interp(self.dar.rect_wave, self.standard_wave,
                            self.standard_flam) / self.clam)
        self.R_wave = self.dar.rect_wave
        size = int(250. / np.diff(self.R_wave).mean())
        if size % 2 == 0:
            size = size + 1
        if size <= 1:
            size += 2
        self.smooth_R = savgol_filter(self.R, size, 1)

    def adjust_ftf_from_twi_to_sky(self):
        '''
        Use pre-calculated relation which converts fiber to fiber calculated
        from the twilight frame to that of the a night sky.  The differences
        are not yet understood
        '''
        self.log.info('Adjusting Fiber to Fiber')
        Cor = np.loadtxt(op.join(self.ftfcor_folder, self.instr_name + '.txt'))
        if hasattr(self, 'ftf') and self.use_twi:
            self.ftf = self.ftf + Cor[:, 1:2]

    def define_good_fibers(self, thresh=0.5):
        if hasattr(self, 'ftf'):
            y = np.nanmedian(self.ftf, axis=1)
            self.goodfibers = np.where(y > thresh)[0]
        else:
            self.goodfibers = np.arange(self.spec.shape[0], dtype=int)

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
        self.inspect_wv = wv
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
        self.inspect = Anew
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
            sol = np.linalg.lstsq(c[xl:xh, :][sel, :], Anew[xl:xh, i][sel])[0]
            smooth[:, i] = np.dot(c, sol)
        ftf = np.zeros(spec.shape)
        for i in np.arange(smooth.shape[1]):
            ftf[i, :] = np.interp(wavelength[i, :], wv, smooth[:, i])
        return ftf

    def write_to_fits(self, hdu, outname):
        '''
        Writing fits file to outname
        '''
        try:
            hdu.writeto(outname, overwrite=True)
        except TypeError:
            hdu.writeto(outname, clobber=True)

    def write_header(self, hdu):
        for key in self.header_dict.keys():
            if key in hdu.header:
                continue
            if ('CCDSEC' in key) or ('DATASEC' in key):
                continue
            hdu.header[key] = self.header_dict[key]
        return hdu

    def save(self, image_list=[], name_list=[]):
        ''' Save the multi*.fits file '''
        fn = op.join(self.path, '%s_%s.fits' % (self.multi_name,
                                                self.side_name))
        fits_list = []
        if len(image_list) != len(name_list):
            return None
        for i, image in enumerate(image_list):
            if i == 0:
                fits_list.append(fits.PrimaryHDU(np.array(getattr(self, image),
                                                          dtype='float32')))
            else:
                fits_list.append(fits.ImageHDU(np.array(getattr(self, image),
                                                        dtype='float32')))

            fits_list[-1].header['EXTNAME'] = name_list[i]
        if fits_list:
            fits_list[0] = self.write_header(fits_list[0])
            hdu = fits.HDUList(fits_list)
            self.write_to_fits(hdu, fn)
