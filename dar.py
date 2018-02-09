""" Panacea - dar.py


1) Measures differential atmospheric refraction

.. moduleauthor:: Greg Zeimann <gregz@astro.as.utexas.edu>

"""

import numpy as np
from astropy.modeling.models import Polynomial2D, Moffat2D
from astropy.modeling.fitting import LevMarLSQFitter
from utils import is_outlier


class Dar:
    ''' Differential atmospheric refraction from bright star '''
    def __init__(self, x, y, spec, wave, wavebinsize=100,
                 polyorder=3, rectified_dlam=None):
        '''
        Parameters:
        x : array
            fiber x position in IFU
        y : array
            fiber y position in IFU
        spec : 2-d array
            spectrum for each fiber
        wave : 2-d array
            wavelength for each fiber
        dar_wavebinsize : float
            Size of the wavelength bin to collapse spectra for measuring
            spatial PSF model including the x,y position in the IFU.
        '''
        self.x = x
        self.y = y
        self.spec = spec
        self.wave = wave
        self.model = Polynomial2D(degree=1) + Moffat2D()
        self.fitter = LevMarLSQFitter()
        self.wavebinsize = wavebinsize
        self.polyorder = polyorder
        self.rectified_dlam = rectified_dlam
        self.rectify()

    def rectify(self):
        ''' Rectify spectra to same "rect_wave" '''
        dlam = np.zeros(self.wave.shape)
        dlam[:, 1:] = np.diff(self.wave, axis=1)
        dlam[:, 0] = dlam[:, 1]
        if self.rectified_dlam is None:
            self.rectified_dlam = np.nanmedian(dlam)
        self.rect_wave = np.arange(int(self.wave.min()), int(self.wave.max()),
                                   self.rectified_dlam)
        self.rect_spec = np.zeros((self.spec.shape[0], len(self.rect_wave)))
        for i in np.arange(self.spec.shape[0]):
            self.rect_spec[i, :] = np.interp(self.rect_wave, self.wave[i],
                                             self.spec[i] / dlam[i],
                                             left=np.nan, right=np.nan)

    def iter_poly_fit(self, x, y, niter=3):
        sel = np.ones(x.shape, dtype=bool)
        sel[0], sel[-1] = (False, False)
        for i in np.arange(niter):
            p0 = np.polyfit(x[sel], y[sel], self.polyorder)
            ym = np.polyval(p0, x)
            sel = ~is_outlier(ym - y)
        p0 = np.polyfit(x[sel], y[sel], self.polyorder)
        return np.polyval(p0, x)

    def measure_dar(self):
        ''' Measure PSF model as a function of wavelength '''
        bins = np.arange(self.wave.min(), self.wave.max()+self.wavebinsize,
                         self.wavebinsize)
        xc, yc, wc, fwhm, gamma, alpha = ([], [], [], [], [], [])
        for i, v in enumerate(bins[:-1]):
            xl = np.searchsorted(self.rect_wave, bins[i], side='left')
            xh = np.searchsorted(self.rect_wave, bins[i+1], side='right')
            self.Z = np.nansum(self.rect_spec[:, xl:xh], axis=1)
            self.model.x_0_1.value = 0.0
            self.model.y_0_1.value = 0.0
            self.model.c1_0_0.value = 0.0
            self.model.c0_1_0.value = 0.0
            self.model.gamma_1.value = 1.
            self.model.alpha_1.value = 1.
            self.model.c0_0_0 = np.nanmedian(self.Z)
            self.model.amplitude_1 = np.nanpercentile(self.Z, 99)
            fit = self.fitter(self.model, self.x, self.y, self.Z)
            self.model = fit
            xc.append(fit.x_0_1.value)
            yc.append(fit.y_0_1.value)
            wc.append((bins[i] + bins[i+1]) / 2.)
            a = self.model.gamma_1.value
            b = self.model.alpha_1.value
            gamma.append(a)
            alpha.append(b)
            fwhm.append(np.abs(a * 2 * np.sqrt(2**(1./b) - 1.)))
        xc, yc, wc, gamma, alpha, fwhm = [np.array(i) for i in [xc, yc, wc,
                                          gamma, alpha, fwhm]]
        self.dar_x = self.iter_poly_fit(wc, xc)
        self.dar_y = self.iter_poly_fit(wc, yc)
        self.dar_fwhm = self.iter_poly_fit(wc, fwhm)
        self.dar_gamma = self.iter_poly_fit(wc, gamma)
        self.dar_alpha = self.iter_poly_fit(wc, alpha)
        self.dar_wave = wc

    def psfextract(self, boxsize_forintegral=40., boxsize_gridlength=100,
                   xoff=0.0, yoff=0.0):
        ''' Extract spectrum using psf model from DAR calculation '''
        xp, yp = (np.linspace(-boxsize_forintegral/2., boxsize_forintegral/2.,
                              int(boxsize_gridlength+1)),
                  np.linspace(-boxsize_forintegral/2., boxsize_forintegral/2.,
                              int(boxsize_gridlength+1)))

        xgrid, ygrid = np.meshgrid(xp, yp)
        self.flux = np.zeros(self.rect_wave.shape)
        for i, v in enumerate(self.rect_wave):
            self.model.x_0_1.value = np.interp(self.rect_wave[i],
                                               self.dar_wave, self.dar_x+xoff)
            self.model.y_0_1.value = np.interp(self.rect_wave[i],
                                               self.dar_wave, self.dar_y+yoff)
            self.model.c1_0_0.value = 0.0
            self.model.c0_1_0.value = 0.0
            self.model.gamma_1.value = np.interp(self.rect_wave[i],
                                                 self.dar_wave, self.dar_gamma)
            self.model.alpha_1.value = np.interp(self.rect_wave[i],
                                                 self.dar_wave, self.dar_alpha)
            pnames = ['x_0_1', 'y_0_1', 'gamma_1', 'alpha_1']
            for pn in pnames:
                setattr(getattr(self.model, pn), 'fixed', True)
            self.model.c0_0_0 = np.nanmedian(self.rect_spec[:, i])
            self.model.amplitude_1 = np.nanpercentile(self.rect_spec[:, i], 99)
            fit = self.fitter(self.model, self.x, self.y, self.rect_spec[:, i])
            self.model = fit
            self.flux[i] = (self.model(xp.ravel(), yp.ravel()).sum() *
                            (boxsize_forintegral / boxsize_gridlength)**2)
