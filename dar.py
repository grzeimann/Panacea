""" Panacea - dar.py


1) Measures differential atmospheric refraction

.. moduleauthor:: Greg Zeimann <gregz@astro.as.utexas.edu>

"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.modeling.models import Polynomial2D, Moffat2D, Gaussian2D
from asymmoffat import AsymMoffat2D
from astropy.modeling.fitting import LevMarLSQFitter
from astropy.visualization import AsinhStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from scipy.interpolate import splev, splrep
from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import RANSACRegressor
from sklearn.pipeline import Pipeline
from input_utils import setup_logging


class Dar:
    ''' Differential atmospheric refraction from bright star '''
    def __init__(self, x, y, spec, wave, wavebinsize=50,
                 polyorder=3, rectified_dlam=None, psfmodel=AsymMoffat2D,
                 backmodel=Polynomial2D, goodfibers=None):
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
        self.model = backmodel(0) + psfmodel()
        if goodfibers is None:
            self.goodfibers = np.arange(self.spec.shape[0], dtype=int)
        else:
            self.goodfibers = goodfibers
        self.Back, self.PSF = self.model
        self.tinker_params = list(self.PSF.param_names)
        if 'theta' in self.tinker_params:
            self.PSF.theta.bounds = (0, 180.)
        del self.tinker_params[0]
        self.PSF.parameters[3] = 1.5
        self.PSF.parameters[4] = 4.0
        self.default_model_parms = self.model.parameters * 1.

        self.fitter = LevMarLSQFitter()
        self.wavebinsize = wavebinsize
        self.polyorder = polyorder
        self.rectified_dlam = rectified_dlam
        self.rectify()
        self.log = setup_logging('dar')

    def rectify(self, minwave=None, maxwave=None):
        ''' Rectify spectra to same "rect_wave" '''
        dlam = np.zeros(self.wave.shape)
        dlam[:, 1:] = np.diff(self.wave, axis=1)
        dlam[:, 0] = dlam[:, 1]
        if self.rectified_dlam is None:
            self.rectified_dlam = np.nanmedian(dlam)

        self.rect_wave = np.arange(self.wave.min(),
                                   self.wave.max() + self.rectified_dlam,
                                   self.rectified_dlam)
        if minwave is not None and maxwave is not None:
            wnew = np.arange(minwave, maxwave + self.rectified_dlam,
                             self.rectified_dlam)
        else:
            wnew = self.rect_wave * 1.
        self.rect_spec = np.zeros((self.spec.shape[0], len(wnew)))
        xs = np.linspace(0, 1, len(self.rect_wave))
        xn = np.interp(wnew, self.rect_wave, xs)
        for i in np.arange(self.spec.shape[0]):
            y = self.spec[i] / dlam[i]
            xp = np.interp(self.wave[i], self.rect_wave, xs)
            tck = splrep(xp, y)
            self.rect_spec[i, :] = splev(xn, tck)
        self.rect_wave = wnew * 1.

    def grab_archive_psf_model(self):
        pass

    def locate_point_source(self, wavebinsize=None, fixed_list=None,
                            bins=None):
        '''
        Use a polynomial background + moffat profile to model collapsed
        spectra.  Useful for DAR and point source extraction.

        Parameters:
        wavebinsize : float/int
            Size of wavelength bin for measuring background+pointsource model
            Over this bin size, spectra will be collapsed to a single value.
        '''
        if wavebinsize is None:
            wavebinsize = self.wavebinsize
        if not hasattr(self, 'rect_wave'):
            self.rectify()
        if bins is None:
            bins = np.arange(self.wave.min(), self.wave.max(), wavebinsize)
        A = np.zeros((len(bins) - 1, len(self.tinker_params) + 2))
        for i, v in enumerate(bins[:-1]):
            self.fixed_params(parameters=None, value=False)
            if fixed_list is not None:
                self.fixed_params(parameters=fixed_list, value=True)
            xl = np.searchsorted(self.rect_wave, bins[i], side='left')
            xh = np.searchsorted(self.rect_wave, bins[i+1], side='right')
            A[i, 0] = (bins[i] + bins[i+1]) / 2.
            self.Z = np.nanmedian(self.rect_spec[:, xl:xh], axis=1)
            # set model to default
            self.model.parameters = self.default_model_parms * 1.
            self.Back.parameters[0] = np.nanmedian(self.Z[self.goodfibers])
            self.PSF.parameters[0] = np.nanpercentile(self.Z[self.goodfibers],
                                                      99)
            inds = np.argsort(self.Z[self.goodfibers])[-8:-2]
            if not getattr(self.PSF, self.tinker_params[0]).fixed:
                self.PSF.parameters[1] = np.mean(self.x[self.goodfibers][inds])
            else:
                setattr(getattr(self.PSF, self.tinker_params[0]), 'value',
                        np.interp(A[i, 0], self.dar_wave,
                                  getattr(self, 'dar_'+self.tinker_params[0])))
            if not getattr(self.PSF, self.tinker_params[1]).fixed:
                self.PSF.parameters[2] = np.mean(self.y[self.goodfibers][inds])
            else:
                setattr(getattr(self.PSF, self.tinker_params[1]), 'value',
                        np.interp(A[i, 0], self.dar_wave,
                                  getattr(self, 'dar_'+self.tinker_params[1])))
            for j in np.arange(2, len(self.tinker_params)):
                if hasattr(self, 'dar_' + self.tinker_params[j]):
                    y = np.interp(A[i, 0], self.dar_wave,
                                  getattr(self, 'dar_' +
                                  self.tinker_params[j]))
                    setattr(getattr(self.PSF, self.tinker_params[j]), 'value',
                            y)
            fit = self.fitter(self.model, self.x[self.goodfibers],
                              self.y[self.goodfibers], self.Z[self.goodfibers])
            self.model = fit
            self.Back, self.PSF = self.model
            for j, parname in enumerate(self.tinker_params):
                A[i, j+1] = getattr(self.PSF, parname).value

            if hasattr(self.PSF, 'fwhm'):
                A[i, -1] = self.PSF.fwhm
            else:
                A[i, -1] = (np.sqrt(self.PSF.parameters[3]**2 +
                                    self.PSF.parameters[4]**2)*2.355)
        return A

    def robust_poly_fit(self, x, y, order=3):
        model = Pipeline([('poly', PolynomialFeatures(degree=order)),
                          ('linear', RANSACRegressor())])
        try:
            model = model.fit(x[:, np.newaxis], y)
            return model.predict(x[:, np.newaxis])
        except:
            self.log.warn('Polynomial fit failed.  Using median value')
            return np.ones(x.shape) * np.nanmedian(y)

    def measure_dar(self, fixed_list=None):
        ''' Measure PSF model as a function of wavelength '''
        self.A = self.locate_point_source(fixed_list=fixed_list)
        pnames = self.tinker_params + ['fwhm']
        self.dar_wave = self.A[:, 0]

        for j, parname in enumerate(pnames):
            setattr(self, 'dar_'+parname,
                    self.robust_poly_fit(self.dar_wave, self.A[:, j+1]))

    def quick_check_psf_fit(self, x, y, z, psf, back):
            ''' Check fit at a given wavelength, index '''
            fig, (top, middle, bottom) = plt.subplots(3, 1, sharex=True,
                                                      figsize=(6, 15))
            im = top.scatter(x, y, c=z, marker='h',
                             s=220,
                             norm=ImageNormalize(stretch=AsinhStretch()))
            top.axis('equal')
            fig.colorbar(im, ax=top)
            im2 = middle.scatter(x, y, c=(psf + back),
                                 marker='h', s=220,
                                 norm=ImageNormalize(stretch=AsinhStretch()))
            middle.axis('equal')
            fig.colorbar(im2, ax=middle)
            im3 = bottom.scatter(x, y, c=(z - psf - back), marker='h', s=220,
                                 norm=ImageNormalize(stretch=AsinhStretch()))
            bottom.axis('equal')
            fig.colorbar(im3, ax=bottom)
            titles = [r'Data',
                      'Model: Poly Background + Moffat Source',
                      'Residual: Data - Model']
            axs = [top, middle, bottom]
            for ax, title in zip(axs, titles):
                ax.set_title(title)
            plt.show()

    def check_psf_fit(self, outname='check_dar.png', index=500,
                      index_range=None):
        ''' Check fit at a given wavelength, index '''
        fig, (top, middle, bottom) = plt.subplots(3, 1, sharex=True,
                                                  figsize=(6, 15))
        if index_range is not None:
            Z = np.median(self.rect_spec[:, (index-index_range):
                                            (index+index_range)], axis=1)
            M = np.median(self.spec_model[:, (index-index_range):
                                             (index+index_range)], axis=1)
            z = np.median(self.rect_spec[:, (index-index_range):
                                            (index+index_range)] -
                          self.spec_model[:, (index-index_range):
                                             (index+index_range)], axis=1)
        else:
            Z = self.rect_spec[:, index]
            M = self.spec_model[:, index]
            z = Z - M

        min1 = np.percentile(M, 05)
        max1 = np.percentile(M, 95)
        ran = max1 - min1
        vmin1 = min1 - 0.2 * ran
        vmax1 = max1 + 0.2 * ran
        im = top.scatter(self.x, self.y, c=Z, marker='h', vmin=vmin1, vmax=vmax1,
                         s=220, norm=ImageNormalize(stretch=AsinhStretch()))
        top.axis('equal')
        fig.colorbar(im, ax=top)
        im2 = middle.scatter(self.x, self.y, c=M, vmin=vmin1, vmax=vmax1,
                             marker='h', s=220,
                             norm=ImageNormalize(stretch=AsinhStretch()))
        middle.axis('equal')
        fig.colorbar(im2, ax=middle)
        min1 = np.percentile(z, 05)
        max1 = np.percentile(z, 95)
        ran = max1 - min1
        vmin1 = 0 - 0.5 * ran
        vmax1 = 0 + 0.5 * ran
        im3 = bottom.scatter(self.x, self.y, c=z, marker='h', s=220,
                             vmin=vmin1, vmax=vmax1,
                             norm=ImageNormalize(stretch=AsinhStretch()))
        bottom.axis('equal')
        titles = [r'Data at $\lambda = $%0.1f' % self.rect_wave[index],
                  'Model: Poly Background + Moffat Source',
                  'Residual: Data - Model']
        axs = [top, middle, bottom]
        for ax, title in zip(axs, titles):
            ax.set_title(title)
        fig.colorbar(im3, ax=bottom)
        fig.savefig(outname)
        plt.close(fig)

    def fixed_params(self, parameters=None, value=True):
        if parameters is None:
            parameters = self.tinker_params
        for pn in parameters:
            setattr(getattr(self.PSF, pn), 'fixed', value)

    def psfextract(self, boxsize_forintegral=20., boxsize_gridlength=200,
                   xoff=0.0, yoff=0.0):
        ''' Extract spectrum using psf model from DAR calculation '''
        xp, yp = (np.linspace(-boxsize_forintegral/2., boxsize_forintegral/2.,
                              int(boxsize_gridlength+1)),
                  np.linspace(-boxsize_forintegral/2., boxsize_forintegral/2.,
                              int(boxsize_gridlength+1)))

        xgrid, ygrid = np.meshgrid(xp, yp)
        self.flux = np.zeros(self.rect_wave.shape)
        self.back = np.zeros(self.rect_wave.shape)
        self.spec_model = np.zeros(self.rect_spec.shape)
        self.back_model = np.zeros(self.rect_spec.shape)
        for i, v in enumerate(self.rect_wave):
            self.model.parameters = self.default_model_parms * 1.
            self.fixed_params()
            for pn in self.tinker_params:
                getattr(self.PSF, pn).value = np.interp(self.rect_wave[i],
                                                        self.dar_wave,
                                                        getattr(self,
                                                                'dar_' + pn))
            self.PSF.x_0.value += xoff
            self.PSF.y_0.value += yoff
            self.Back.parameters[0] = np.nanmedian(self.rect_spec[:, i])
            self.PSF.parameters[0] = np.nanpercentile(self.rect_spec[:, i], 99)
            fit = self.fitter(self.model, self.x[self.goodfibers],
                              self.y[self.goodfibers],
                              self.rect_spec[self.goodfibers, i])
            self.model = fit
            self.Back, self.PSF = self.model
            self.spec_model[:, i] = self.model(self.x, self.y)
            self.back_model[:, i] = self.Back(self.x, self.y)
            xo = self.PSF.parameters[1]
            yo = self.PSF.parameters[2]
            self.flux[i] = (self.PSF(xgrid.ravel()+xo,
                                     ygrid.ravel()+yo).sum() *
                            (boxsize_forintegral / boxsize_gridlength)**2 /
                            (0.295**2 * np.pi))
            self.back[i] = (self.Back(xo, yo) * self.flux[i] /
                            self.PSF.parameters[0])
