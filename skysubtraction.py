""" Panacea - skysubtraction.py


1) Build telluric absorption model

.. moduleauthor:: Greg Zeimann <gregz@astro.as.utexas.edu>

"""

import matplotlib.pyplot as plt
import numpy as np
from fiber_utils import find_maxima, bspline_x0
from utils import biweight_location, matrixCheby2D_7
from scipy.signal import savgol_filter, medfilt
from scipy.ndimage.filters import percentile_filter
from astropy.stats import biweight_midvariance
from astropy.convolution import Gaussian1DKernel, convolve
from scipy.interpolate import splrep, splev
from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import RANSACRegressor
from sklearn.pipeline import Pipeline


class Sky:
    ''' Wrapper for reduction routines with processed data, multi*.fits '''
    def __init__(self, wave, spec, rect_wave, rect_spec, trace,
                 goodfibers, skyline_file=None):
        '''
        DESCRIPTION

        Parameters:
        wave : 2-d array
            wavelength for each fiber
        spec : 2-d array
            fiber to fiber corrected spectrum for each fiber
        rect_wave : 1-d array
            Common wavelength for rectified spectra
        rect_spec : 2-d array
            Rectified spectra corrected for fiber to fiber
        '''
        self.wave = wave * 1.
        self.spec = spec * 1.
        self.trace = trace * 1.
        self.goodfibers = goodfibers
        self.rect_wave = rect_wave * 1.
        self.rect_spec = rect_spec * 1.
        self.skyline_file = skyline_file

    def convert_vac_to_air(self):
        s2 = (1e4 / self.skyline[:, 0])**2
        n = (1 + 0.0000834254 + 0.02406147 / (130 - s2) + 0.00015998 /
             (38.9 - s2))
        self.skyline[:, 0] = self.skyline[:, 0] / n

    def get_skyline_file(self):
        V = np.loadtxt(self.skyline_file)
        V[:, 0] = V[:, 0] * 1e4
        sel = np.where((V[:, 0] > self.rect_wave.min()) *
                       (V[:, 0] < self.rect_wave.max()))[0]
        self.skyline = V[sel, :]
        self.convert_vac_to_air()

    def fit_continuum_sky(self, wv, sky, fil_len=95, func=np.array):
        skym_s = 1. * sky
        sky_sm = savgol_filter(skym_s, fil_len, 1)
        for i in np.arange(5):
            mad = np.median(np.abs(sky - sky_sm))
            outlier = func(sky - sky_sm) > 1.5 * mad
            skym_s = 1.*sky
            skym_s[outlier] = np.interp(wv[outlier], wv[~outlier],
                                        sky_sm[~outlier])
            sky_sm = savgol_filter(skym_s, fil_len, 1)
        return sky_sm

    def make_skyline_model(self, skylines, norm, kernel_size=2.1):
        G = Gaussian1DKernel(kernel_size)
        pixsize = self.rect_wave[1] - self.rect_wave[0]
        skymodel = np.zeros(self.rect_wave.shape)
        for line in skylines:
            x0 = int(np.argmin(np.abs(self.rect_wave - line[0])))
            xl = x0 - len(G.array) / 2 - 3
            xh = x0 + len(G.array) / 2 + 4
            xl = np.max([0, xl])
            xh = np.min([len(self.rect_wave), xh])
            xlg = line[0] - (len(G.array) / 2) * pixsize
            xhg = line[0] + (len(G.array) / 2 + 0.5) * pixsize
            fac = line[1] * norm
            skymodel[xl:xh] += fac * np.interp(self.rect_wave[xl:xh],
                                               np.arange(xlg, xhg, pixsize),
                                               G.array, left=0.0, right=0.0)
        return skymodel

    def make_sky_data(self, lowfib, highfib, perc=50):
        sky = np.nanpercentile(self.rect_spec[lowfib:highfib, :], perc, axis=0)
        return sky

    def wavelength_from_sky(self, window_size=8, order=3):
        num = self.rect_spec.shape[0]
        mastersky = self.make_sky_data(0, num)
        lowsky = percentile_filter(mastersky, 3, 40)
        sky_cont = self.fit_continuum_sky(self.rect_wave, lowsky,
                                          func=np.array)
        mad = np.sqrt(np.median(sky_cont))
        pmr, pmh = find_maxima(self.rect_wave, mastersky, y_window=30,
                               interp_window=20, repeat_length=10)
        sel = (pmh - np.interp(pmr, self.rect_wave, sky_cont)) > (5. * mad)
        pmr, pmh = (pmr[sel], pmh[sel])
        self.A = -999.*np.ones((num, len(pmr)))
        self.skyline_loc = pmr
        xn = np.arange(self.rect_wave.shape[0])
        xi = np.arange(self.wave.shape[1])
        for i in np.arange(0, num, window_size):
            if i < num/2:
                mx = np.min([num/2, i + window_size + 1])
                mn = np.max([0, i - window_size])
            else:
                mx = np.min([num, i + window_size + 1])
                mn = np.max([num/2, i - window_size])
            fit_spec = self.make_sky_data(mn, mx)
            pr, ph = find_maxima(xn, fit_spec, y_window=30,
                                 interp_window=20, repeat_length=10)
            wv = np.interp(pr, xn, self.rect_wave)
            if len(wv):
                for j, p in enumerate(pmr):
                    if np.min(np.abs(wv - p)) < 5.:
                        ind = np.argmin(np.abs(wv - p))
                        self.A[i, j] = np.interp(wv[ind], self.wave[i, :], xi)
                    else:
                        self.A[i, j] = -999.
        xp, yp, wp = ([], [], [])
        n = self.wave.shape[1] - 1.
        for i in np.arange(self.A.shape[1]):
            sel = np.where(self.A[:, i] > -999.)[0]
            try:
                xloc = int(np.median(self.A[sel, i]))
            except:
                print(np.median(self.A[sel, i]), sel)    
            yp.append(self.trace[sel, xloc])
            xp.append(self.A[sel, i])
            wp.append(np.ones((len(sel),))*self.skyline_loc[i])
        x, y, w = [np.hstack(i) for i in [xp, yp, wp]]
        self.xp, self.yp, self.wp = (x, y, w)
        V = matrixCheby2D_7(x / n, y / n)
        w_sol = np.linalg.lstsq(V, w)[0]
        find, xind = np.indices(self.trace.shape)
        bigV = matrixCheby2D_7(xind.ravel() / n,
                               self.trace.ravel() / n)
        newwave = np.dot(bigV, w_sol)
        newwave = newwave.reshape(self.trace.shape)

        return newwave

    def adjust_wave(self, wave, buff=8, order=3):
        wv = wave * 1.
        trace = self.trace * 1.
        for i in np.arange(trace.shape[1]):
            wv[:, i] = self.robust_poly_fit(trace[:, i], wave[:, i],
                                            order=order)
        return wv

    def robust_poly_fit(self, x, y, order=3):
        model = Pipeline([('poly', PolynomialFeatures(degree=order)),
                          ('linear', RANSACRegressor())])
        model = model.fit(x[:, np.newaxis], y)
        return model.predict(x[:, np.newaxis])

    def iter_poly_fit(self, x, y, order=3, filt_len=7):
        s = savgol_filter(y, filt_len, 1)
        m = medfilt(y, filt_len)
        std = np.median(np.abs(s-m))
        inlier = np.where(np.abs(s - m) < (3*std))[0]
        p0 = np.polyfit(x[inlier], y[inlier], 1)
        ym = np.polyval(p0, x)
        std = np.percentile(np.abs(y-ym), 30)
        inlier = np.where(np.abs(s - m) < (3*std))[0]
        p0 = np.polyfit(x[inlier], y[inlier], order)
        return p0

    def make_sky_data2(self, wavelength, spectrum, fac=10):
        mn = wavelength.min()
        mx = wavelength.max()
        N, D = spectrum.shape
        xs = np.linspace(0, 1, D*fac)
        A = np.zeros((len(xs), len(self.goodfibers)))
        i = 0
        if self.goodfibers is None:
            goodfibers = np.arange(N)
        else:
            goodfibers = self.goodfibers
        for j, i in enumerate(goodfibers):
            y = spectrum[i, :]
            xp = np.interp(wavelength[i, :],
                           np.linspace(mn, mx, D*fac),
                           xs, left=0.0, right=0.0)
            tck = splrep(xp, y)
            A[:, j] = splev(xs, tck)
        ys = biweight_location(A, axis=(1,))
        masterwave = np.linspace(mn, mx, D*fac)
        B, c = bspline_x0(masterwave, nknots=D)
        sol = np.linalg.lstsq(c, ys, rcond=None)[0]
        mastersky = np.dot(c, sol)
        return masterwave, mastersky
    
    def subtract_sky(self, order=1, mask=None):
        if mask is None:
            usefibers = self.goodfibers
        else:
            usefibers = np.setdiff1d(self.goodfibers, mask)
        x0 = np.arange(self.rect_spec.shape[0])
        x = x0[usefibers]
        self.sky = np.zeros(self.rect_spec.shape)
        self.rect_spec[np.isnan(self.rect_spec)] = 0.0
        for i in np.arange(self.rect_spec.shape[1]):
            y = self.rect_spec[usefibers, i]
            self.sky[self.goodfibers, i] = np.polyval(np.polyfit(x, y, order),
                                                      x0[self.goodfibers])
        self.sky_subtracted = self.rect_spec - self.sky

            
