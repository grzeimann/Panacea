""" Panacea - telluricabs.py


1) Build telluric absorption model

.. moduleauthor:: Greg Zeimann <gregz@astro.as.utexas.edu>

"""

import telfit
import numpy as np
import matplotlib.pyplot as plt
from astropy.convolution import Gaussian1DKernel, convolve
from scipy.interpolate import splrep, splev

try:
    import seaborn as sns
    sns.set_context('notebook', font_scale=1.5)
    sns.set_style('whitegrid')
except ImportError:
    pass


class TelluricAbs:
    ''' Telluric absorption through telfit.py '''
    def __init__(self, wave, spec, RH, T, P, ZD):
        '''
        Build the fitter and data set for telluric absorption fitting

        Parameters:
        wave : array
            wavelength
        spec : array
            A-star (or standard) spectrum for telluric model
        RH : float
            Relative Humidity
        T : float
            Temperature ()
        P : float
            Pressure (mm Hg)
        ZD : float
            Zenith Distance (angle in degrees)
        '''
        self.wave = wave
        self.spec = spec
        self.data = telfit.DataStructures.xypoint(x=self.wave/10.,
                                                  y=self.spec)
        self.fitter = telfit.TelluricFitter()
        self.fitter.SetObservatory('McDonald')
        self.modeler = telfit.Modeler()
        self.RH = RH
        self.T = T
        self.P = P
        self.ZD = ZD
        self.fitter.AdjustValue(dict(wavestart=self.data.x[0]-5.0,
                                     waveend=self.data.x[-1]+5.0))
        self.set_weather_values()
        self.set_resolution()

    def get_model(self, lowwave, highwave, wave, humidity=50., o2=180000.,
                  resolution=0.2):
        lowfreq = 1e8 / highwave
        highfreq = 1e8 / lowwave
        m = self.modeler.MakeModel(lowfreq=lowfreq, highfreq=highfreq,
                                   humidity=humidity, o2=o2)
        m.x = m.x * 10.
        tck = splrep(m.x, m.y)
        newwave = np.arange(m.x.min(), m.x.max(), np.min(np.diff(m.x)))
        ntrans = splev(newwave, tck)
        kernel = resolution / np.min(np.diff(m.x))
        G = Gaussian1DKernel(kernel)
        smooth = convolve(ntrans, G)
        tck = splrep(newwave, smooth)
        lowres = splev(wave, tck)
        return lowres

    def set_weather_values(self):
        ''' set weather data for different observations '''
        self.weather_dict = dict(temperature=(self.T + 273.15),
                                 pressure=self.P, angle=self.ZD)
        self.fitter.AdjustValue(self.weather_dict)

    def set_resolution(self, fixed=True, res_in_pix=4.5):
        resolution = (np.median(self.wave) /
                      (res_in_pix * np.median(np.diff(self.wave))))
        self.fitter.AdjustValue(dict(resolution=resolution))
        self.fitter.SetBounds(dict(resolution=(0.7*resolution,
                                               1.3*resolution)))
        if not fixed:
            self.fitter.FitVariabe(dict(resolution=resolution))

    def fit_telluric_abs(self):
        self.fitter.FitVariable(dict(h2o=self.RH))
        self.fitter.SetBounds(dict(h2o=(1, 99)))
        self.source, self.model = self.fitter.Fit(data=self.data,
                                                  resolution_fit_mode='gauss',
                                                  adjust_wave='data',
                                                  continuum_fit_order=7,
                                                  wavelength_fit_order=3,
                                                  air_wave=True,
                                                  fit_source=True)

    def plot_telluric_fit(self, outname):
        newdata = self.fitter.data
        fig, (top, bottom) = plt.subplots(2, 1, sharex=True, figsize=(13, 10),
                                          gridspec_kw=dict(height_ratios=(3,
                                                                          1)))

        top.plot(newdata.x, newdata.y/newdata.cont, 'k-', label='Data')
        top.plot(self.model.x, self.model.y, 'r-', label='Telluric Model')
        bottom.plot(newdata.x, newdata.y/newdata.cont / self.model.y, 'k-')

        bottom.set_xlabel('Wavelength (nm)')
        bottom.set_ylabel('Residuals')
        top.set_ylabel('Normalized Flux')

        top.legend(loc='best', fancybox=True)

        bottom.set_ylim((0.9, 1.1))
        fig.savefig(outname)
        plt.close(fig)
