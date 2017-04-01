# coding: utf-8

# Python 3.x compatibility
from __future__ import absolute_import as _absolute_import
from __future__ import division as _division
from __future__ import print_function as _print_function

# the Python Package Index
import fmflow as fm
import numpy as np
from astropy import units as u
from scipy.ndimage import filters
from scipy.optimize import curve_fit

# imported items
__all__ = ['GaussianLines', 'GaussianFilter']


class GaussianLines(object):
    def __init__(self, smooth=100, threshold=10.0, gain=0.5):
        self.info = {
            'gain': gain,
            'smooth': smooth,
            'threshold': threshold,
        }

    def retrievefrom(self, fmarray, weights=None, mode='normal'):
        freq = fm.models.getfreq(fmarray, 'GHz')
        spec = fm.models.getspec(fmarray, 'K', weights)
        fwhm0 = 5.0 * (1e-9*fmarray.info['chanwidth']) # GHz

        if mode == 'normal':
            model = self._fit(freq, spec, fwhm0)
        elif mode == 'diff':
            model = self._dfit(freq, spec, fwhm0)
        else:
            raise fm.utils.FMFlowError('invalid mode')

        fmmodel  = fm.zeros_like(fmarray)
        fmmodel_ = fm.demodulate(fmmodel)
        fmmodel_ += model
        fmmodel  = fm.modulate(fmmodel_)
        return fmmodel

    def subtractfrom(self, fmarray, weights=None, mode='normal'):
        fmarray -= self.retrievefrom(fmarray, weights, mode)

    @staticmethod
    def _sn(spec):
        return spec / fm.utils.mad(spec)

    def _fit(self, freq, spec, fwhm0):
        model = np.zeros_like(spec)
        resid = spec.copy()

        while np.max(self._sn(resid)) > self.info['threshold']:
            cent0 = freq[np.argmax(resid)]
            ampl0 = np.max(resid)

            p0 = [cent0, fwhm0, ampl0]
            bs = ([np.min(freq), 0.0, 0.0], [np.max(freq), 10.0*fwhm0, 10.0*ampl0])
            popt, pcov = curve_fit(fm.utils.gaussian, freq, resid, p0, bounds=bs)

            model += self.info['gain'] * fm.utils.gaussian(freq, *popt)
            resid -= self.info['gain'] * fm.utils.gaussian(freq, *popt)

        return model

    def _dfit(self, freq, spec, fwhm0):
        dspec  = np.gradient(spec)
        dspec -= filters.gaussian_filter(dspec, self.info['smooth'])
        dmodel = np.zeros_like(dspec)
        dresid = dspec.copy()

        while np.max(self._sn(dresid)) > self.info['threshold']:
            fmax = freq[np.argmax(dresid)]
            smax = np.max(dresid)

            cent0 = fmax + fwhm0/np.sqrt(8*np.log(2))
            ampl0 = smax * np.exp(0.5) * fwhm0/np.sqrt(8*np.log(2))

            p0 = [cent0, fwhm0, ampl0]
            bs = ([np.min(freq), 0.0, 0.0], [np.max(freq), 10.0*fwhm0, 10.0*ampl0])
            popt, pcov = curve_fit(fm.utils.dgaussian, freq, dresid, p0, bounds=bs)

            dmodel += self.info['gain'] * fm.utils.dgaussian(freq, *popt)
            dresid -= self.info['gain'] * fm.utils.dgaussian(freq, *popt)

        model = np.cumsum(dmodel)
        return model


class GaussianFilter(object):
    def __init__(self, smooth=100):
        self.info = {'smooth': smooth}

    def retrievefrom(self, fmarray, weights=None, mode='normal'):
        spec = fm.models.getspec(fmarray, 'K', weights)

        if mode == 'normal':
            model = self._fit(spec)
        elif mode == 'diff':
            model = self._dfit(spec)
        else:
            raise fm.utils.FMFlowError('invalid mode')

        fmmodel  = fm.zeros_like(fmarray)
        fmmodel_ = fm.demodulate(fmmodel)
        fmmodel_ += model
        fmmodel  = fm.modulate(fmmodel_)
        return fmmodel

    def subtractfrom(self, fmarray, weights=None, mode='normal'):
        fmarray -= self.retrievefrom(fmarray, weights, mode)

    def _fit(self, spec):
        return fm.utils.fmgf(spec, self.info['smooth'])

    def _dfit(self, spec):
        dspec  = np.gradient(spec)
        return fm.utils.fmgf(dspec, self.info['smooth'])
