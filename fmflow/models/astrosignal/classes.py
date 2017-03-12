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
__all__ = ['GaussianModel']


class GaussianModel(object):
    def __init__(self, threshold=5.0, smooth=50, gain=0.75):
        self.info = {
            'gain': gain,
            'smooth': smooth,
            'threshold': threshold,
            'tuned': False
        }

    def tune(self, fmarray):
        self.info['fwhm0'] = 5.0 * (1e-9*fmarray.info['chanwidth']) # GHz
        self.info['tuned'] = True

    def retrievefrom(self, fmarray, weights=None, mode='normal'):
        if not self.info['tuned']:
            raise fm.utils.FMFlowError('not tuned yet')

        freq = fm.model.getfreq(fmarray, 'GHz')
        spec = fm.model.getspec(fmarray, 'K', weights)

    def subtractfrom(self, fmarray, weights=None, mode='normal'):
        if not self.info['tuned']:
            raise fm.utils.FMFlowError('not tuned yet')

        fmarray -= self.retrievefrom(fmarray, weights, mode)

    def _snr(self, spec):
        return spec / fm.utils.mad(spec)

    def _fit(self, freq, spec):
        model = np.zeros_like(spec)
        resid = spec.copy()

        while np.max(self._snr(resid)) > self.info['threshold']:
            cent0 = freq[np.argmax(resid)]
            fwhm0 = self.info['fwhm0']
            ampl0 = np.max(resid)

            p0 = [cent0, fwhm0, ampl0]
            bs = (0.0, [10.0*np.abs(ampl0), np.max(freq), 10.0*fwhm0])
            popt, pcov = curve_fit(fm.utils.gaussian, freq, resid, p0, bounds=bs)

            model += self.info['gain'] * fm.utils.gaussian(freq, *popt)
            resid -= self.info['gain'] * fm.utils.gaussian(freq, *popt)

        return model

    def _dfit(self, freq, spec):
        dspec  = np.gradient(spec)
        dspec -= filters.gaussian_filter(dspec, self.info['smooth'])
        dmodel = np.zeros_like(dspec)
        dresid = dspec.copy()

        while np.max(self._snr(dresid)) > self.info['threshold']:
            fmax = freq[np.argmax(dresid)]
            smax = np.max(dresid)

            cent0 = fmax + fwhm0/np.sqrt(8*np.log(2))
            fwhm0 = self.info['fwhm0']
            ampl0 = smax * np.exp(0.5) * fwhm0/np.sqrt(8*np.log(2))

            p0 = [cent0, fwhm0, ampl0]
            bs = (0.0, [10.0*np.abs(ampl0), np.max(freq), 10.0*fwhm0])
            popt, pcov = curve_fit(fm.utils.dgaussian, freq, dresid, p0, bounds=bs)

            dmodel += self.info['gain'] * fm.utils.dgaussian(freq, *popt)
            dresid -= self.info['gain'] * fm.utils.dgaussian(freq, *popt)

        model = np.cumsum(dmodel)
        return model

