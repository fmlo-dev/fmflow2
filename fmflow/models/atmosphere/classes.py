# coding: utf-8

# Python 3.x compatibility
from __future__ import absolute_import as _absolute_import
from __future__ import division as _division
from __future__ import print_function as _print_function

# the Python standard library
import os
from subprocess import Popen, PIPE

# the Python Package Index
import yaml
import fmflow as fm
import numpy as np
from astropy import units as u
from astropy import constants as consts
from scipy.interpolate import interp1d
from scipy.ndimage import filters
from scipy.optimize import curve_fit

# imported items
__all__ = ['OzoneLines']

# constants
C          = consts.c.value # spped of light in vacuum
DIR_MODULE = os.path.dirname(os.path.realpath(__file__)) # module directory
AM_MODELS  = {
    'midaltitude': 'model_midaltitude.yaml',
    'chajnantor': 'model_chajnantor.yaml'
} # am models


class OzoneLines(object):
    def __init__(self, model='midaltitude', smooth=50):
        self.info = {
            'model': model,
            'smooth': smooth,
            'tuned': False
        }

        with open(os.path.join(DIR_MODULE, AM_MODELS[model])) as f:
            d = yaml.load(f)
            self.amc    = d['amc']
            self.layers = d['layers']

    def tune(self, fmarray):
        freq  = fm.model.getfreq(fmarray, 'GHz') # GHz
        step  = 1e-6 * fmarray.info['chanwidth'] # MHz

        fmin   = np.floor(np.min(freq))
        fmax   = np.ceil(np.max(freq))
        fstep  = float('{:.0e}'.format(0.5*step))
        params = {'fmin': fmin, 'fmax': fmax, 'fstep': fstep}

        mfreq = None
        mtaus, mTbs = [], []
        for i in range(len(self.layers)):
            frac = (i+1)/len(self.layers)
            fm.utils.progressbar(frac)

            params.update(**self.layers[i])
            amc = self.amc.format(**params)

            proc = Popen(['am', '-'], stdin=PIPE, stdout=PIPE, stderr=PIPE)
            stdout, stderr = proc.communicate(amc)
            output = np.loadtxt(stdout.split('\n'))

            if mfreq is None:
                mfreq = output[:,0]

            mtaus.append(output[:,1])
            mTbs.append(output[:,2])

        self.mfreq = mfreq
        self.mtaus = np.array(mtaus)
        self.mTbs  = np.array(mTbs)
        self.info['tuned'] = True

    def retrievefrom(self, fmarray, weights=None, mode='normal'):
        if not self.info['tuned']:
            raise fm.utils.FMFlowError('not tuned yet')

        freq = fm.model.getfreq(fmarray, 'GHz')
        spec = fm.model.getspec(fmarray, 'K', weights)
        vrad = np.mean(fmarray.vrad) # m/s
        frad = 1e-9 * fmarray.info['restfreq'] * vrad/C # GHz

        if mode == 'normal':
            model = self._fit(freq-frad, spec)
        elif mode == 'diff':
            model = self._dfit(freq-frad, spec)
        else:
            raise fm.utils.FMFlowError('invalid mode')

        fmmodel  = fm.zeros_like(fmarray)
        fmmodel_ = fm.demodulate(fmmodel)
        fmmodel_ += model
        fmmodel  = fm.modulate(fmmodel_)
        return fmmodel

    def subtractfrom(self, fmarray, weights=None, mode='normal'):
        if not self.info['tuned']:
            raise fm.utils.FMFlowError('not tuned yet')

        fmarray -= self.retrievefrom(fmarray, weights, mode)

    def _fit(self, freq, spec):
        Tbs = interp1d(self.mfreq, self.mTbs, axis=1)(freq)

        def func(freq, *coeffs):
            coeffs = np.asarray(coeffs)[:,np.newaxis]
            return np.sum(coeffs * Tbs, 0)

        p0 = np.zeros(len(Tbs))
        bs = (0.0, 1.0)
        popt, pcov = curve_fit(func, freq, spec, p0, bounds=bs)

        model = func(freq, *popt)
        if np.max(model) < 3.0*fm.utils.mad(spec):
            model[:] = 0.0

        return model

    def _dfit(self, freq, spec):
        Tbs  = interp1d(self.mfreq, self.mTbs, axis=1)(freq)
        dTbs = np.gradient(Tbs, axis=1)

        dspec = np.gradient(spec)
        dspec -= filters.gaussian_filter(dspec, self.info['smooth'])

        def func(freq, *coeffs):
            coeffs = np.asarray(coeffs)[:,np.newaxis]
            return np.sum(coeffs * dTbs, 0)

        p0 = np.zeros(len(dTbs))
        bs = (0.0, 1.0)
        popt, pcov = curve_fit(func, freq, dspec, p0, bounds=bs)

        dmodel = func(freq, *popt)
        if np.max(dmodel) < 3.0*fm.utils.mad(dspec):
            dmodel[:] = 0.0

        model = np.cumsum(dmodel)
        return model

