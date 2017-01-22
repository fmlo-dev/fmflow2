# coding: utf-8

# Python 3.x compatibility
from __future__ import absolute_import as __absolute_import
from __future__ import division as __division
from __future__ import print_function as __print_function

# the Python Package Index¬
import numpy as np
import numpy.ma as ma
import astropy.units as u
from fmflow import utils as ut

# imported items
__all__ = ['getfrequency', 'getspectrum', 'getfreq', 'getspec']


def getfrequency(array, unit='GHz', **kwargs):
    if array.ismodulated:
        raise ut.FMFlowError('array should be demodulated')

    info = array.info.copy()
    info.update(kwargs)
    rest = info['restfreq']
    step = info['chanwidth']
    fmindex = info['fmindex']

    start = rest - step * (0.5*(np.diff(fmindex)[0]-1)+fmindex[0])
    end   = start + step * array.shape[1]

    freq = np.arange(start, end, step) * u.Hz
    freq = freq.to(getattr(u, unit)).value
    return freq


def getspectrum(array, unit='K', weights=None):
    if array.ismodulated:
        raise ut.FMFlowError('array should be demodulated')

    spec = ma.average(array, 0, weights) * u.K
    spec = spec.to(getattr(u, unit)).value
    return spec


# aliases
getfreq = getfrequency
getspec = getspectrum
