# coding: utf-8

# Python 3.x compatibility
from __future__ import absolute_import as _absolute_import
from __future__ import division as _division
from __future__ import print_function as _print_function

# the Python Package Index
import numpy as np
import numpy.ma as ma
import fmflow as fm
from astropy import units as u

# imported items
__all__ = ['getfrequency', 'getspectrum', 'getfreq', 'getspec']


def getfrequency(fmarray, unit='GHz', **kwargs):
    if fmarray.ismodulated:
        fmarray = fm.demodulate(fmarray)

    info = fmarray.info.copy()
    info.update(kwargs)
    fmindex = info['fmindex']
    rest = info['restfreq']
    step = info['chanwidth']

    start = rest - step*(0.5*(np.diff(fmindex)[0]-1)+fmindex[0])
    end   = start + step*fmarray.shape[1]

    freq = np.arange(start, end, step) * u.Hz
    freq = freq.to(getattr(u, unit)).value
    return freq


def getspectrum(fmarray, unit='K', weights=None):
    if fmarray.ismodulated:
        fmarray = fm.demodulate(fmarray)

    if weights is not None:
        if weights.ismodulated:
            weights = fm.demodulate(weights)

    spec = ma.average(fmarray, 0, weights) * u.K
    spec = spec.to(getattr(u, unit)).value
    return spec


# aliases
getfreq = getfrequency
getspec = getspectrum
