# coding: utf-8

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

# dependent libraries¬
import numpy as np
import numpy.ma as ma
import astropy.units as u
from fmflow.utils import exceptions as e


def getfrequency(array, unit='GHz', **kwargs):
    if array.ismodulated:
        raise e.FMflowError('array should be demodulated')

    info = array.info.copy()
    info.update(kwargs)
    fmid = info['fmid']
    rest = info['restfreq']
    step = info['chanwidth']
    
    start = rest - step * (0.5*(fmid[1]-1)+fmid[0])
    end   = start + step * array.shape[1]
    
    freq = np.arange(start, end, step) * u.Hz
    freq = freq.to(getattr(u, unit)).value
    return freq

def getspectrum(array, unit='K', weights=None):
    if array.ismodulated:
        raise e.FMflowError('array should be demodulated')

    spec = ma.average(array, 0, weights) * u.K
    spec = spec.to(getattr(u, unit)).value
    return spec
