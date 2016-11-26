# coding: utf-8

from __future__ import absolute_import as _absolute_import
from __future__ import division as _division
from __future__ import print_function as _print_function

# dependent packages
import numpy as np
from scipy import signal
from scipy import interpolate

# submodules
from ..core import *

# imported items
__all__ = ['hmedfilt', 'lmeds_spline']


@fmfunc
@timechunk
def hmedfilt(array_in, kernel=3):
    array_out = signal.medfilt(array_in, (1, kernel))
    return array_out


@fmfunc
@timechunk
def lmeds_spline(array_in, nsample=50, niter=1000, **kwargs):
    if array_in.ndim == 1:
        array_in = array_in[np.newaxis]

    freq = np.arange(array_in.shape[1])
    array_out = np.zeros_like(array_in)

    for i in range(array_in.shape[0]):
        data = array_in[i]
        models = []
        medians = []
        j = 0
        while j < niter:
            index = np.random.choice(freq, nsample)
            index.sort()

            spline = interpolate.UnivariateSpline(index, data[index], s=0, **kwargs)
            model = spline(freq)
            if np.all(np.isnan(model)):
                continue

            medians.append(np.median((data-model)**2))
            models.append(model)
            j += 1
            
        array_out[i] = models[np.argmin(medians)]
    
    array_out = np.squeeze(array_out)
    return array_out
