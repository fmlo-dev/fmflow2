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
def hmedfilt(fmarray_in, kernel=3):
    fmarray_out = signal.medfilt(fmarray_in, (1, kernel))
    return fmarray_out


@fmfunc
@timechunk
def lmeds_spline(fmarray_in, nsample=50, niter=1000, **kwargs):
    if fmarray_in.ndim == 1:
        fmarray_in = fmarray_in[np.newaxis]

    freq = np.arange(fmarray_in.shape[1])
    fmarray_out = np.zeros_like(fmarray_in)

    for i in range(fmarray_in.shape[0]):
        data = fmarray_in[i]
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
            
        fmarray_out[i] = models[np.argmin(medians)]
    
    fmarray_out = np.squeeze(fmarray_out)
    return fmarray_out
