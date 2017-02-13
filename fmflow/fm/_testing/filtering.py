# coding: utf-8

# Python 3.x compatibility
from __future__ import absolute_import as __absolute_import
from __future__ import division as __division
from __future__ import print_function as __print_function

# the Python standard library
from functools import partial

# the Python Package Index
import numpy as np
from scipy import signal
from scipy import interpolate
from scipy.ndimage import filters

# FMFlow submodules
from .._array import *

# imported items
__all__ = ['hmedfilt', 'lmeds_spline', 'fme_gaussian_filter']


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


def fme_gaussian_filter(fmarray_in, sigma=10, m=10, dz=None):
    """Apply Gaussian filter to an array.

    Args:
        fmarray_in (array)
        sigma (float)
        dz (float)
        m (int)

    Returns:
        fmarray_out (array)

    """
    fmarray_in = np.asarray(fmarray_in).copy()
    x, z = np.arange(len(fmarray_in)), fmarray_in
    sigma = float(sigma)
    m = int(m)

    if dz is None:
        dz = np.ptp(z) / (len(z)/5.0)
    else:
        dz = float(dz)

    g = lambda x: np.exp(-x**2/(2*sigma**2)) / np.sqrt(2*np.pi*sigma**2)
    c = partial(filters.convolve1d, weights=np.ones(m), axis=0)

    rho = z / dz
    binrho = np.arange(int(min(rho))-m, int(max(rho))+m)
    binidx = np.digitize(rho, binrho)

    w = np.zeros([len(binrho), len(x)])
    for i in range(len(x)):
        w[binidx[i]] += g(x-x[i])

    wc = c(c(c(w)))
    q = lambda i: wc[np.argmax(wc,0)+i,x]
    m = lambda x: (x[np.argmax(wc,0)]+x[np.argmax(wc,0)-1])/2.0

    fmarray_out = dz * (m(binrho) + (q(-1)-q(0)) / (q(-1)-2.0*q(0)+q(1)))
    return fmarray_out
