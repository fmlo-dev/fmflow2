# coding: utf-8

# Python 3.x compatibility
from __future__ import absolute_import as _absolute_import
from __future__ import division as _division
from __future__ import print_function as _print_function

# the Python standard library
from functools import partial

# the Python Package Index
import fmflow as fm
import numpy as np
from scipy.ndimage import filters

# imported items
__all__ = ['fmgf']


def fmgf(array, sigma):
    x, y = np.arange(len(array)), array.copy()
    yg = filters.gaussian_filter(y, sigma)
    y -= yg

    # digitizing
    m = 101
    dy = 6.0*fm.utils.mad(y) / m
    ybin = np.arange(np.min(y)-5*dy, np.max(y)+5*dy+dy, dy)
    z = np.zeros([len(ybin), len(x)])
    z[np.digitize(y, ybin), x] = 1.0

    # filtering
    g = partial(filters.gaussian_filter, sigma=(0,sigma))
    c = partial(filters.convolve1d, weights=np.ones(m), axis=0)    
    zf = c(c(c(g(z))))

    # estimates
    ym1, y0, yp1 = [ybin[np.argmax(zf,0)+i] for i in (-1,0,1)]
    zm1, z0, zp1 = [zf[np.argmax(zf,0)+i, x] for i in (-1,0,1)]
    t = (zm1-z0) / (zm1-2*z0+zp1)

    filtered = yg + ((1-t)**2)*ym1 + (2*t*(1-t))*y0 + (t**2)*yp1
    return filtered
