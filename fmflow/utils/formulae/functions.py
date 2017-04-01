# coding: utf-8

# Python 3.x compatibility
from __future__ import absolute_import as _absolute_import
from __future__ import division as _division
from __future__ import print_function as _print_function

# the Python Package Index
import fmflow as fm
import numpy as np
import numpy.ma as ma

# imported items
__all__ = [
    'gaussian', 'lorentzian', 'pseudovoigt',
    'dgaussian', 'dlorentzian', 'dpseudovoigt', 'normaldist', 'mad'
]


def gaussian(x, x0=0.0, fwhm=1.0, ampl=1.0):
    return ampl * np.exp(-4*np.log(2)*((x-x0)/fwhm)**2)


def lorentzian(x, x0=0.0, fwhm=1.0, ampl=1.0):
    return ampl * (1+4*((x-x0)/fwhm)**2)**(-1)


def pseudovoigt(x, x0=0.0, fwhm=1.0, ampl=1.0, frac=0.5):
    return frac*gaussian(x, x0, fwhm, ampl) + (1-frac)*lorentzian(x, x0, fwhm, ampl)


def dgaussian(x, x0=0.0, fwhm=1.0, ampl=1.0):
    return -8*np.log(2)*(x-x0)/fwhm**2 * gaussian(x, x0, fwhm, ampl)


def dlorentzian(x, x0=0.0, fwhm=1.0, ampl=1.0):
    return -8*(x-x0)/fwhm**2 * lorentzian(x, x0, fwhm, ampl)


def dpseudovoigt(x, x0=0.0, fwhm=1.0, ampl=1.0, frac=0.5):
    return frac*dgaussian(x, x0, fwhm, ampl) + (1-frac)*dlorentzian(x, x0, fwhm, ampl)


def normaldist(x, mean=0.0, variance=1.0):
    return np.sqrt(2*np.pi*variance)**(-1) * np.exp(-(x-mean)**2/(2*variance))


def mad(array, axis=None, keepdims=False):
    """Compute the median absolute deviation (MAD) along the given axis.

    Args:
        array (array): An input array.
        axis (int, optional): Axis along which the MADs are computed.
            The default is to compute the MAD along a flattened version of the array.
        keepdims (bool, optional): If True, the axes which are reduced are left
            in the result as dimensions with size one.

    Returns:
        mad (array): A new array holding the result.

    """
    ad = ma.abs(array - ma.median(array, axis, keepdims=True))
    mad = ma.median(ad, axis, keepdims=keepdims)
    return mad
