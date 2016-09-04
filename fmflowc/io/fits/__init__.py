# coding: utf-8

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

# dependent libraries
from astropy.io import fits

# sub modules/functions
from .aste.merge_data import fromaste


def open(*args, **kwargs):
    return fits.open(*args, **kwargs)


def getheader(*args, **kwargs):
    return fits.getheader(*args, **kwargs)


def getdata(*args, **kwargs):
    return fits.getdata(*args, **kwargs)
