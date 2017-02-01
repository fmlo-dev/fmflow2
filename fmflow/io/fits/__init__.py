# coding: utf-8

# Python 3.x compatibility
from __future__ import absolute_import as __absolute_import
from __future__ import division as __division
from __future__ import print_function as __print_function

# the Python Package Index
from astropy.io import fits

# FMFlow submodules
from .array import *
from ._aste import *

# imported items
__all__ = ['open', 'getdata', 'getheader']
__all__ += array.__all__
__all__ += _aste.__all__


open = fits.open
getdata = fits.getdata
getheader = fits.getheader
