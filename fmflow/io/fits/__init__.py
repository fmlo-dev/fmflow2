# coding: utf-8

from __future__ import absolute_import as _absolute_import
from __future__ import division as _division
from __future__ import print_function as _print_function

# dependent packages
from astropy.io import fits
open = fits.open
getdata = fits.getdata
getheader = fits.getheader

# submodules
from .array import *
from .aste import *

# imported items
__all__ = ['open', 'getdata', 'getheader']
__all__ += array.__all__
__all__ += aste.__all__