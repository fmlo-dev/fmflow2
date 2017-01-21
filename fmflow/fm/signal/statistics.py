# coding: utf-8

from __future__ import absolute_import as __absolute_import
from __future__ import division as __division
from __future__ import print_function as __print_function

# dependent packages
import numpy as np
import numpy.ma as ma

# imported items
__all__ = ['mad']


def mad(array, axis=None):
    """Compute the median absolute deviation (MAD) along the given axis.

    Args
    ----

    - array (array-like): An input array.
    - axis (int, None): Axis along which the MADs are computed.
      The default is to compute the MAD along a flattened version of the array.
    
    Returns
    -------

    - mad (array): 
    """
    if axis == 1:
        mad = ma.median(ma.abs(array - ma.median(array, 1)[:,np.newaxis]), 1)
    else:
        mad = ma.median(ma.abs(array - ma.median(array, axis)), axis)
    
    return mad