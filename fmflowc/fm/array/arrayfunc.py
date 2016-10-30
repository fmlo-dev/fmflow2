# coding: utf-8

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

# dependent libraries
import numpy as np
import numpy.ma as ma
from fmflowc.utils import exceptions as e

# sub modules/functions
from .array import FMData, FMArray


def array(array, fmch=None, coord=None, info=None):
    return FMArray(array, fmch, coord, info)

def asarray(array):
    return FMArray(array)

def asmaskedarray(array):
    return array.tomaskedarray()

def zeros(shape, dtype=float):
    return FMArray(np.zeros(shape, dtype))

def ones(shape, dtype=float):
    return FMArray(np.ones(shape, dtype))

def zeros_like(array, dtype=float, keepmeta=True):
    if keepmeta:
        return np.zeros_like(array, dtype)
    else:
        return FMArray(np.zeros_like(array.tomaskedarray(), dtype))

def ones_like(array, dtype=float, keepmeta=True):
    if keepmeta:
        return np.ones_like(array, dtype)
    else:
        return FMArray(np.ones_like(array.tomaskedarray(), dtype))

def concatenate(arrays):
    if not (type(arrays)==tuple or type(arrays)==list):
        raise e.FMflowError('arrays must be tuple or list of arrays')

    if len(arrays) > 2:
        array_0 = arrays[0]
        array_1 = concat(arrays[1:])
    else:
        array_0 = arrays[0]
        array_1 = arrays[1]
    
    array = ma.concatenate(map(asmaskedarray, [array_0, array_1]), 0)
    table = np.concatenate([array_0.table, array_1.table])
    info = array_0.info
    info.update(array_1.info)
    
    return FMArray.frommaskedarray(array, table, info)
