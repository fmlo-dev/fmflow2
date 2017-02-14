# coding: utf-8

"""Module for functions of basic fmarray operations.

fmarray is created and operated with the functions defined in this module.
These functions can be used as fm.<functions> like::

    >>> import fmflow as fm
    >>> fmarray = fm.array(array, table, info) # create an fmarray
    >>> fdarray = fm.demodulate(fmarray) # demodulate an fmarray

For more information, see (https://github.com/snoopython/fmflow/wiki/fmarray).

"""

# Python 3.x compatibility
from __future__ import absolute_import as _absolute_import
from __future__ import division as _division
from __future__ import print_function as _print_function

# the Python standard library
import uuid

# the Python Package Index
import astropy.units as u
import numpy as np
import numpy.ma as ma
import fmflow as fm

# imported items
__all__ = [
    'array', 'asarray', 'toarray', 'tomaskedarray', 'getarray',
    'demodulate', 'modulate', 'zeros', 'ones', 'zeros_like', 'ones_like',
    'concatenate', 'save', 'mad', 'median', 'load',
    'getfrequency', 'getspectrum', 'getfreq', 'getspec',
]


def array(array, fmch=None, coord=None, info=None):
    """Create a modulated fmarray from each argument.

    Args:
        array (masked array): A timestream array (mask is optional).
        fmch (array): Array of modulation frequencies in units of channel.
        coord (array): Array of observed coordinates in units of degrees.
        info (dict): Information of the array, observation, etc.

    Returns:
        fmarray (FMArray): A modulated fmarray.

    """
    fmarray = fm.FMArray.fromeach(array, fmch, coord, info)
    return fmarray


def asarray(array):
    """Convert the input array to a modulated fmarray.

    Args:
        array (array-like): An input array-like object.

    Returns:
        fmarray (FMArray): An output modulated fmarray.

    """
    if type(array) == fm.FMArray:
        fmarray = array
    else:
        fmarray = fm.FMArray(array)

    return fmarray


def toarray(fmarray):
    """Convert the fmarray to a NumPy ndarray.

    Args:
        fmarray (FMArray): An input fmarray

    Returns:
        array (array): An output NumPy ndarray.

    """
    array = fmarray.toarray()
    return array


def tomaskedarray(fmarray):
    """Convert the fmarray to a NumPy masked array.

    Args:
        fmarray (FMArray): An input fmarray

    Returns:
        array (masked array): An output NumPy masked array.

    """
    array = fmarray.tomaskedarray()
    return array


def getarray(fitsname, arrayid, scantype):
    """Create a modulated fmarray from a FMFITS.

    This functions is equivalent to the fmflow.io.fits.getarray function
    for the convenience (users can use it w/o importing fmflow.io.fits).

    Args:
        fitsname (str): File name of a FMFITS.
        arrayid (str): An array ID with which the output fmarray is created.
        scantype (str): A scan type with which the output fmarray is created.

    Returns:
        fmarray (FMArray): An output fmarray of the spacified array ID and scan type.

    """
    import fmflow as fm
    fmarray = fm.io.fits.getarray(fitsname, arrayid, scantype)
    return fmarray


def demodulate(fmarray_in, reverse=False):
    """Create a demodulated fmarray from the modulated one.

    This function is only available when the input fmarray is modulated.

    Args:
        fmarray_in (FMArray): An input modulated fmarray.
        reverse (bool): If True, the input fmarray is reverse-demodulated
            (i.e. fmch * -1 is used for demodulation). Default is False.

    Returns:
        fmarray_out (FMArray): An output demodulated fmarray.

    """
    fmarray_out = fmarray_in.demodulate(reverse)
    return fmarray_out


def modulate(fmarray_in):
    """Create an modulated fmarray from the demodulated one.

    This function is only available when the input fmarray is demodulated.

    Args:
        fmarray_in (FMArray): An input demodulated fmarray.

    Returns:
        fmarray_out (FMArray): An output modulated fmarray.

    """
    fmarray_out = fmarray_in.modulate()
    return fmarray_out


def zeros(shape, dtype=float, order='C', **kwargs):
    """Create a modulated fmarray of given shape and type, filled with zeros.

    Args:
        shape (int or sequence of ints): Shape of the fmarray.
        dtype (data-type): The desired data-type for the fmarray.
        order ('C' or 'F'): Whether to store multidimensional data
            in C- or Fortran-contiguous (row- or column-wise) order in memory.
        kwargs (optional): Other arguments of fmarray (fmch, coord, and info).

    Returns:
        fmarray (FMArray): An output modulated fmarray of zeros.

    """
    fmarray = fm.FMArray.fromeach(np.zeros(shape, dtype, order), **kwargs)
    return fmarray


def ones(shape, dtype=float, order='C', **kwargs):
    """Create a modulated fmarray of given shape and type, filled with ones.

    Args:
        shape (int or sequence of ints): Shape of the fmarray.
        dtype (data-type): The desired data-type for the fmarray.
        order ('C' or 'F'): Whether to store multidimensional data
            in C- or Fortran-contiguous (row- or column-wise) order in memory.
        kwargs (optional): Other arguments of fmarray (fmch, coord, and info).

    Returns:
        fmarray (FMArray): An output modulated fmarray of ones.

    """
    fmarray = fm.FMArray.fromeach(np.ones(shape, dtype, order), **kwargs)
    return fmarray


def zeros_like(array, dtype=None, order='K', keepmeta=True):
    """Create a fmarray of zeros with the same shape and type as an input array.

    Args:
        array (array-like): The shape and data-type of it define
            these same attributes of the output fmarray.
        dtype (data-type): If spacified, this function overrides
            the data-type of the output fmarray.
        order ('C', 'F', 'A', or 'K'):  If spacified, this function overrides
            the memory layout of the result. 'C' means C-order, 'F' means F-order,
            'A' means 'F' if `a` is Fortran contiguous, 'C' otherwise.
            'K' means match the layout of `a` as closely as possible.
        keepmeta (bool): Whether table and info are kept in the output fmarray
            if the input array is fmarray. Default is True.

    Returns:
        fmarray (FMArray): An output fmarray of zeros.

    """
    fmarray = np.zeros_like(array, dtype, order, subok=True)

    if not keepmeta:
        fmarray = fm.FMArray(fmarray)

    return fmarray


def ones_like(array, dtype=None, order='K', keepmeta=True):
    """Create a fmarray of ones with the same shape and type as an input array.

    Args:
        array (array-like): The shape and data-type of it define
            these same attributes of the output fmarray.
        dtype (data-type): If spacified, this function overrides
            the data-type of the output fmarray.
        order ('C', 'F', 'A', or 'K'):  If spacified, this function overrides
            the memory layout of the result. 'C' means C-order, 'F' means F-order,
            'A' means 'F' if `a` is Fortran contiguous, 'C' otherwise.
            'K' means match the layout of `a` as closely as possible.
        keepmeta (bool): Whether table and info are kept in the output fmarray
            if the input array is fmarray. Default is True.

    Returns:
        fmarray (FMArray): An output fmarray of ones.

    """
    fmarray = np.ones_like(array, dtype, order, subok=True)

    if not keepmeta:
        fmarray = fm.FMArray(fmarray)

    return fmarray


def concatenate(fmarray_ins):
    """Join a sequence of fmarrays along time axis.

    This function joins a sequence of fmarrays along time axis.
    All input fmarrays must have the same width (i.e. same shape[1]).
    Tables of fmarrays are also joined and infos are updated one after another.

    Args:
        fmarray_ins (sequence): An input sequence of fmarrays

    Returns:
        fmarray_out (FMArray): An output concatenated fmarray.

    """
    if type(fmarray_ins) not in (tuple, list):
        import fmflow as fm
        raise fm.utils.FMFlowError('fmarray_ins must be sequence of fmarrays.')

    if len(fmarray_ins) > 2:
        fmarray_0 = fmarray_ins[0]
        fmarray_1 = concatenate(fmarray_ins[1:])
    else:
        fmarray_0 = fmarray_ins[0]
        fmarray_1 = fmarray_ins[1]

    array = ma.concatenate(map(tomaskedarray, [fmarray_0, fmarray_1]), 0)
    table = np.concatenate([fmarray_0.table, fmarray_1.table])
    info = fmarray_0.info
    info.update(fmarray_1.info)

    fmarray_out = fm.FMArray(array, table, info)
    return fmarray_out


def mad(fmarray, axis=None, keepdims=False):
    """Compute the median absolute deviation (MAD) along the given axis.

    Args:
        fmarray (FMArray): An input fmarray.
        axis (int, optional): Axis along which the MADs are computed.
            The default is to compute the MAD along a flattened version of the array.
        keepdims (bool, optional): If True, the axes which are reduced are left
            in the result as dimensions with size one.

    Returns:
        mad (FMArray): A new array holding the result.

    """
    ad = ma.abs(fmarray - ma.median(fmarray, axis, keepdims=True))
    mad = ma.median(ad, axis, keepdims=keepdims)
    return mad


def median(fmarray, axis=None, out=None, overwrite_input=False, keepdims=False):
    """Compute the median along the spacified axis.

    It is equivalent to the np.ma.median function.
    i.e. fm.median(x, **kwargs) <=> np.ma.median(x, **kwargs)

    Args:
        fmarray (FMArray): An input fmarray.
        axis (int, optional):  Axis or axes along which the medians are computed.
        out (array, optional): Alternative output array in which to place the result.
        overwrite_input (bool, optional): If True, then allow use of memory of input
            fmarray for calculations. The input array will be modified by the call to median.
        keepdims (bool, optional): If True, the axes which are reduced are left
            in the result as dimensions with size one.

    Returns:
        median (FMArray): A new array holding the result.

    """
    median = ma.median(fmarray, axis, out, overwrite_input, keepdims)
    return median


def save(fmarray, filename=None):
    """Save a fmarray into a single file in uncompressed npz format.

    Args:
        fmarray (FMArray): An fmarray to be saved.
        filename (str): A file name (used as <filename>.npz).
            If not spacified, random 8-character name will be used.

    """
    data = fmarray.data
    mask = fmarray.mask
    table = fmarray.table
    info  = fmarray.info

    if filename is None:
        filename = uuid.uuid4().hex[:8]

    np.savez(filename, data=data, mask=mask, table=table, info=info)


def load(filename):
    """Load a fmarray from a npz file.

    Args:
        filename (str): A file name (*.npz).

    Returns:
        fmarray (FMArray): An output fmarray.

    """
    d = np.load(filename)
    array = ma.MaskedArray(d['data'], d['mask'])
    table = d['table']
    info  = d['info'].item()

    fmarray = fm.FMArray(array, table, info)
    return fmarray


def getfrequency(array, unit='GHz', **kwargs):
    if array.ismodulated:
        raise fm.utils.FMFlowError('array should be demodulated')

    info = array.info.copy()
    info.update(kwargs)
    rest = info['restfreq']
    step = info['chanwidth']
    fmindex = info['fmindex']

    start = rest - step * (0.5*(np.diff(fmindex)[0]-1)+fmindex[0])
    end   = start + step * array.shape[1]

    freq = np.arange(start, end, step) * u.Hz
    freq = freq.to(getattr(u, unit)).value
    return freq


def getspectrum(array, unit='K', weights=None):
    if array.ismodulated:
        raise fm.utils.FMFlowError('array should be demodulated')

    spec = ma.average(array, 0, weights) * u.K
    spec = spec.to(getattr(u, unit)).value
    return spec


# aliases
getfreq = getfrequency
getspec = getspectrum
