# coding: utf-8

"""Module for reading binary data in FMFlow.

Structure is an ordered dictionary (collections.OrderedDict),
which defines each name of parameter and corresponding (format, shape)
and will be used with the fmflow.utils.readbinary function.
Format must be compatible with the Python's struct module.
For example, '>i' (int with big endian), or '10s' (10 chars).
For more information, see (http://docs.python.jp/2/library/struct.html).

"""

# Python 3.x compatibility
from __future__ import absolute_import as _absolute_import
from __future__ import division as _division
from __future__ import print_function as _print_function

# the Python standard library
import re
from collections import OrderedDict
from struct import calcsize, unpack

# the Python Package Index
import numpy as np

# imported items
__all__ = ['readbinary', 'getfitsformat']


def readbinary(f, structure):
    """Sequentially read a data structure on a binary file object.

    Args:
        f (file): A binary file object of the logging.
        structure (OrderedDict): An ordered dictionary defining the data structure.

    Returns:
        readdata (OrderedDict): An ordered dictionary which stores the read data.

    """
    readdata = OrderedDict()
    for key in structure:
        readdata[key] = _readstructure(f, *structure[key])

    return readdata


def getfitsformat(structure):
    """Convert a data structure to the corresponding FITS format.

    Args:
        structure (OrderedDict): An ordered dictionary defining the data structure.

    Returns:
        formats (OrderedDict): An ordered dictionary which stores the FITS formats.

    """
    formats = OrderedDict()
    for key in structure:
        formats[key] = _convertformat(*structure[key])

    return formats


def _readstructure(f, structfmt, shape=1):
    try:
        bytesize = calcsize(structfmt)
    except:
        raise ValueError, structfmt

    try:
        length = np.prod(shape)
    except:
        raise ValueError, shape

    readdata = []
    for i in range(length):
        bindata = f.read(bytesize)
        if re.search('s', structfmt):
            readdata.append(unpack(structfmt, bindata)[0].split('\x00')[0])
        else:
            readdata.append(unpack(structfmt, bindata)[0])

    readdata = np.reshape(readdata, shape).tolist()
    readdata = readdata[0] if length == 1 else readdata

    return readdata


def _convertformat(structfmt, shape=1):
    if re.search('s', structfmt):
        code = 'A' + re.findall('\d+', structfmt)[0]
    elif re.search('i', structfmt):
        code = 'J'
    elif re.search('d', structfmt):
        code = 'D'
    else:
        raise ValueError, structfmt

    try:
        length = np.prod(shape)
    except:
        raise ValueError, shape

    fitsfmt = code if length == 1 else str(length) + code

    return fitsfmt
