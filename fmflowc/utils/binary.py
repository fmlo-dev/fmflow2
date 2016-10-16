# coding: utf-8

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

# standard libraries
import re
from struct import calcsize, unpack
from collections import OrderedDict

# dependent libraries
import numpy as np


def dumpstruct(f, defobj):
    '''Sequentially read and dump a data structure on a file object.

    Args:
    - f (file): a file object of the logging.
    - defobj (OrderedDict): ordered dict that defines the data structure.

    Returns:
    - struct (OrderedDict): ordered dict in which the data structure is stored.
    '''
    struct = OrderedDict()
    for key in defobj:
        struct[key] = dumpdata(f, *defobj[key])

    return struct


def dumpdata(f, fmt, shape=1):
    '''Sequentially read and unpack binary data on a file object.

    Args:
    - f (file): a file object of the logging.
    - fmt (str): format of the data that is compatible with the Python's
        struct module. for example: '>i' (int with big endian), '10s' (10 chars).
        for more information, see http://docs.python.jp/2/library/struct.html.
    - shape (int or tuple of int): data shape (if the data type is array).

    Returns:
    - data (object): the data on the file object.
    '''
    try:
        bytesize = calcsize(fmt)
    except:
        raise ValueError, fmt

    try:
        length = np.prod(shape)
    except:
        raise ValueError, shape

    data = []
    for i in range(length):
        bindata = f.read(bytesize)
        if re.search('s', fmt):
            data.append(unpack(fmt, bindata)[0].split('\x00')[0])
        else:
            data.append(unpack(fmt, bindata)[0])

    data = np.reshape(data, shape).tolist()
    data = data[0] if length == 1 else data

    return data


def fitsformat(defobj):
    '''Return the FITS format corresponding the date structure.

    Args:
    - defobj (OrderedDict): ordered dict that defines the data structure.

    Returns:
    - struct (OrderedDict): ordered dict in which the FITS format
        of the data structure is stored.
    '''
    struct = OrderedDict()
    for key in defobj:
        struct[key] = parseformat(*defobj[key])

    return struct


def parseformat(fmt, shape=1):
    '''Return a FITS format code from fmt and shape.

    Args:
    - fmt (str): format of the data that is compatible with the Python's
        struct module. for example: '>i' (int with big endian), '10s' (10 chars).
        for more information, see http://docs.python.jp/2/library/struct.html.
    - shape (int or tuple of int): data shape (if the data type is array).

    Returns:
    - fitsfmt (str): format of the data that is compatible with FITS.
    '''
    if re.search('s', fmt):
        code = 'A' + re.findall('\d+', fmt)[0]
    elif re.search('i', fmt):
        code = 'J'
    elif re.search('d', fmt):
        code = 'D'
    else:
        raise ValueError, fmt

    try:
        length = np.prod(shape)
    except:
        raise ValueError, shape

    fitsfmt = code if length == 1 else str(length) + code

    return fitsfmt
