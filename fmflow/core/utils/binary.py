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
import json
import re
from collections import OrderedDict
from struct import Struct

# the Python Package Index
import numpy as np

# imported items
__all__ = ['CStructReader']


class CStructReader(object):
    def __init__(self, fields, ignored=None, byteorder=None):
        self.fields = fields
        self.ignored = ignored or '$.'
        self.byteorder = byteorder or '@'
        self.formats, self.shapes = self._parsefields()
        self.struct = Struct(self.joinedformat)
        self._data = self._initdata()

    def read(self, f):
        bindata = f.read(self.struct.size)
        unpdata = list(self.struct.unpack(bindata))
        for key in self.shapes:
            shape = self.shapes[key]
            count = np.prod(shape)
            datum = [unpdata.pop(0) for i in range(count)]
            self._data[key].append(np.asarray(datum))

    @property
    def data(self):
        data = OrderedDict()
        for key in self.shapes:
            if re.search(self.ignored, key):
                continue

            shape = [len(self._data[key])] + self.shapes[key]
            datum = np.reshape(self._data[key], shape)
            if np.prod(shape) == 1:
                data[key] = np.squeeze(datum).item()
            else:
                data[key] = np.squeeze(datum)

        return data

    @property
    def jsondata(self):
        data = self.data
        for key in data:
            if isinstance(data[key], np.ndarray):
                data[key] = data[key].tolist()

        jsondata = json.dumps(data)
        return jsondata

    @property
    def fitsformats(self):
        fitsformats = OrderedDict()
        for key in self.formats:
            fmt = self.formats[key]
            shape = self.shapes[key]
            count = np.prod(shape)

            if re.search('s', fmt):
                code = 'A'
                code += re.findall('\d+', fmt)[0]
            elif re.search('i', fmt):
                code = 'J'
            elif re.search('d', fmt):
                code = 'D'
            else:
                raise ValueError, fmt

            if count == 1:
                fitsfmt = code
            else:
                fitsfmt = str(count) + code

            fitsformats[key] = fitsfmt

        return fitsformats

    @property
    def joinedformat(self):
        joinedformat = self.byteorder
        for key in self.formats:
            fmt = self.formats[key]
            count = np.prod(self.shapes[key])
            joinedformat += fmt * count

        return joinedformat

    def _parsefields(self):
        formats = OrderedDict()
        fitsformats = OrderedDict()
        shapes = OrderedDict()
        for field in self.fields:
            if len(field) == 2:
                key, fmt = field
                shape = [1]
            elif len(field) == 3:
                key, fmt, shape = field
                if type(shape) == int:
                    shape = [shape]

            formats[key] = fmt
            shapes[key]  = shape

        return formats, shapes

    def _initdata(self):
        _data = OrderedDict()
        for key in self.formats:
            _data[key] = []

        return _data
