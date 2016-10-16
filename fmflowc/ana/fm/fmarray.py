# coding: utf-8

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

# standard libraries
from copy import deepcopy

# dependent libraries¬
import numpy as np
import numpy.ma as ma
from numpy.ma import MaskedArray
from fmflowc.utils import exceptions as e


class FMDataArray(np.ndarray):
    def __new__(cls, array, fmch=None, info=None):
        obj = np.asarray(array).view(cls)

        if obj.ndim <= 1:
            obj.fmch = cls._init_param(fmch, np.zeros(1, int))
        elif obj.ndim == 2:
            obj.fmch = cls._init_param(fmch, np.zeros(obj.shape[0], int))
        else:
            raise e.FMflowError('dimension of array must be <= 2')

        obj.info = cls._init_param(info, {})
        return obj

    @staticmethod
    def _init_param(value, default):
        return default if value is None else deepcopy(value) 

    def __array_finalize__(self, obj):
        if obj is not None:
            self.fmch = getattr(obj, 'fmch', None)
            self.info = getattr(obj, 'info', None)

    def __repr__(self):
        repr_string = 'FMDataArray({})'.format(self.__str__())
        repr_string = repr_string.replace('\n', '\n' + ' '*12)
        return repr_string


class FMArray(MaskedArray, FMDataArray):
    def __new__(cls, array, fmch=None, info=None):
        array = ma.asarray(array)
        data = FMDataArray(array.data, fmch, info)
        mask = array.mask

        obj = MaskedArray.__new__(cls, data, mask)
        obj.fmch = data.fmch
        obj.info = data.info
        return obj

    def demodulate(self):
        from fmflowc.ana import fm

        fmid  = (-np.min(self.fmch), self.shape[1])
        shape = (self.shape[0], self.shape[1]+np.ptp(self.fmch))
        array = ma.zeros(shape, dtype=float)
        array.mask = True
        array[:,:self.shape[1]] = self.tomaskedarray()

        for i in range(len(array)):
            array[i] = np.roll(array[i], fmid[0]+self.fmch[i])

        return fm.FDArray(array, self.fmch, fmid, info=self.info)

    def tomaskedarray(self):
        data = deepcopy(np.asarray(self.data))
        mask = deepcopy(np.asarray(self.mask))
        return MaskedArray(data, mask)

    def _parse_index_row(self, index):
        if type(index) == tuple:
            return self._parse_index_row(index[0])
        elif type(index) == int:
            return slice(index, index+1, 1)
        else:
            return index

    def __array_finalize__(self, obj):
        if obj is not None:
            MaskedArray.__array_finalize__(self, obj)
            FMDataArray.__array_finalize__(self, obj)

    def __repr__(self):
        repr_string = 'FMArray({})'.format(self.__str__())
        repr_string = repr_string.replace('\n', '\n' + ' '*8)
        return repr_string

    def __getitem__(self, index):
        array = self.tomaskedarray()[index]
        fmch  = self.fmch[self._parse_index_row(index)]
        return FMArray(array, fmch, self.info)

    def __getslice__(self, i, j):
        return self.__getitem__(slice(i, j))
