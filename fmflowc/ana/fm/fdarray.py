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


class FDDataArray(np.ndarray):
    def __new__(cls, array, fmch=None, fmid=None, cutch=None, info=None):
        obj = np.asarray(array).view(cls)

        if obj.ndim <= 1:
            obj.fmch = cls._init_param(fmch, np.zeros(1, int))
            obj.fmid = cls._init_param(fmid, (0, len(obj)))
        elif obj.ndim == 2:
            obj.fmch = cls._init_param(fmch, np.zeros(obj.shape[0], int))
            obj.fmid = cls._init_param(fmid, (0, obj.shape[1]))
        else:
            raise e.FMflowError('dimension of array must be <= 2')

        obj.cutch = cls._init_param(cutch, (0, 0))
        obj.info  = cls._init_param(info, {})
        return obj

    @staticmethod
    def _init_param(value, default):
        return default if value is None else deepcopy(value) 

    def __array_finalize__(self, obj):
        if obj is not None:
            self.fmch  = getattr(obj, 'fmch', None)
            self.fmid  = getattr(obj, 'fmid', None)
            self.cutch = getattr(obj, 'cutch', None)
            self.info  = getattr(obj, 'info', None)

    def __repr__(self):
        repr_string = 'FDDataArray({})'.format(self.__str__())
        repr_string = repr_string.replace('\n', '\n' + ' '*12)
        return repr_string


class FDArray(MaskedArray, FDDataArray):
    def __new__(cls, array, fmch=None, fmid=None, cutch=None, info=None):
        array = ma.asarray(array)
        data = FDDataArray(array.data, fmch, fmid, cutch, info)
        mask = array.mask

        obj = MaskedArray.__new__(cls, data, mask)
        obj.fmch  = data.fmch
        obj.fmid  = data.fmid
        obj.cutch = data.cutch
        obj.info  = data.info
        return obj

    def modulate(self):
        from fmflowc.ana import fm

        shape = (self.shape[0], self.shape[1]+np.sum(self.cutch))
        array = ma.zeros(shape, dtype=float)
        array.mask = True
        lcutch, rcutch = self.cutch
        array[:,lcutch:array.shape[1]-rcutch] = self.tomaskedarray()

        for i in range(len(array)):
            array[i] = np.roll(array[i], -self.fmid[0]-self.fmch[i])

        array = array[:,:self.fmid[1]]
        return fm.FMArray(array, self.fmch, self.info)

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

    def _parse_index_col(self, index):
        if type(index) != tuple:
            return (0, 0)

        index = index[1]
        if type(index) == int:
            return (index, self.shape[1]-(index+1))
        elif type(index) == slice:
            start, stop, stride = index.indices(self.shape[1])
            return (start, self.shape[1]-stop)
        elif type(index) == np.ndarray:
            cutch = lambda index: np.sum(np.cumsum(index)==0)
            return (cutch(index), cutch(index[::-1]))
        else:
            return (0, 0)

    def __array_finalize__(self, obj):
        if obj is not None:
            MaskedArray.__array_finalize__(self, obj)
            FDDataArray.__array_finalize__(self, obj)

    def __repr__(self):
        repr_string = 'FDArray({})'.format(self.__str__())
        repr_string = repr_string.replace('\n', '\n' + ' '*8)
        return repr_string

    def __getitem__(self, index):
        array  = self.tomaskedarray()[index]
        fmch   = self.fmch[self._parse_index_row(index)]
        dl, dr = self._parse_index_col(index)
        cutch  = (self.cutch[0]+dl, self.cutch[1]+dr)
        return FDArray(array, fmch, self.fmid, cutch, self.info)

    def __getslice__(self, i, j):
        return self.__getitem__(slice(i, j))
