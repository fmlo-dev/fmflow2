# coding: utf-8

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

# standard libraries
from copy import deepcopy

# dependent libraries¬
import numpy as np
import numpy.ma as ma
from fmflowc.utils import exceptions as e


class FMData(np.ndarray):
    def __new__(cls, array, table=None, info=None):
        '''
        
        '''
        # initialize array
        obj = np.asarray(array).copy()
        if obj.ndim == 0:
            return obj.item()
        elif obj.ndim == 1:
            return np.asarray(obj)

        obj = obj.view(cls)
        
        # initialize table
        _table = np.zeros(len(array), [('fmch','i8'), ('coord','2f8')])
        _table = table.copy() if table is not None else _table
        _table = _table.view(np.recarray)

        if not set(_table.dtype.names) <= {'fmch', 'coord'}:
            raise e.FMflowError('table must contain both fmch and coord column')

        # initialize info
        _info = {'modulated': True, 'cutch': (0,0), 'fmid': (0,0)}
        if info is not None:
            _info.update(deepcopy(info))

        # add attributes
        obj.table = _table
        obj.info  = _info
        return obj

    def _row(self, index):
        if type(index) == tuple:
            return self._row(index[0])
        elif type(index) == int:
            return slice(index, index+1, 1)
        else:
            return index

    def _cutch(self, index):
        if type(index) != tuple:
            return (0, 0)

        index = index[1]
        if type(index) == int:
            return (index, self.shape[1]-(index+1))
        elif type(index) == slice:
            start, stop, stride = index.indices(self.shape[1])
            return (start, self.shape[1]-stop)
        elif type(index) == np.ndarray:
            f = lambda index: np.sum(np.cumsum(index)==0)
            return (f(index), f(index[::-1]))
        else:
            return (0, 0)

    def __array_finalize__(self, obj):
        if obj is not None:
            self.table = getattr(obj, 'table', None)
            self.info  = getattr(obj, 'info', None)

    def __array_wrap__(self, obj):
        return np.ndarray.__array_wrap__(self, obj)

    def __getitem__(self, index):
        array = np.asarray(self)[index]
        table = self.table[self._row(index)]
        info  = deepcopy(self.info)

        cutch  = np.array(info['cutch'])
        dcutch = np.array(self._cutch(index))
        info['cutch'] = tuple(cutch+dcutch)
        return FMData(array, table, info)

    def __getslice__(self, i, j):
        return self.__getitem__(slice(i, j))

    def __repr__(self):
        string = 'fmdata({})'.format(self.__str__())
        return string.replace('\n', '\n'+' '*7)


class FMArray(ma.MaskedArray):
    def __new__(cls, array, fmch=None, coord=None, info=None):
        '''
        
        '''
        array = ma.asarray(array)
        fmdata, mask = FMData(array.data), array.mask
        obj = ma.MaskedArray(fmdata, mask).copy()

        if obj.ndim == 0:
            return obj.item()
        elif obj.ndim == 1:
            return np.asarray(obj)

        obj = obj.view(cls)

        if fmch is not None:
            obj.table.fmch = fmch.copy()

        if coord is not None:
            obj.table.coord = coord.copy()

        if info is not None:
            obj.update(deepcopy(info))

        return obj

    @property
    def fmch(self):
        return self.table.fmch

    @fmch.setter
    def fmch(self, value):
        self.table.fmch = value

    @property
    def coord(self):
        return self.table.coord

    @coord.setter
    def coord(self, value):
        self.table.coord = value

    @property
    def ismodulated(self):
        return self.info['modulated']

    @property
    def isdemodulated(self):
        return not self.info['modulated']

    @classmethod
    def fromfmdata(cls, fmdata, mask=None):
        '''
        
        '''
        obj = ma.MaskedArray(fmdata, mask).copy()

        if obj.ndim == 0:
            return obj.item()
        elif obj.ndim == 1:
            return np.asarray(obj)
        else:
            return obj.view(cls)

    @classmethod
    def frommaskedarray(cls, array, table, info):
        '''
        
        '''
        array = ma.asarray(array)
        fmdata, mask = FMData(array.data, table, info), array.mask
        obj = ma.MaskedArray(fmdata, mask).copy()

        if obj.ndim == 0:
            return obj.item()
        elif obj.ndim == 1:
            return np.asarray(obj)
        else:
            return obj.view(cls)

    def tofmdata(self):
        '''
        
        '''
        fmdata = self.data.copy()
        return fmdata

    def tomaskedarray(self):
        '''
        
        '''
        data, mask = np.asarray(self), self.mask
        array = ma.MaskedArray(data, mask).copy()
        return array

    def demodulate(self):
        '''
        
        '''
        if self.isdemodulated:
            raise e.FMflowError('this array is already demodulated')

        array = ma.zeros((self.shape[0], self.shape[1]+np.ptp(self.fmch)))
        array.mask = True
        array[:,:self.shape[1]] = self.tomaskedarray()

        for i in range(len(array)):
            array[i] = np.roll(array[i], -np.min(self.fmch)+self.fmch[i])

        table = self.table
        info  = deepcopy(self.info)
        fmid  = (-np.min(self.fmch), self.shape[1])
        info.update({'modulated': False, 'cutch': (0,0), 'fmid': fmid})
        return FMArray.frommaskedarray(array, table, info)

    def modulate(self):
        '''
        
        '''
        if self.ismodulated:
            raise e.FMflowError('this array is already modulated')

        lcutch, rcutch = self.info['cutch']

        array = ma.zeros((self.shape[0], self.shape[1]+lcutch+rcutch))
        array.mask = True
        array[:,lcutch:array.shape[1]-rcutch] = self.tomaskedarray()

        for i in range(len(array)):
            array[i] = np.roll(array[i], -self.info['fmid'][0]-self.fmch[i])

        array = array[:,:self.info['fmid'][1]]
        
        table = self.table
        info  = deepcopy(self.info)
        info.update({'modulated': True, 'cutch': (0,0), 'fmid': (0,0)})
        return FMArray.frommaskedarray(array, table, info)

    def __array_finalize__(self, obj):
        if obj is not None:
            ma.MaskedArray.__array_finalize__(self, obj)

    def __array_wrap__(self, obj):
        return ma.MaskedArray.__array_wrap__(self, obj)

    def __getitem__(self, index):
        fmdata = self.data[index]
        mask = self.mask[index] if self.mask.ndim else self.mask
        return FMArray.fromfmdata(fmdata, mask)

    def __getslice__(self, i, j):
        return self.__getitem__(slice(i, j))

    def __repr__(self):
        string = 'fmarray({})'.format(self.__str__())
        return string.replace('\n', '\n'+' '*8)

