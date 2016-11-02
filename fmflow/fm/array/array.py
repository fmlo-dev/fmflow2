# coding: utf-8

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

# dependent libraries¬
import numpy as np
import numpy.ma as ma
from fmflow.utils import exceptions as e


class FMArray(ma.MaskedArray):
    def __new__(cls, array, table=None, info=None):
        # masked array
        array = ma.asarray(array).copy()
        data, mask = array.data, array.mask
        obj = super(FMArray, cls).__new__(cls, data, mask)

        if obj.ndim == 0:
            return obj.item()

        # table
        _shape = 1 if obj.ndim<=1 else len(obj)
        _table = np.zeros(_shape, [('fmch', 'i8'), ('coord', '2f8')])
        _table = table.copy() if table is not None else _table
        _table = _table.view(np.recarray)

        if not set(_table.dtype.names) <= {'fmch', 'coord'}:
            raise e.FMflowError('table must contain fmch and coord columns')

        # info
        _info = {'modulated': True, 'cutch': (0,0), 'fmid': (0,0)}
        _info.update(info.copy() if info is not None else {})
        
        # add table and info to array
        obj._optinfo['table'] = _table
        obj._optinfo['info']  = _info
        return obj

    @classmethod
    def fromeach(cls, array=None, fmch=None, coord=None, info=None):
        obj = cls.__new__(cls, array, info=info)

        # update table of array
        if fmch is not None:
            obj.fmch = fmch.copy()

        if coord is not None:
            obj.coord = coord.copy()

        return obj

    def demodulate(self):
        '''Return a demodulated array.'''
        if self.isdemodulated:
            raise e.FMflowError('this array is already demodulated')

        fmid  = (-np.min(self.fmch), self.shape[1])

        # masked array
        array = ma.zeros((self.shape[0], self.shape[1]+np.ptp(self.fmch)))
        array.mask = True
        array[:,:self.shape[1]] = self.asmaskedarray()

        for i in range(len(array)):
            array[i] = np.roll(array[i], fmid[0]+self.fmch[i])

        # table
        table = self.table

        # info
        info  = self.info.copy()
        info.update({'modulated': False, 'cutch': (0,0), 'fmid': fmid})

        return FMArray(array, table, info)

    def modulate(self):
        '''Return a modulated array.'''
        if self.ismodulated:
            raise e.FMflowError('this array is already modulated')

        fmid = self.info['fmid']
        lcutch, rcutch = self.info['cutch']

        # masked array
        array = ma.zeros((self.shape[0], self.shape[1]+lcutch+rcutch))
        array.mask = True
        array[:,lcutch:array.shape[1]-rcutch] = self.asmaskedarray()

        for i in range(len(array)):
            array[i] = np.roll(array[i], -fmid[0]-self.fmch[i])

        array = array[:,:fmid[1]]

        # table
        table = self.table

        # info
        info  = self.info.copy()
        info.update({'modulated': True, 'cutch': (0,0), 'fmid': (0,0)})

        return FMArray(array, table, info)

    def asmaskedarray(self):
        '''Return a masked array (table and info are deleted).'''
        data, mask = np.asarray(self), self.mask
        array = ma.MaskedArray(data, mask).copy()
        return array

    @property
    def fmch(self):
        return self._optinfo['table'].fmch

    @property
    def coord(self):
        return self._optinfo['table'].coord

    @fmch.setter
    def fmch(self, value):
        self._optinfo['table'].fmch = value

    @coord.setter
    def coord(self, value):
        self._optinfo['table'].coord = value

    @property
    def table(self):
        return self._optinfo['table']

    @property
    def info(self):
        return self._optinfo['info']

    @property
    def ismodulated(self):
        return self._optinfo['info']['modulated']

    @property
    def isdemodulated(self):
        return not self._optinfo['info']['modulated']

    def _row(self, index):
        if type(index) in (tuple, list):
            return self._row(index[0])
        elif type(index) == int:
            return slice(index, index+1, 1)
        else:
            return index

    def _cutch(self, index):
        try:
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
        except:
            return (0, 0)

    def __getitem__(self, index):
        # masked array
        array = self.asmaskedarray()[index]

        # table
        try:
            table = self.table[self._row(index)]
        except:
            table = None

        # info
        try:
            info = self.info.copy()
            cutch = np.array(info['cutch'])
            dcutch = np.array(self._cutch(index))
            info['cutch'] = tuple(cutch+dcutch)
        except:
            info = None

        return FMArray(array, table, info)

    def __getslice__(self, i, j):
        return self.__getitem__(slice(i, j))

    def __repr__(self):
        string = 'fmarray({})'.format(self.__str__())
        return string.replace('\n', '\n'+' '*8)
