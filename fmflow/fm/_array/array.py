# coding: utf-8

"""Module for fmarray.

fmarray is one of the fundamental data formats in FMFlow,
which is created as an instance of FMArray class.

FMArray is subclass of NumPy masked array, which handles timestream array
together with timestream table (fmch, coord, etc) and info dictionary.

Normally, fmarray is created and operated with the functions in fmflow.fm
(e.g. array, (de)modulate, etc) and FMArray is thus not touched by users.

Note:
    FMArray is imported only in the fmflow.fm.arrayfunc module.

"""

# Python 3.x compatibility
from __future__ import absolute_import as __absolute_import
from __future__ import division as __division
from __future__ import print_function as __print_function

# the Python Package Index
import numpy as np
import numpy.ma as ma
from fmflow import utils as ut

# imported items
__all__ = ['FMArray']

# constants
TABLE = lambda shape: np.zeros(shape, [('fmch', 'i8'), ('coord', '2f8')])
INFO  = lambda: {'fmstatus': 'FM', 'fmindex': (0,0), 'fmcutcol': (0,0)}


class FMArray(ma.MaskedArray):
    """A class for fmarray in FMFlow.

    It is subclass of NumPy masked array, which handles timestream array
    together with timestream table (fmch, coord, etc) and info dictionary.

    Normally, fmarray is created and operated with the functions in fmflow.fm
    (e.g. array, (de)modulate, etc) and FMArray is thus not touched by users.

    Attributes:
        table: A timestream table that stores fmch, coord, etc.
        fmch: An array of modulation frequencies in units of channel.
        coord: An array of observed coordinates in units of degrees.
        info: A dictionary that stores information of the fmarray.
        isdemodulate: A boolean that indicates whether the fmarray is demodulated.
        ismodulate: A boolean that indicates whether the fmarray is modulated.

    """

    def __new__(cls, array, table=None, info=None):
        """Create an modulated fmarray from (masked) array, table, and info.

        Normally, this method is for the internal use.
        Users should create an fmarray with the fm.array function.

        Args:
            array (masked array): A timestream array (mask is optional).
            table (record array): A timestream table that stores fmch, coord, etc.
            info (dict): A dictionary that stores information of the fmarray.

        Returns:
            obj (FMArray): A modulated fmarray.

        """
        # array
        array = ma.asarray(array).copy()
        data, mask = array.data, array.mask
        obj = super(FMArray, cls).__new__(cls, data, mask)

        if obj.ndim == 0:
            return obj.item()

        # table
        _shape = 1 if obj.ndim<=1 else len(obj)
        _table = TABLE(_shape) if table is None else table.copy()
        _table = _table.view(np.recarray)

        # info
        _info = INFO()
        _info.update({} if info is None else info.copy())

        # finally
        obj._optinfo['table'] = _table
        obj._optinfo['info']  = _info
        return obj

    @classmethod
    def fromeach(cls, array=None, fmch=None, coord=None, info=None):
        """Create a modulated fmarray from each argument.

        Normally, this method is for the internal use.
        It is equivalent to the fm.array function (recommended to use).
        i.e. FMArray.fromeach(...) <=> fm.array(...)

        Args:
            array (masked array): A timestream array (mask is optional).
            fmch (array): Array of modulation frequencies in units of channel.
            coord (array): Array of observed coordinates in units of degrees.
            info (dict): Information of the fmarray, observation, etc.

        Returns:
            obj (FMArray): A modulated fmarray.

        """
        obj = cls.__new__(cls, array, info=info)

        # update table of array
        if fmch is not None:
            obj.fmch = fmch.copy()

        if coord is not None:
            obj.coord = coord.copy()

        return obj

    def demodulate(self, reverse=False):
        """Create a demodulated fmarray from the modulated one.

        This method is only available when the original fmarray is modulated.
        It is equivalent to the fm.demodulate function (recommended to use).
        i.e. x.demodulate(reverse) <=> fm.demodulate(x, reverse)

        Args:
            reverse (bool): If True, the original fmarray is reverse-demodulated
                (i.e. fmch * -1 is used for demodulation). Default is False.

        Returns:
            array (FMArray): A demodulated fmarray.

        """
        if self.isdemodulated:
            raise ut.FMFlowError('this fmarray is already demodulated')

        # table, info
        table = self.table.copy()
        info  = self.info.copy()

        if reverse:
            table.fmch *= -1
            fmstatus = '-FD'
        else:
            fmstatus = '+FD'

        fmindex0 = -np.min(table.fmch)
        fmindex  = (fmindex0, fmindex0+self.shape[1])

        # array
        array = ma.zeros((self.shape[0], self.shape[1]+np.ptp(table.fmch)))
        array.mask = True
        array[:,:self.shape[1]] = self.asmaskedarray()

        for i in range(len(array)):
            array[i] = np.roll(array[i], fmindex[0]+table.fmch[i])

        # finally
        info.update({'fmstatus': fmstatus, 'fmindex': fmindex, 'fmcutcol': (0,0)})
        return FMArray(array, table, info)

    def modulate(self):
        """Create a modulated fmarray from the demodulated one.

        This method is only available when the original fmarray is demodulated.
        It is equivalent to the fm.modulate function (recommended to use).
        i.e. x.modulate() <=> fm.modulate(x)

        Returns:
            array (FMArray): A modulated fmarray

        """
        if self.ismodulated:
            raise ut.FMFlowError('this fmarray is already modulated')

        # table, info
        table = self.table.copy()
        info  = self.info.copy()

        fmstatus = info['fmstatus']
        fmindex  = info['fmindex']
        fmcutcol = info['fmcutcol']

        # array
        array = ma.zeros((self.shape[0], self.shape[1]+np.sum(fmcutcol)))
        array.mask = True
        array[:,fmcutcol[0]:array.shape[1]-fmcutcol[1]] = self.asmaskedarray()

        for i in range(len(array)):
            array[i] = np.roll(array[i], -(fmindex[0]+table.fmch[i]))

        array = array[:,:np.diff(fmindex)[0]]

        # finally
        if fmstatus == '-FD':
            table.fmch *= -1

        info.update({'fmstatus': 'FM', 'fmindex': (0,0), 'fmcutcol': (0,0)})
        return FMArray(array, table, info)

    def asndarray(self):
        """Convert the fmarray to a NumPy ndarray.
        
        It is equivalent to the fm.asndarray function (recommended to use).
        i.e. x.asndarray() <=> fm.asndarray(x)
        
        """
        array = self.data.copy()
        return array

    def asmaskedarray(self):
        """Convert the fmarray to a NumPy masked array.

        It is equivalent to the fm.asmaskedarray function (recommended to use).
        i.e. x.asmaskedarray() <=> fm.asmaskedarray(x)

        """
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
        return 'FM' in self._optinfo['info']['fmstatus']

    @property
    def isdemodulated(self):
        return 'FD' in self._optinfo['info']['fmstatus']

    def _tableindex(self, index):
        """Convert the index for array to one for table."""
        if type(index) in (tuple, list):
            return self._tableindex(index[0])
        elif type(index) == int:
            return slice(index, index+1, 1)
        else:
            return index

    def _fmcutcol(self, index):
        """Calculate channels of edge-cut column from an index."""
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
        """x.__getitem__(y) <=> x[y]"""
        # array
        array = self.asmaskedarray()[index]

        # table
        try:
            tableindex = self._tableindex(index)
            table = self.table.copy()[tableindex]
        except:
            table = None

        # info
        try:
            info = self.info.copy()
            fmcutcol = np.array(info['fmcutcol'])
            dfmcutcol = np.array(self._fmcutcol(index))
            info.update({'fmcutcol': tuple(fmcutcol+dfmcutcol)})
        except:
            info = None

        return FMArray(array, table, info)

    def __getslice__(self, i, j):
        """x.__getslice__(i, j) <=> x[i:j]"""
        return self.__getitem__(slice(i, j))

    def __repr__(self):
        """x.__repr__() <=> repr(x)"""
        string = 'fmarray({})'.format(self.__str__())
        return string.replace('\n', '\n'+' '*8)