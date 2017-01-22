# coding: utf-8

"""Module for multiprocessing in FMFlow.

This module provides helper classes for multiprocessing in FMFlow.

"""

# Python 3.x compatibility
from __future__ import absolute_import as __absolute_import
from __future__ import division as __division
from __future__ import print_function as __print_function

# the Python standard library
import multiprocessing as mp

# the Python Package Index
import numpy as np

# importing items
__all__ = ['Pool']


class Pool(object):
    """Return a process pool object.

    Attributes:
        mpcompatible (bool): Whether your NumPy/SciPy is compatible with multiprocessing.

    """

    def __init__(self, processes=None):
        """Initialize a process pool object.

        At this moment, whether your NumPy/SciPy is compatible with multiprocessing,
        is automatically checked and the result is stored in self.mpcompatible.

        Args:
            processes (int): The number of processes to be created. Default is
                <CPU count of your machine> -1 (one thread is saved for backup).

        """
        self.processes = processes or mp.cpu_count() - 1

        # check multiprocessing compatibility
        lapack_opt_info = np.__config__.lapack_opt_info

        if 'libraries' in lapack_opt_info:
            libraries = lapack_opt_info['libraries']

            if any([('mkl' in lib) for lib in libraries]):
                self.mpcompatible = True
            elif any([('blas' in lib) for lib in libraries]):
                self.mpcompatible = True
            elif any([('atlas' in lib) for lib in libraries]):
                self.mpcompatible = True
            else:
                self.mpcompatible = False
        else:
            self.mpcompatible = False

    def map(self, func, sequence):
        """Return a list of the results of applying the function to the sequence.

        If self.mpcompatible is True, mapping is multiprocessed with the spacified
        number of processes (default is <CPU count> - 1). If False, mapping is
        singleprocessed (equivalent to the bulitin map function).

        Args:
            func (function): Applying function.
            sequence (list): List of items to which function is applied.

        Returns:
            ret (list): The results of applying the function to the sequence.

        """
        if self.mpcompatible:
            return self._mpmap(func, sequence)
        else:
            return map(func, sequence)

    def _mpmap(self, func, sequence):
        def pipefunc(conn, arg):
            conn.send(map(func, arg))
            conn.close()

        index = np.linspace(0, len(sequence), self.processes+1, dtype=int)
        args = [sequence[index[i]:index[i+1]] for i in range(len(index)-1)]

        ret = []
        plist, clist = [], []
        for arg in args:
            pconn, cconn = mp.Pipe()
            plist.append(mp.Process(target=pipefunc, args=(cconn, arg)))
            clist.append(pconn)

        for p in plist:
            p.start()

        for c in clist:
            ret += c.recv()

        for p in plist:
            p.join()

        return ret
