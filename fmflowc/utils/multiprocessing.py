# coding: utf-8

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

# standard libraries
import multiprocessing as mp

# dependent libraries
import numpy as np


class Pool(object):
    def __init__(self, processes=None):
        self.processes = processes or mp.cpu_count()-1

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

    def map(self, func, iterable):
        if self.mpcompatible:
            return self.mpmap(func, iterable)
        else:
            return map(func, iterable)

    def mpmap(self, func, iterable):
        def pipefunc(conn, arg):
            conn.send(map(func, arg))
            conn.close()

        index = np.linspace(0, len(iterable), self.processes+1, dtype=int)
        args = [iterable[index[i]:index[i+1]] for i in range(len(index)-1)]

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

