# coding: utf-8

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

# dependent libraries
import numpy as np
from ..array.arrayfunc import *


@fmfunc
@timechunk
def pca(array_in, n_pc=0.9):
    if not (type(n_pc) == float or type(n_pc) == int):
        raise ValueError, 'n_pc must be float or int'

    # PCA using SVD
    U, d, Vt = np.linalg.svd(array_in, full_matrices=False)
    sumvar = np.cumsum(d**2)
    sumvar /= sumvar[-1]

    # reconstruction
    if type(n_pc) == float:
        n_pc = np.searchsorted(sumvar, n_pc) + 1

    array_out = np.dot(U[:,:n_pc], np.dot(np.diag(d[:n_pc]), Vt[:n_pc]))
    return array_out

