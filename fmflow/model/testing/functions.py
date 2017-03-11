# coding: utf-8

# Python 3.x compatibility
from __future__ import absolute_import as _absolute_import
from __future__ import division as _division
from __future__ import print_function as _print_function

# the Python standard library
from functools import wraps
from warnings import warn

# the Python Package Index
import fmflow as fm
import numpy as np
from scipy.special import gammaln

# imported items
__all__ = ['pca', 'ppca']


def testing(func):
    @wraps(func)
    def wrapper(*args, **kwargs):
        warn(fm.utils.FMFlowWarning('this function is under testing'))
        return func(*args, **kwargs)
    
    return wrapper


@testing
@fm.arrayfunc
@fm.utils.timechunk
def pca(fmarray_in, n_pc=0.9):
    if not (type(n_pc) == float or type(n_pc) == int):
        raise ValueError, 'n_pc must be float or int'

    # PCA using SVD
    U, d, Vt = np.linalg.svd(fmarray_in, full_matrices=False)
    sumvar = np.cumsum(d**2)
    sumvar /= sumvar[-1]

    # reconstruction
    if type(n_pc) == float:
        n_pc = np.searchsorted(sumvar, n_pc) + 1

    fmarray_out = np.dot(U[:,:n_pc], np.dot(np.diag(d[:n_pc]), Vt[:n_pc]))
    return fmarray_out


@testing
@fm.arrayfunc
@fm.utils.timechunk
def ppca(fmarray_in, max_pc=None):
    # PCA using SVD
    U, d, Vt = np.linalg.svd(fmarray_in, full_matrices=False)

    # estimate n_pc
    N, D = fmarray_in.shape
    max_pc = N if max_pc is None else max_pc

    l_org = np.zeros(D)
    l_org[:len(d)] = d**2 / len(d)

    def laplace(k):
        assert 0 < k < N
        m = D*k - k*(k+1)/2.0
        v_opt = np.sum(l_org[k:]) / (D-k)
        l_opt = np.zeros_like(l_org)
        l_opt[:k] = l_org[:k]
        l_opt[k:] = v_opt

        def _Az():
            l_org_v = l_org[:k].reshape((k,1))
            l_org_h = l_org
            l_opt_v = l_opt[:k].reshape((k,1))
            l_opt_h = l_opt

            Az_map = np.log10((1/l_opt_h-1/l_opt_v) * (l_org_v-l_org_h) * N)
            Az = np.sum(np.ma.masked_invalid(np.triu(Az_map)))
            return Az

        def _pU():
            i  = np.arange(k)+1.0
            q  = (D-i+1.0)/2.0

            pU = 0.0
            pU += -k * np.log10(2)
            pU += np.sum(-q * np.log10(np.pi))
            pU += np.sum((gammaln(q)/np.log(10)))
            return pU

        prob_k = _pU()
        prob_k += (-N/2.0) * np.sum(np.log10(l_org[:k]))
        prob_k += (-N*(D-k)/2.0) * np.log10(v_opt)
        prob_k += ((m+k)/2.0) * np.log10(2.0*np.pi)
        prob_k += -0.5 * _Az()
        prob_k += -(k/2.0) * np.log10(N)
        return prob_k

    probs = map(laplace, np.arange(1, max_pc, dtype=int))
    n_pc = np.argmax(probs) + 1

    # reconstruction
    fmarray_out = np.dot(U[:,:n_pc], np.dot(np.diag(d[:n_pc]), Vt[:n_pc]))
    return fmarray_out
