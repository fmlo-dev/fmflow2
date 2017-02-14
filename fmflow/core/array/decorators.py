# coding: utf-8

"""Module for decorators of fmarray functions.

"""

# Python 3.x compatibility
from __future__ import absolute_import as _absolute_import
from __future__ import division as _division
from __future__ import print_function as _print_function

# the Python standard library
from functools import partial, wraps
from inspect import getargspec

# the Python Package Index
import numpy as np
import numpy.ma as ma
import fmflow as fm

# imported items
__all__ = ['fmfunc', 'timechunk']


def fmfunc(func):
    """Make a function compatible with fmarray.

    This function is used as a decorator like::

        >>> @fm.fmfunc
        >>> def func(fmarray):
        ...     return fmarray # do nothing
        >>>
        >>> result = func(fmarray)

    Args:
        func (function): A function to be wrapped.
            The first argument of the function must be `fmarray`.

    Returns:
        wrapper (function): A wrapped function.

    """
    @wraps(func)
    def wrapper(*args, **kwargs):
        fmarray_in  = kwargs.pop('fmarray', args[0])
        fmarray_out = fmarray_in.copy()

        if type(fmarray_in) == fm.FMArray:
            array_in = fmarray_in.toarray()
        else:
            array_in = np.asarray(fmarray_in)

        array_out = func(array_in, *args[1:], **kwargs)
        fmarray_out[:] = array_out
        return fmarray_out

    return wrapper


def timechunk(func):
    """Make a function compatible with multicore time-chunk processing.

    This function is used as a decorator like::

        >>> @fm.fmfunc
        >>> @fm.timechunk
        >>> def func(fmarray):
        ...     return fmarray # do nothing
        >>>
        >>> result = func(fmarray, chunklen=100)

    Args:
        func (function): A function to be wrapped.
            The first argument of the function must be fmarray.

    Returns:
        wrapper (function): A wrapped function.

    """
    import fmflow as fm
    
    @wraps(func)
    def wrapper(*args, **kwargs):
        argnames = getargspec(func).args
        for i in range(len(args)):
            kwargs[argnames[i]] = args[i]

        array_in  = kwargs.pop('fmarray')
        chunk_len = kwargs.pop('chunklen', len(array_in))
        chunk_num = round(len(array_in)/chunk_len)

        index = np.linspace(0, len(array_in), chunk_num+1, dtype=int)
        subs_in = [array_in[index[i]:index[i+1]] for i in range(len(index)-1)]

        p = fm.utils.Pool()
        pfunc = partial(func, **kwargs)
        subs_out = p.map(pfunc, subs_in)
        array_out = np.concatenate(subs_out)
        return array_out

    return wrapper
