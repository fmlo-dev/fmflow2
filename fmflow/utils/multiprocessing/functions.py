# coding: utf-8

# Python 3.x compatibility
from __future__ import absolute_import as _absolute_import
from __future__ import division as _division
from __future__ import print_function as _print_function

# the Python standard library
from functools import partial, wraps
from inspect import getargspec

# the Python Package Index
import fmflow as fm
import numpy as np

# imported items
__all__ = ['numchunk', 'timechunk']


def numchunk(func):
    """Make a function compatible with multicore numchunk processing.

    This function is used as a decorator like::

        >>> @fm.utils.numchunk
        >>> def func(array):
        ...     return array # do nothing
        >>>
        >>> result = func(array, numchunk=10)

    Args:
        func (function): A function to be wrapped. The first argument
            of the function must be an array to be num-chunked.

    Returns:
        wrapper (function): A wrapped function.

    """
    @wraps(func)
    def wrapper(*args, **kwargs):
        array_in = args[0]

        p = fm.utils.MPPool()
        N = kwargs.pop('numchunk', p.processes)

        argnames = getargspec(func).args
        if len(args) > 1:
            for i in range(1, len(args)):
                kwargs[argnames[i]] = args[i]

        indxs = np.linspace(0, len(array_in), N+1, dtype=int)
        subarrays_in = [array_in[indxs[i]:indxs[i+1]] for i in range(N)]

        pfunc = partial(func, **kwargs)
        subarrays_out = p.map(pfunc, subarrays_in)
        array_out = np.concatenate(subarrays_out)
        return array_out

    return wrapper


def timechunk(func):
    """Make a function compatible with multicore timechunk processing.

    This function is used as a decorator like::

        >>> @fm.utils.timechunk
        >>> def func(array):
        ...     return array # do nothing
        >>>
        >>> result = func(array, timechunk=100)

    Args:
        func (function): A function to be wrapped. The first argument
            of the function must be an array to be time-chunked.

    Returns:
        wrapper (function): A wrapped function.

    """
    @wraps(func)
    def wrapper(*args, **kwargs):
        array_in = args[0]

        p = fm.utils.MPPool()
        T = kwargs.pop('timechunk', len(array_in))
        N = int(round(len(array_in) / T))

        argnames = getargspec(func).args
        if len(args) > 1:
            for i in range(1, len(args)):
                kwargs[argnames[i]] = args[i]

        indxs = np.linspace(0, len(array_in), N+1, dtype=int)
        subarrays_in = [array_in[indxs[i]:indxs[i+1]] for i in range(N)]

        pfunc = partial(func, **kwargs)
        subarrays_out = p.map(pfunc, subarrays_in)
        array_out = np.concatenate(subarrays_out)
        return array_out

    return wrapper

