# coding: utf-8

# Python 3.x compatibility
from __future__ import absolute_import as _absolute_import
from __future__ import division as _division
from __future__ import print_function as _print_function

# the Python Standard Library
from collections import defaultdict
from copy import deepcopy

# the Python Package Index
import numpy as np
import fmflow as fm
from sklearn import decomposition

# imported items
__all__ = ['reducedim']

# constants
PARAMS = defaultdict(dict)
PARAMS['KernelPCA'] = {'fit_inverse_transform': True}


@fm.fmfunc
@fm.timechunk
def reducedim(fmarray, decomposer='TruncatedSVD', **kwargs):
    """Compute a dimension-reduced fmarray via a decomposition algorithm.

    Args:
        fmarray (FMArray): An input fmarray.
        decomposer (str): A name of decomposition class
            which sklearn.decomposition provides.
        kwargs (dict): Parameters for the spacified algorithm such as `n_components`.

    Returns:
        result (FMArray): An output dimension-reduced fmarray.

    Example:
        To compute a fmarray reconstructed from top two principal components:

        >>> result = fm.model.reduceim(fmarray, 'PCA', n_components=2)

    """
    AlgorithmClass = getattr(decomposition, decomposer)
    params = deepcopy(PARAMS[decomposer])
    params.update(kwargs)

    model = AlgorithmClass(**params)
    fit = model.fit_transform(fmarray)
    
    if hasattr(model, 'inverse_transform'):
        result = model.inverse_transform(fit)
    elif hasattr(model, 'components_'):
        result = np.dot(fit, model.components_)
    else:
        raise fm.utils.FMFlowError('cannot decompose with the spacified algorithm')

    return result
