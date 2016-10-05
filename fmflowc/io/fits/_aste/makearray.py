# coding: utf-8

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

# dependent libraries
import numpy as np
from astropy.io import fits

# sub modules/functions
from fmflowc.ana import fm


def astearray(f, array_id, scan_type):
    '''
    
    '''
    if type(f) == str:
        f = fits.open(f)
