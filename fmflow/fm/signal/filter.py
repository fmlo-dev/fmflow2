# coding: utf-8

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

# dependent libraries
import numpy as np
from scipy import signal
from ..array.arrayfunc import *


@fmfunc
@timechunk
def hmedfilt(array_in, kernel=3):
        array_out = signal.medfilt(array_in, (1, kernel))
            return array_out
