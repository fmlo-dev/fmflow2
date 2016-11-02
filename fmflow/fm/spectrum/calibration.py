# coding: utf-8

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

# dependent libraries¬
import numpy as np
import numpy.ma as ma
import astropy.units as u
from fmflow.utils import exceptions as e


def gettsys(r, sky, T_amb=293.0):
    tsys = T_amb / (r/sky-1)
    return tsys

def calibrate(on, off, r, sky, T_amb=293.0):
    array = T_amb * (on-off) / (r-sky)
    return array
