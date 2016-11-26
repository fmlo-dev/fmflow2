# coding: utf-8

from __future__ import absolute_import as _absolute_import
from __future__ import division as _division
from __future__ import print_function as _print_function

# imported items
__all__ = ['gettsys', 'calibrate']

# constants
T_AMB = 293.0 # ambient temperature


def gettsys(r, sky, T_AMB=293.0):
    tsys = T_AMB / (r/sky-1)
    return tsys


def calibrate(on, off, r, sky, T_AMB=293.0):
    array = T_AMB * (on-off) / (r-sky)
    return array
