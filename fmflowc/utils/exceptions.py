# coding: utf-8

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function


class FMflowError(Exception):
    def __init__(self, message):
        self.message = message

    def __str__(self):
        return repr(self.message)


class FMflowWarning(Warning):
    def __init__(self, message):
        self.message = message

    def __str__(self):
        return repr(self.message)
