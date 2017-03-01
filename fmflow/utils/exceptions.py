# coding: utf-8

# Python 3.x compatibility
from __future__ import absolute_import as _absolute_import
from __future__ import division as _division
from __future__ import print_function as _print_function

# importing items
__all__ = ['FMFlowError', 'FMFlowWarning']


class FMFlowError(Exception):
    """Error class of FMFlow.

    """

    def __init__(self, message):
        self.message = message

    def __str__(self):
        return repr(self.message)


class FMFlowWarning(Warning):
    """Warning class of FMFlow.

    """

    def __init__(self, message):
        self.message = message

    def __str__(self):
        return repr(self.message)
