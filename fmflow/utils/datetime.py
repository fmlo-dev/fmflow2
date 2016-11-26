# coding: utf-8

"""Module for date and time format in FMFlow.

This module provides converter functions of datetime <==> string
dealing with ISO 8601 format (YYYY-mm-ddTHH:MM:SS.ssssss).

Available classes:
- DatetimeConverter: Convert a datetime string in given format to one in ISO format.
"""

from __future__ import absolute_import as _absolute_import
from __future__ import division as _division
from __future__ import print_function as _print_function

# importing items
__all__ = ['DatetimeConverter']

# the standard library
from datetime import datetime, timedelta

# constants
ISO_8601 = '%Y-%m-%dT%H:%M:%S.%f'


class DatetimeConverter(object):
    """Convert a datetime string in given format to one in ISO format.
    
    Public methods:
    - __call__: c.__call__(string) is equivalent to c(string).
    """
    def __init__(self, fmt):
        self.fmt = fmt
    
    def __call__(self, string):
        """c.__call__(string) is equivalent to c(string)."""
        return datetime.strptime(string, self.fmt).strftime(ISO_8601)
