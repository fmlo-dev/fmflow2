# coding: utf-8

"""Module for date and time format in FMFlow.

This module provides converter functions of datetime <==> string
dealing with ISO 8601 format (YYYY-mm-ddTHH:MM:SS.ssssss).

Available classes:
- DatetimeConverter: Convert a datetime string in given format to one in ISO format.
"""

# Python 3.x compatibility
from __future__ import absolute_import as __absolute_import
from __future__ import division as __division
from __future__ import print_function as __print_function

# the Python standard library
from datetime import datetime, timedelta

# importing items
__all__ = ['DatetimeConverter']

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
