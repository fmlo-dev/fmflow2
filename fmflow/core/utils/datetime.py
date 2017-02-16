# coding: utf-8

"""Module for date and time format in FMFlow.

This module provides parser class of datetime <==> string
dealing with ISO 8601 format (YYYY-mm-ddTHH:MM:SS.ssssss).

"""

# Python 3.x compatibility
from __future__ import absolute_import as _absolute_import
from __future__ import division as _division
from __future__ import print_function as _print_function

# the Python standard library
from datetime import datetime, timedelta

# importing items
__all__ = ['DatetimeParser']

# constants
ISO_8601 = '%Y-%m-%dT%H:%M:%S.%f'
PATTERNS = [
    '%Y%m%d%H%M%S',
    '%Y%m%d%H%M%S.',
    '%Y%m%d%H%M%S.%f',
    '%y%m%d%H%M%S',
    '%y%m%d%H%M%S.',
    '%y%m%d%H%M%S.%f',
    ISO_8601,
]


class DatetimeParser(object):
    """Convert a datetime string to one in ISO format.

    """

    def __init__(self, rounddown=True):
        self.rounddown = rounddown
        self.pattern = ISO_8601
    
    def __call__(self, string):
        try:
            dt = datetime.strptime(string, self.pattern)
        except:
            self.setpattern(string)
            dt = datetime.strptime(string, self.pattern)
        
        dt_iso = dt.strftime(ISO_8601)
        if self.rounddown:
            dt_iso = dt_iso[:-5] + '00000'
            
        return dt_iso

    def setpattern(self, string):
        for pattern in PATTERNS:
            try:
                dt = datetime.strptime(string, pattern)
                self.pattern = pattern
            except ValueError:
                pass

