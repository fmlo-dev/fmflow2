# coding: utf-8

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

# standard libraries
from datetime import datetime, timedelta


ISO_8601 = '%Y-%m-%dT%H:%M:%S.%f'

def from_isoformat(string):
	return datetime.strptime(string, ISO_8601)

def to_isoformat(string, fmt):
	return datetime.strptime(string, fmt).strftime(ISO_8601)
