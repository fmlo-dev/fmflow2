# coding: utf-8

"""Module for defining structures of OTF logging (common).

Structure is an ordered dictionary (collections.OrderedDict),
which defines each name of parameter and corresponding (format, shape)
and will be used with the fmflow.utils.readbinary function.
Format must be compatible with the Python's struct module.
For example, '>i' (int with big endian), or '10s' (10 chars).
For more information, see http://docs.python.jp/2/library/struct.html.

Available attributes:
- HEAD: Definition of record header info.
- CTL:  Definition of control info.
"""

from __future__ import absolute_import as __absolute_import
from __future__ import division as __division
from __future__ import print_function as __print_function

# the standard library
from collections import OrderedDict

# imported items
__all__ = ['HEAD', 'CTL']

# constants
ORDER = '<' # little endian


# definition of record header info
HEAD = OrderedDict()
HEAD['crec_type'] = (ORDER+'4s',)
HEAD['irec_len']  = (ORDER+'i',)


# definition of control info
CTL  = OrderedDict()
CTL['cversion']  = (ORDER+'8s',)
CTL['clog_name'] = (ORDER+'60s',)
CTL['cbe_type']  = (ORDER+'8s',)
