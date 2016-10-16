# coding: utf-8

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

# standard libraries
from collections import OrderedDict


# record header info
HEAD = OrderedDict()
HEAD['crec_type'] = ('<4s',)
HEAD['irec_len']  = ('<i',)

# control info
CTL = OrderedDict()
CTL['cversion']  = ('<8s',)
CTL['clog_name'] = ('<60s',)
CTL['cbe_type']  = ('<8s',)

