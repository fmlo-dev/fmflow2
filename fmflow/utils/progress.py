# coding: utf-8

# Python 3.x compatibility
from __future__ import absolute_import as __absolute_import
from __future__ import division as __division
from __future__ import print_function as __print_function

# the Python standard library
import sys

# imported items
__all__ = ['inprogress', 'progressbar']


def inprogress(message='in progress'):
    i = 0
    while True:
        status = '\r{} {}'.format(message, '|/-\\'[i%4])
        sys.stdout.write(status)
        sys.stdout.flush()
        i += 1
        yield


def progressbar(message='progress', N=100, barwidth=50):
    for i in range(N):
        prog = float(i+1)/N
        fill = '*' * int(barwidth*prog)
        space = ' ' * (barwidth-int(barwidth*prog))
        status = '\r{} [{}{}] {:.0%}'.format(message, fill, space, prog)
        sys.stdout.write(status)
        sys.stdout.flush()
        yield
