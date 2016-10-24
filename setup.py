# coding: utf-8

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from setuptools import setup, find_packages


setup(
    name = 'fmflowc',
    version = '0.0.1',
    author = 'Akio Taniguchi',
    author_email = 'taniguchi@ioa.s.u-tokyo.ac.jp',
    description = 'Core Python Package for FMflow',
    url = 'https://github.com/snoopython/fmflow',
    keywords = 'astronomy radio fmlo',
    license = 'MIT',
    packages = find_packages(),
    install_requires = [
        'numpy>=1.10',
        'scipy>=0.15',
        'astropy>=1.0'
    ]
)
