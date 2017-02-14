# coding: utf-8

# the Python Package Index
from setuptools import setup, find_packages

# local fmflow
from . import fmflow as fm

setup(
    name = 'fmflow',
    version = fm.__version__,
    author = 'Akio Taniguchi',
    author_email = 'taniguchi@ioa.s.u-tokyo.ac.jp',
    description = 'Core Python Package for FMFlow',
    url = 'https://github.com/snoopython/fmflow',
    keywords = 'astronomy radio fmlo',
    license = 'MIT',
    packages = find_packages(),
    install_requires = [
        'numpy>=1.11',
        'scipy>=0.18',
        'astropy>=1.2'
    ]
)
