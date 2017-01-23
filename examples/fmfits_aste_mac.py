# coding: utf-8

"""An example of creating a FMFITS from ASTE/MAC FMLO loggings.

Author: Akio Taniguchi
FMFlow version: 0.1
"""

from fmflow.io import fits

fmlolog    = 'FMLOLOG.yyyymmddhhmmss'
antennalog = 'ANTENNALOG.yyyymmddhhmmss'
backendlog = 'ACG.<obstable>.<user>.<proj>.yyyymmddhhmmss.s'

f = fits.fromaste(fmlolog, backendlog, antennalog)
f.writeto('example.fits')
