# coding: utf-8

"""An example of creating a FMFITS from ASTE/MAC FMLO loggings.

Author: Akio Taniguchi
FMFlow version: 0.2
"""

import fmflow as fm

fmlolog    = 'FMLOLOG.yyyymmddhhmmss' # replace it with yours
antennalog = 'ANTENNALOG.yyyymmddhhmmss' # replace it with yours
backendlog = 'ACG.obstable.user.proj.yyyymmddhhmmss.s' # replace it with yours
fitsname   = 'saved_fmfits.fits' # replace it with yours

f = fm.io.fits.fromaste(fmlolog, backendlog, antennalog)
f.writeto(fitsname)
