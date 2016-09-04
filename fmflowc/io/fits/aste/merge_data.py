# coding: utf-8

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

# standard libraries
import re
import json
from datetime import datetime as dt
from struct import calcsize, unpack
from collections import OrderedDict

# dependent libraries
import numpy as np
from astropy.io import fits
from fmflowc.utils import binary


def fromaste(speclog, fmlolog, antlog=None, fitsname=None, spectrometer=None):
    '''Merge data taken with ASTE and save them in a FMFITS.

    Args:
    - speclog (str): name of spectrometer (MAC or WHSF) logging file.
    - fmlolog (str): name of FMLO logging file.
    - antlog (str): name of antenna logging file (optional).
    - fitsname (str): name of output fits file (optional).
        if not spacified, then <speclog>.fits will be created.
    - spectrometer (str): name of spectromter used in the observation.
        options: MAC or WHSF . if not spacified, then the function will try
        to detect the spectrometer type from speclog (i.e. ACG or FFX).

    Returns:
    - None (NoneType): this function returns nothing.

    '''
    if spectrometer is None:
        prefix = speclog.split('.')[0]
        if prefix == 'ACG':
            spectrometer = 'MAC'
        elif prefix == 'FFX':
            spectrometer = 'WHSF'
        else:
            message = 'spectrometer type is not detected from {}'.format(speclog)
            raise ValueError, message
    else:
        if spectrometer not in ['MAC', 'WHSF']:
            raise ValueError, spectrometer

    if fitsname is None:
        fitsname = speclog + '.fits'

    # hdu list
    hdulist = fits.HDUList()

    # primary hdu
    header = fits.Header()
    header['FITSTYPE'] = 'FMFITSv0'
    header['MERGEDAT'] = dt.now().isoformat()
    header['TELESCOP'] = 'ASTE'
    header['SPECTROM'] = spectrometer
    hdu = fits.PrimaryHDU(header=header)
    hdulist.append(hdu)

    # speclog hdu
    if spectrometer == 'MAC':
        from .mac_otf_def import HEAD, CTL, OBS, DAT
    elif spectrometer == 'WHSF':
        from .whsf_otf_def import HEAD, CTL, OBS, DAT

    hdu = _load_speclog(speclog, HEAD, CTL, OBS, DAT)
    hdulist.append(hdu)

    # fmlolog hdu
    hdu = _load_fmlolog(fmlolog)
    hdulist.append(hdu)

    # antlog hdu (optional)
    if antlog is not None:
        hdu = _load_antlog(antlog)
        hdulist.append(hdu)

    # output FITS
    try:
        hdulist.writeto(fitsname)
    except IOError as e:
        message = '{}\nDo you want to overwrite it? [y/n] '.format(e)
        if raw_input(message).lower() == 'y':
            hdulist.writeto(fitsname, clobber=True)


def _load_speclog(speclog, HEAD, CTL, OBS, DAT):
    '''

    '''
    header = fits.Header()
    header['EXTNAME'] = 'SPECLOG'
    header['ORGFILE'] = speclog

    with open(speclog, 'rb') as f:
        # control info
        binary.dumpstruct(f, HEAD)
        struct = binary.dumpstruct(f, CTL)
        header['CTLINFO'] = json.dumps(struct)

        # observation info
        binary.dumpstruct(f, HEAD)
        struct = binary.dumpstruct(f, OBS)
        header['OBSINFO'] = json.dumps(struct)

        # data info
        data = OrderedDict()
        for key in DAT:
            if not re.search('dmy', key):
                data[key] = []

        def EOF(f):
            struct = binary.dumpstruct(f, HEAD)
            return struct['crec_type'] == 'ED'

        while not EOF(f):
            struct = binary.dumpstruct(f, DAT)
            for key in data:
                data[key].append(struct[key])

    columns = []
    struct = binary.fitsformat(DAT)
    for key in data:
        columns.append(fits.Column(key, struct[key], array=data[key]))

    hdu = fits.BinTableHDU.from_columns(columns, header)
    return hdu


def _load_fmlolog(fmlolog):
    '''

    '''
    header = fits.Header()
    header['EXTNAME'] = 'FMLOLOG'
    header['ORGFILE'] = fmlolog

    datetime = np.loadtxt(fmlolog, dtype=str, usecols=(0,), skiprows=1)
    scantype = np.loadtxt(fmlolog, dtype=str, usecols=(1,), skiprows=1)
    fmfreq   = np.loadtxt(fmlolog, dtype=float, usecols=(2,), skiprows=1)
    lofreq   = np.loadtxt(fmlolog, dtype=float, usecols=(3,), skiprows=1)
    vrad     = np.loadtxt(fmlolog, dtype=float, usecols=(4,), skiprows=1)

    columns = []
    columns.append(fits.Column('datetime', 'A19', '%Y%m%d%H%M%S.%f', array=datetime))
    columns.append(fits.Column('scantype', 'A4', '', array=scantype))
    columns.append(fits.Column('fmfreq', 'D', 'Hz', array=fmfreq))
    columns.append(fits.Column('lofreq', 'D', 'Hz', array=lofreq))
    columns.append(fits.Column('vrad', 'D', 'm/s', array=vrad))

    hdu = fits.BinTableHDU.from_columns(columns, header)

    return hdu


def _load_antlog(antlog):
    '''

    '''
    header = fits.Header()
    header['EXTNAME'] = 'ANTLOG'
    header['ORGFILE'] = antlog

    datetime = np.loadtxt(antlog, dtype=str, usecols=(0,), skiprows=1)
    ra       = np.loadtxt(antlog, dtype=float, usecols=(1,), skiprows=1)
    dec      = np.loadtxt(antlog, dtype=float, usecols=(2,), skiprows=1)
    azprog   = np.loadtxt(antlog, dtype=float, usecols=(3,), skiprows=1)
    elprog   = np.loadtxt(antlog, dtype=float, usecols=(4,), skiprows=1)
    azcoll   = np.loadtxt(antlog, dtype=float, usecols=(5,), skiprows=1)
    elcoll   = np.loadtxt(antlog, dtype=float, usecols=(6,), skiprows=1)
    azopterr = np.loadtxt(antlog, dtype=float, usecols=(7,), skiprows=1)
    elopterr = np.loadtxt(antlog, dtype=float, usecols=(8,), skiprows=1)

    columns = []
    columns.append(fits.Column('datetime', 'A15', '%y%m%d%H%M%S.%f', array=datetime))
    columns.append(fits.Column('ra', 'D', 'deg', array=ra))
    columns.append(fits.Column('dec', 'D', 'deg', array=dec))
    columns.append(fits.Column('azprog', 'D', 'deg', array=azprog))
    columns.append(fits.Column('elprog', 'D', 'deg', array=elprog))
    columns.append(fits.Column('azcoll', 'D', 'deg', array=azcoll))
    columns.append(fits.Column('elcoll', 'D', 'deg', array=elcoll))
    columns.append(fits.Column('azopterr', 'D', 'deg', array=azopterr))
    columns.append(fits.Column('elopterr', 'D', 'deg', array=elopterr))

    hdu = fits.BinTableHDU.from_columns(columns, header)

    return hdu
