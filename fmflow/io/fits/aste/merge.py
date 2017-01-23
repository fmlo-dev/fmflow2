# coding: utf-8

"""Module for merging loggings of ASTE.

"""

# Python 3.x compatibility
from __future__ import absolute_import as __absolute_import
from __future__ import division as __division
from __future__ import print_function as __print_function

# the Python standard library
import re
import json
from collections import OrderedDict

# the Python Package Index
import numpy as np
from astropy.io import fits
from astropy.coordinates import Angle
from fmflow import utils as ut

# imported items
__all__ = ['fromaste']

# constants
LAT = Angle('-22d58m17.69447s').deg # latitude of ASTE
EFF = 0.92 # exposure time / interval time of Agilent 8257D


def fromaste(fmlolog, backendlog, antennalog=None):
    """Read logging data of ASTE and merge them into a FITS object.

    Args:
        fmlolog (str): File name of FMLO logging.
        backendlog (str): File name of backend logging.
        antennalog (str): File name of antenna logging (optional).

    Returns:
        fitsobj (HDUlist): HDU list containing the merged data.

    """
    # HDU list
    fitsobj = fits.HDUList()

    # PRIMARY HDU
    hdu = fits.PrimaryHDU()
    fitsobj.append(hdu)

    # FMLOINFO HDU
    hdu = _read_fmlolog(fmlolog)
    fitsobj.append(hdu)

    # BACKEND HDU
    backend = _check_backend(backendlog)

    if backend == 'AC45':
        hdu = _read_backendlog_mac(backendlog)
    elif backend == 'FFX':
        # hdu = _read_backendlog_whsf(backendlog)
        raise ut.FMFlowError('WHSF logging is not supported yet')
    else:
        raise ut.FMFlowError('invalid logging type')

    fitsobj.append(hdu)

    # ANTENNA HDU
    if antennalog is not None:
        hdu = _read_antennalog(antennalog)
        fitsobj.append(hdu)

    # OBSINFO HDU
    hdu = _make_obsinfo(fitsobj)
    fitsobj.insert(1, hdu)

    return fitsobj


def _make_obsinfo(fitsobj):
    """Make a OBSINFO HDU from FITS object.

    Args:
        fitsobj (HDUList): FITS object containing FMLOINFO, BACKEND HDUs.

    Returns:
        hdu (BinTableHDU): OBSINFO HDU containing the formatted observational info.

    """
    d_ctl = json.loads(fitsobj['BACKEND'].header['CTLINFO'])
    d_obs = json.loads(fitsobj['BACKEND'].header['OBSINFO'])

    N = d_obs['iary_num']
    flag = np.array(d_obs['iary_usefg'], dtype=bool)
    backend = d_ctl['cbe_type']

    header = fits.Header()
    header['EXTNAME']  = 'OBSINFO'
    header['FITSTYPE'] = 'FMFITSv1'
    header['TELESCOP'] = 'ASTE'
    header['FRONTEND'] = np.unique(np.array(d_obs['cfe_type'])[flag])[0]
    header['BACKEND']  = backend

    arrayid   = np.array(['A{}'.format(i+1) for i in range(len(flag))])[flag]
    sideband  = np.array(d_obs['csid_type'])[flag]
    interval  = np.tile(d_obs['diptim'], N)
    exposure  = np.tile(d_obs['diptim']*EFF, N)
    restfreq  = np.array(d_obs['dcent_freq'])[flag]
    intmfreq  = np.array(d_obs['dflif'])[flag]
    bandwidth = np.array(d_obs['dbebw'])[flag]
    chanwidth = np.array(d_obs['dbechwid'])[flag]

    if backend == 'AC45':
        numofchan = np.tile(d_obs['ichanel'], N)
    elif backend == 'FFX':
        numofchan = np.array(d_obs['ichanel'])

    # bintable HDU
    columns = []
    columns.append(fits.Column('ARRAYID', 'A3', '', array=arrayid))
    columns.append(fits.Column('SIDEBAND', 'A3', '', array=sideband))
    columns.append(fits.Column('INTERVAL', 'D', 's', array=interval))
    columns.append(fits.Column('EXPOSURE', 'D', 's', array=exposure))
    columns.append(fits.Column('RESTFREQ', 'D', 'Hz', array=restfreq))
    columns.append(fits.Column('INTMFREQ', 'D', 'Hz', array=intmfreq))
    columns.append(fits.Column('BANDWIDTH', 'D', 'Hz', array=bandwidth))
    columns.append(fits.Column('CHANWIDTH', 'D', 'Hz', array=chanwidth))
    columns.append(fits.Column('NUMOFCHAN', 'J', 'ch', array=numofchan))

    hdu = fits.BinTableHDU.from_columns(columns, header)
    return hdu


def _read_fmlolog(fmlolog):
    """Read a FMLO logging of ASTE.

    Args:
        fmlolog (str): File name of FMLO logging.

    Returns:
        hdu (BinTableHDU): HDU containing the read FMLO logging.

    """
    header = fits.Header()
    header['EXTNAME'] = 'FMLOLOG'
    header['FILENAME'] = fmlolog

    # datetime converter
    c = ut.DatetimeConverter('%Y%m%d%H%M%S.%f')

    starttime = np.loadtxt(fmlolog, str, usecols=(0,), skiprows=1, converters={0: c})
    scantype  = np.loadtxt(fmlolog, str, usecols=(1,), skiprows=1)
    fmfreq    = np.loadtxt(fmlolog, float, usecols=(2,), skiprows=1)
    lofreq    = np.loadtxt(fmlolog, float, usecols=(3,), skiprows=1)
    vrad      = np.loadtxt(fmlolog, float, usecols=(4,), skiprows=1)

    # fmflag
    fmflag = (scantype == 'ON')

    # bintable HDU
    columns = []
    columns.append(fits.Column('STARTTIME', 'A26', 'ISO 8601', array=starttime))
    columns.append(fits.Column('SCANTYPE',  'A4', '', array=scantype))
    columns.append(fits.Column('FMFLAG',    'L', 'bool', array=fmflag))
    columns.append(fits.Column('FMFREQ',    'D', 'Hz', array=fmfreq))
    columns.append(fits.Column('LOFREQ',    'D', 'Hz', array=lofreq))
    columns.append(fits.Column('VRAD',      'D', 'm/s', array=vrad))

    hdu = fits.BinTableHDU.from_columns(columns, header)
    return hdu


def _read_antennalog(antennalog):
    """Read an antenna logging of ASTE.

    Args:
        antennalog (str): File name of antenna logging.

    Returns:
        hdu (BinTableHDU): HDU containing the read antenna logging.

    """
    header = fits.Header()
    header['EXTNAME'] = 'ANTENNA'
    header['FILENAME'] = antennalog

    # datetime converter
    c = ut.DatetimeConverter('%y%m%d%H%M%S.%f')

    starttime = np.loadtxt(antennalog, str, usecols=(0,), skiprows=1, converters={0: c})
    ra_prog   = np.loadtxt(antennalog, float, usecols=(1,), skiprows=1)
    dec_prog  = np.loadtxt(antennalog, float, usecols=(2,), skiprows=1)
    az_prog   = np.loadtxt(antennalog, float, usecols=(3,), skiprows=1)
    el_prog   = np.loadtxt(antennalog, float, usecols=(4,), skiprows=1)
    az_real   = np.loadtxt(antennalog, float, usecols=(5,), skiprows=1)
    el_real   = np.loadtxt(antennalog, float, usecols=(6,), skiprows=1)
    az_error  = np.loadtxt(antennalog, float, usecols=(7,), skiprows=1)
    el_error  = np.loadtxt(antennalog, float, usecols=(8,), skiprows=1)

    # RA,Dec real
    degsin = lambda deg: np.sin(np.deg2rad(deg))
    degcos = lambda deg: np.cos(np.deg2rad(deg))
    q = -np.arcsin(degsin(az_prog)*degcos(LAT)/degcos(dec_prog))

    ra_error  = -np.cos(q)*az_error + np.sin(q)*el_error
    dec_error = np.sin(q)*az_error + np.cos(q)*el_error
    ra_real   = ra_prog - ra_error
    dec_real  = dec_prog - ra_error

    # bintable HDU
    columns = []
    columns.append(fits.Column('STARTTIME', 'A26', 'ISO 8601', array=datetime))
    columns.append(fits.Column('RA',        'D', 'deg', array=ra_real))
    columns.append(fits.Column('DEC',       'D', 'deg', array=dec_real))
    columns.append(fits.Column('AZ',        'D', 'deg', array=az_real))
    columns.append(fits.Column('EL',        'D', 'deg', array=el_real))
    columns.append(fits.Column('RA_PROG',   'D', 'deg', array=ra_prog))
    columns.append(fits.Column('DEC_PROG',  'D', 'deg', array=dec_prog))
    columns.append(fits.Column('AZ_PROG',   'D', 'deg', array=az_prog))
    columns.append(fits.Column('EL_PROG',   'D', 'deg', array=el_prog))
    columns.append(fits.Column('RA_ERROR',  'D', 'deg', array=ra_error))
    columns.append(fits.Column('DEC_ERROR', 'D', 'deg', array=dec_error))
    columns.append(fits.Column('AZ_ERROR',  'D', 'deg', array=az_error))
    columns.append(fits.Column('EL_ERROR',  'D', 'deg', array=el_error))

    hdu = fits.BinTableHDU.from_columns(columns, header)
    return hdu


def _check_backend(backendlog):
    """Check backend type from a backend logging of ASTE.

    Args:
        backendlog (str): File name of backend logging.

    Returns:
        backend (str): Backend type.

    """
    from .otflog_common import HEAD, CTL

    with open(backendlog, 'rb') as f:
        ut.readbinary(f, HEAD)
        readdata = ut.readbinary(f, CTL)

    backend = readdata['cbe_type']
    return backend


def _read_backendlog_mac(backendlog):
    """Read a backend logging of ASTE/MAC.

    Args:
        backendlog (str): File name of backend logging.

    Returns:
        hdu (BinTableHDU): HDU containing the read backend logging.

    """
    header = fits.Header()
    header['EXTNAME'] = 'BACKEND'
    header['FILENAME'] = backendlog

    from .otflog_common import HEAD, CTL
    from .otflog_mac import OBS, DAT

    prog = ut.inprogress('reading the backendlog', 50)

    def EOF(f):
        prog.next()
        readdata = ut.readbinary(f, HEAD)
        return readdata['crec_type'] == 'ED'

    # read data
    with open(backendlog, 'rb') as f:
        # control info
        ut.readbinary(f, HEAD)
        readdata = ut.readbinary(f, CTL)
        header['CTLINFO'] = json.dumps(readdata)

        # observation info
        ut.readbinary(f, HEAD)
        readdata = ut.readbinary(f, OBS)
        header['OBSINFO'] = json.dumps(readdata)

        # data info
        data = OrderedDict()
        for key in DAT:
            if not re.search('dmy', key):
                data[key] = []

        while not EOF(f):
            readdata = ut.readbinary(f, DAT)
            for key in data:
                data[key].append(readdata[key])

    # read formats
    fmts = ut.getfitsformat(DAT)

    # edit data and formats
    for key in data:
        data[key] = np.asarray(data[key])

    data['STARTTIME'] = data.pop('cint_sttm')
    data['ARRAYID']   = data.pop('cary_name')
    data['SCANTYPE']  = data.pop('cscan_type')
    data['ARRAYDATA'] = data.pop('iary_data')

    fmts['STARTTIME'] = fmts.pop('cint_sttm')
    fmts['ARRAYID']   = fmts.pop('cary_name')
    fmts['SCANTYPE']  = fmts.pop('cscan_type')
    fmts['ARRAYDATA'] = fmts.pop('iary_data')

    # edit starttime
    c = ut.DatetimeConverter('%Y%m%d%H%M%S.%f')
    data['STARTTIME'] = np.asarray([c(t[:-2]) for t in data['STARTTIME']])
    fmts['STARTTIME'] = 'A26' # 26 characters in FITS

    # edit arraydata
    ons  = ut.where(data['SCANTYPE'] == 'ON')
    rs   = ut.where(data['SCANTYPE'] == 'R')
    skys = ut.where(data['SCANTYPE'] == 'SKY')

    arraydata = data['ARRAYDATA'].astype(float)
    arraydata *= data['dary_scf'][:,np.newaxis]
    arraydata += data['dary_offset'][:,np.newaxis]

    for on in ons:
        for arrayid in np.unique(data['ARRAYID']):
            flag = (data['ARRAYID'][on] == arrayid)
            arraydata[on][flag] *= np.mean(data['dalpha'][on][flag])

    for (r, sky) in zip(rs, skys):
        arraydata[r] *= data['dbeta'][sky][:,np.newaxis]

    data['ARRAYDATA'] = arraydata
    fmts['ARRAYDATA'] = re.findall('\d+', fmts['ARRAYDATA'])[0] + 'D'

    # bintable HDU
    columns = []
    for key in data:
        columns.append(fits.Column(key, fmts[key], array=data[key]))

    hdu = fits.BinTableHDU.from_columns(columns, header)
    return hdu

