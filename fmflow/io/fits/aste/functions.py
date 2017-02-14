# coding: utf-8

"""Module for merging loggings of ASTE.

"""

# Python 3.x compatibility
from __future__ import absolute_import as _absolute_import
from __future__ import division as _division
from __future__ import print_function as _print_function

# the Python standard library
import re
import json

# the Python Package Index
import numpy as np
import yaml
import fmflow as fm
from astropy.io import fits
from astropy.coordinates import Angle

# imported items
__all__ = ['fromaste']

# constants
LAT_ASTE   = Angle('-22d58m17.69447s').deg
EFF_8257D  = 0.92 # exposure / interval time of Agilent 8257D


def fromaste(fmlolog, backendlog, antennalog=None, byteorder='<'):
    """Read logging data of ASTE and merge them into a FITS object.

    Args:
        fmlolog (str): File name of FMLO logging.
        backendlog (str): File name of backend logging.
        antennalog (str): File name of antenna logging (optional).
        byteorder (str):

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
    backend = _check_backend(backendlog, byteorder)

    if backend == 'AC45':
        hdu = _read_backendlog_mac(backendlog, byteorder)
    elif backend == 'FFX':
        # hdu = _read_backendlog_whsf(backendlog, byteorder)
        raise fm.utils.FMFlowError('WHSF logging is not supported yet')
    else:
        raise fm.utils..FMFlowError('invalid logging type')

    fitsobj.append(hdu)

    # ANTENNA HDU
    if antennalog is not None:
        hdu = _read_antennalog(antennalog)
        fitsobj.append(hdu)

    # OBSINFO HDU
    hdu = _make_obsinfo(fitsobj)
    fitsobj.insert(1, hdu)

    return fitsobj


def _read_fmlolog(fmlolog):
    """Read a FMLO logging of ASTE.

    Args:
        fmlolog (str): File name of FMLO logging.

    Returns:
        hdu (BinTableHDU): HDU containing the read FMLO logging.

    """
    # datetime converter
    c = fm.utils.DatetimeConverter('%Y%m%d%H%M%S.%f')
    starttime = np.loadtxt(fmlolog, str, usecols=(0,), skiprows=1, converters={0: c})
    scantype  = np.loadtxt(fmlolog, str, usecols=(1,), skiprows=1)
    fmfreq    = np.loadtxt(fmlolog, float, usecols=(2,), skiprows=1)
    lofreq    = np.loadtxt(fmlolog, float, usecols=(3,), skiprows=1)
    vrad      = np.loadtxt(fmlolog, float, usecols=(4,), skiprows=1)

    # fmflag
    fmflag = (scantype == 'ON')

    # bintable HDU
    header = fits.Header()
    header['EXTNAME'] = 'FMLOLOG'
    header['FILENAME'] = fmlolog

    columns = []
    columns.append(fits.Column('STARTTIME', 'A26', 'ISO 8601', array=starttime))
    columns.append(fits.Column('SCANTYPE',  'A4', '', array=scantype))
    columns.append(fits.Column('FMFLAG', 'L', 'bool', array=fmflag))
    columns.append(fits.Column('FMFREQ', 'D', 'Hz', array=fmfreq))
    columns.append(fits.Column('LOFREQ', 'D', 'Hz', array=lofreq))
    columns.append(fits.Column('VRAD', 'D', 'm/s', array=vrad))

    hdu = fits.BinTableHDU.from_columns(columns, header)
    return hdu


def _read_antennalog(antennalog):
    """Read an antenna logging of ASTE.

    Args:
        antennalog (str): File name of antenna logging.

    Returns:
        hdu (BinTableHDU): HDU containing the read antenna logging.

    """
    # datetime converter
    c = fm.utils.DatetimeConverter('%y%m%d%H%M%S.%f')
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
    q = -np.arcsin(degsin(az_prog)*degcos(LAT_ASTE)/degcos(dec_prog))

    ra_error  = -np.cos(q)*az_error + np.sin(q)*el_error
    dec_error = np.sin(q)*az_error + np.cos(q)*el_error
    ra_real   = ra_prog - ra_error
    dec_real  = dec_prog - ra_error

    # bintable HDU
    header = fits.Header()
    header['EXTNAME'] = 'ANTENNA'
    header['FILENAME'] = antennalog

    columns = []
    columns.append(fits.Column('STARTTIME', 'A26', 'ISO 8601', array=starttime))
    columns.append(fits.Column('RA', 'D', 'deg', array=ra_real))
    columns.append(fits.Column('DEC', 'D', 'deg', array=dec_real))
    columns.append(fits.Column('AZ', 'D', 'deg', array=az_real))
    columns.append(fits.Column('EL', 'D', 'deg', array=el_real))
    columns.append(fits.Column('RA_PROG', 'D', 'deg', array=ra_prog))
    columns.append(fits.Column('DEC_PROG', 'D', 'deg', array=dec_prog))
    columns.append(fits.Column('AZ_PROG', 'D', 'deg', array=az_prog))
    columns.append(fits.Column('EL_PROG', 'D', 'deg', array=el_prog))
    columns.append(fits.Column('RA_ERROR', 'D', 'deg', array=ra_error))
    columns.append(fits.Column('DEC_ERROR', 'D', 'deg', array=dec_error))
    columns.append(fits.Column('AZ_ERROR', 'D', 'deg', array=az_error))
    columns.append(fits.Column('EL_ERROR', 'D', 'deg', array=el_error))

    hdu = fits.BinTableHDU.from_columns(columns, header)
    return hdu


def _check_backend(backendlog, byteorder):
    """Check backend type from a backend logging of ASTE.

    Args:
        backendlog (str): File name of backend logging.
        byteorder (str):

    Returns:
        backend (str): Backend type.

    """
    with open('struct_common.yaml') as f:
         d = yaml.load(f)
         head = fm.utils.CStructReader(d['head'])
         ctl  = fm.utils.CStructReader(d['ctl'])

    with open(backendlog, 'rb') as f:
        head.read(f)
        ctl.read(f)
        backend = ctl.data['cbe_type']

    return backend


def _read_backendlog_mac(backendlog):
    """Read a backend logging of ASTE/MAC.

    Args:
        backendlog (str): File name of backend logging.
        byteorder (str):

    Returns:
        hdu (BinTableHDU): HDU containing the read backend logging.

    """
    with open('struct_common.yaml') as f:
         d = yaml.load(f)
         head = fm.utils.CStructReader(d['head'])
         ctl  = fm.utils.CStructReader(d['ctl'])

    with open('struct_mac.yaml') as f:
         d = yaml.load(f)
         obs = fm.utils.CStructReader(d['obs'])
         dat = fm.utils.CStructReader(d['dat'])

    prog = fm.utils.inprogress('reading backendlog', 100)

    def EOF(f):
        prog.next()
        head.read(f)
        flag = head._data['crec_type'][-1][0]
        return (flag == 'ED')

    # read data
    with open(backendlog, 'rb') as f:
        ## control info
        EOF(f)
        ctl.read(f)

        ## observation info
        EOF(f)
        obs.read(f)

        ## data info
        while not EOF(f):
            dat.read(f)

    # edit data
    data = dat.data
    data['STARTTIME'] = data.pop('cint_sttm')
    data['ARRAYID']   = data.pop('cary_name')
    data['SCANTYPE']  = data.pop('cscan_type')
    data['ARRAYDATA'] = data.pop('iary_data')

    ## starttime
    c = fm.utils.DatetimeConverter('%Y%m%d%H%M%S.%f')
    data['STARTTIME'] = np.array([c(t[:-2]) for t in data['STARTTIME']])

    ## arraydata
    ons  = fm.utils.where(data['SCANTYPE'] == 'ON')
    rs   = fm.utils.where(data['SCANTYPE'] == 'R')
    skys = fm.utils.where(data['SCANTYPE'] == 'SKY')

    ## apply scaling factor and offset
    arraydata = data['ARRAYDATA'].astype(float)
    arraydata *= data['dary_scf'][:,np.newaxis]
    arraydata += data['dary_offset'][:,np.newaxis]

    ## apply coeff. for ON
    for on in ons:
        for arrayid in np.unique(data['ARRAYID']):
            flag = (data['ARRAYID'][on] == arrayid)
            arraydata[on][flag] *= np.mean(data['dalpha'][on][flag])

    ## apply coeff. for R
    for (r, sky) in zip(rs, skys):
        arraydata[r] *= data['dbeta'][sky][:,np.newaxis]

    ## reverse array (if USB)
    usefg = np.array(obsinfo['iary_usefg'], dtype=bool)
    isusb = np.array(obsinfo['csid_type'])[usefg] == 'USB'
    for arrayid in np.unique(data['ARRAYID'])[isusb]:
        flag = (data['ARRAYID'] == arrayid)
        arraydata[flag] = arraydata[flag,::-1]

    data['ARRAYDATA'] = arraydata

    # read and edit formats
    fmts = dat.fitsformats
    fmts['STARTTIME'] = fmts.pop('cint_sttm')
    fmts['ARRAYID']   = fmts.pop('cary_name')
    fmts['SCANTYPE']  = fmts.pop('cscan_type')
    fmts['ARRAYDATA'] = fmts.pop('iary_data')
    fmts['STARTTIME'] = 'A26' # YYYY-mm-ddTHH:MM:SS.ssssss
    fmts['ARRAYDATA'] = '{[0]}E'.format(re.findall('\d+', fmts['ARRAYDATA']))

    # bintable HDU
    header = fits.Header()
    header['EXTNAME'] = 'BACKEND'
    header['FILENAME'] = backendlog
    header['CTLINFO'] = ctl.jsondata
    header['OBSINFO'] = obs.jsondata

    columns = []
    for key in data:
        columns.append(fits.Column(key, fmts[key], array=data[key]))

    hdu = fits.BinTableHDU.from_columns(columns, header)
    return hdu


def _make_obsinfo(fitsobj):
    """Make a OBSINFO HDU from FITS object.

    Args:
        fitsobj (HDUList): FITS object containing FMLOINFO, BACKEND HDUs.

    Returns:
        hdu (BinTableHDU): OBSINFO HDU containing the formatted observational info.

    """
    ctlinfo = json.loads(fitsobj['BACKEND'].header['CTLINFO'])
    obsinfo = json.loads(fitsobj['BACKEND'].header['OBSINFO'])
    datinfo = fitsobj['BACKEND'].data

    N = obsinfo['iary_num']
    flag = np.array(obsinfo['iary_usefg'], dtype=bool)

    fitstype  = 'FMFITS{}'.format(fm.__version__)
    frontend  = np.unique(np.array(obsinfo['cfe_type'])[flag])[0]
    backend   = ctlinfo['cbe_type']
    arrayid   = np.unique(datinfo['ARRAYID'])
    sideband  = np.array(obsinfo['csid_type'])[flag]
    interval  = np.tile(obsinfo['diptim'], N)
    exposure  = np.tile(obsinfo['diptim']*EFF_8257D, N)
    restfreq  = np.array(obsinfo['dcent_freq'])[flag]
    intmfreq  = np.array(obsinfo['dflif'])[flag]
    bandwidth = np.array(obsinfo['dbebw'])[flag]
    chanwidth = np.array(obsinfo['dbechwid'])[flag]

    if backend == 'AC45':
        numofchan = np.tile(obsinfo['ichanel'], N)
    elif backend == 'FFX':
        numofchan = np.array(obsinfo['ichanel'])

    # bintable HDU
    header = fits.Header()
    header['EXTNAME']  = 'OBSINFO'
    header['TELESCOP'] = 'ASTE'
    header['FITSTYPE'] = fitstype
    header['FRONTEND'] = frontend
    header['BACKEND']  = backend

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

