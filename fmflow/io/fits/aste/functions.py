# coding: utf-8

"""Module for merging loggings of ASTE.

"""

# Python 3.x compatibility
from __future__ import absolute_import as _absolute_import
from __future__ import division as _division
from __future__ import print_function as _print_function

# the Python standard library
import os
import re
import json

# the Python Package Index
import yaml
import numpy as np
import fmflow as fm
from astropy import units as u
from astropy import constants as consts
from astropy.io import fits
from astropy.coordinates import Angle

# imported items
__all__ = ['fromaste']

# constants
C           = consts.c.value # spped of light in vacuum
D_ASTE      = (10.0 * u.m).value # diameter of the ASTE
LAT_ASTE    = Angle('-22d58m17.69447s').deg # latitude of the ASTE
EFF_8257D   = 0.92 # exposure / interval time of Agilent 8257D
DIR_MODULE  = os.path.dirname(os.path.realpath(__file__)) # module directory
IGNORED_KEY = '^[a-z]dmy([^_]|$)' # cdmy, cdmy2, ..., except for idmy_flag


def fromaste(fmlolog, backendlog, antennalog=None, byteorder='<'):
    """Read logging data of ASTE and merge them into a FITS object.

    Args:
        fmlolog (str): File name of FMLO logging.
        backendlog (str): File name of backend logging.
        antennalog (str): File name of antenna logging (optional).
        byteorder (str): format string that represents byte order
            of the backendlog. Default is '<' (little-endian).
            If the data in the returned FITS seems to be wrong,
            try to spacify '>' (big-endian).

    Returns:
        fitsobj (HDUlist): HDU list containing the merged data.

    See Also:
        https://docs.python.org/2/library/struct.html
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
        raise fm.utils.FMFlowError('invalid logging type')

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
    p = fm.utils.DatetimeParser()
    starttime = np.loadtxt(fmlolog, str, usecols=(0,), skiprows=1, converters={0: p})
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
    p = fm.utils.DatetimeParser()
    starttime = np.loadtxt(antennalog, str, usecols=(0,), skiprows=1, converters={0: p})
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
        byteorder (str): format string that represents byte order
            of the backendlog. Default is '<' (little-endian).
            If the data in the returned FITS seems to be wrong,
            try to spacify '>' (big-endian).

    Returns:
        backend (str): Backend type.

    """
    with open(os.path.join(DIR_MODULE, 'struct_common.yaml')) as f:
        d = yaml.load(f)
        head = fm.utils.CStructReader(d['head'], IGNORED_KEY, byteorder)
        ctl  = fm.utils.CStructReader(d['ctl'], IGNORED_KEY, byteorder)

    with open(backendlog, 'rb') as f:
        head.read(f)
        ctl.read(f)
        backend = ctl.data['cbe_type']

    return backend


def _read_backendlog_mac(backendlog, byteorder):
    """Read a backend logging of ASTE/MAC.

    Args:
        backendlog (str): File name of backend logging.
        byteorder (str): format string that represents byte order
            of the backendlog. Default is '<' (little-endian).
            If the data in the returned FITS seems to be wrong,
            try to spacify '>' (big-endian).

    Returns:
        hdu (BinTableHDU): HDU containing the read backend logging.

    """
    with open(os.path.join(DIR_MODULE, 'struct_common.yaml')) as f:
        d = yaml.load(f)
        head = fm.utils.CStructReader(d['head'], IGNORED_KEY, byteorder)
        ctl  = fm.utils.CStructReader(d['ctl'], IGNORED_KEY, byteorder)

    with open(os.path.join(DIR_MODULE, 'struct_mac.yaml')) as f:
        d = yaml.load(f)
        obs = fm.utils.CStructReader(d['obs'], IGNORED_KEY, byteorder)
        dat = fm.utils.CStructReader(d['dat'], IGNORED_KEY, byteorder)

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

        print('done.')

    # edit data
    print('post processing', end=' ')

    data = dat.data
    data['STARTTIME'] = data.pop('cint_sttm')
    data['ARRAYID']   = data.pop('cary_name')
    data['SCANTYPE']  = data.pop('cscan_type')
    data['ARRAYDATA'] = data.pop('iary_data')

    ## starttime
    p = fm.utils.DatetimeParser()
    data['STARTTIME'] = np.array([p(t) for t in data['STARTTIME']])

    ## scantype (bug?)
    data['SCANTYPE'][data['SCANTYPE']=='R\x00RO'] = 'R'

    ## arraydata
    arraydata = data['ARRAYDATA'].astype(float)

    ## apply scaling factor and offset
    arraydata *= data['dary_scf'][:,np.newaxis]
    arraydata += data['dary_offset'][:,np.newaxis]

    ## slices of each scantype of arraydata
    ons  = fm.utils.slicewhere(data['SCANTYPE'] == 'ON')
    rs   = fm.utils.slicewhere(data['SCANTYPE'] == 'R')
    skys = fm.utils.slicewhere(data['SCANTYPE'] == 'SKY')
    zero = fm.utils.slicewhere(data['SCANTYPE'] == 'ZERO')[0]

    ## apply ZERO and coeff. to ON data
    for on in ons:
        for aid in np.unique(data['ARRAYID']):
            aid_on   = (data['ARRAYID'][on] == aid)
            aid_zero = (data['ARRAYID'][zero] == aid)
            arraydata[on][aid_on] -= arraydata[zero][aid_zero]
            arraydata[on][aid_on] *= np.mean(data['dalpha'][on][aid_on])

    ## apply ZERO and coeff. to R data
    for (r, sky) in zip(rs, skys):
        for i, aid in enumerate(np.unique(data['ARRAYID'])):
            aid_r = (data['ARRAYID'][r] == aid)
            aid_sky  = (data['ARRAYID'][sky] == aid)
            aid_zero = (data['ARRAYID'][zero] == aid)
            arraydata[r][aid_r] -= arraydata[zero][aid_zero]
            arraydata[r][aid_r] *= data['dbeta'][sky][aid_sky]

    ## apply ZERO to SKY data
    for sky in skys:
        for aid in np.unique(data['ARRAYID']):
            aid_sky  = (data['ARRAYID'][sky] == aid)
            aid_zero = (data['ARRAYID'][zero] == aid)
            arraydata[sky][aid_sky] -= arraydata[zero][aid_zero]

    ## reverse array (if USB)
    usefg = np.array(obs.data['iary_usefg'], dtype=bool)
    isusb = np.array(obs.data['csid_type'])[usefg] == 'USB'
    for arrayid in np.unique(data['ARRAYID'])[isusb]:
        flag = (data['ARRAYID'] == arrayid)
        arraydata[flag] = arraydata[flag,::-1]

    ## finally
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

    # bintable HDU header
    p = fm.utils.DatetimeParser()

    header = fits.Header()
    header['EXTNAME']  = 'OBSINFO'
    header['FITSTYPE'] = 'FMFITS{}'.format(fm.__version__)
    header['TELESCOP'] = 'ASTE'
    header['DATE-OBS'] = p(obsinfo['clog_id'])[:-3]
    header['OBSERVER'] = obsinfo['cobs_user']
    header['OBJECT']   = obsinfo['cobj_name']
    header['RA']       = obsinfo['dsrc_pos'][0][0]
    header['DEC']      = obsinfo['dsrc_pos'][1][0]
    header['EQUINOX']  = float(re.findall('\d+', obsinfo['cepoch'])[0])

    # bintable HDU data
    N = obsinfo['iary_num']
    flag = np.array(obsinfo['iary_usefg'], dtype=bool)

    arrayid   = np.unique(datinfo['ARRAYID'])
    sideband  = np.array(obsinfo['csid_type'])[flag]
    restfreq  = np.array(obsinfo['dcent_freq'])[flag]
    intmfreq  = np.array(obsinfo['dflif'])[flag]
    beamsize  = np.rad2deg(1.2*C/D_ASTE) / restfreq
    bandwidth = np.array(obsinfo['dbebw'])[flag]
    chanwidth = np.array(obsinfo['dbechwid'])[flag]
    interval  = np.tile(obsinfo['diptim'], N)
    integtime = np.tile(obsinfo['diptim']*EFF_8257D, N)
    frontend  = np.array(obsinfo['cfe_type'])[flag]
    backend   = np.tile(ctlinfo['cbe_type'], N)

    if backend[0] == 'AC45':
        numofchan = np.tile(obsinfo['ichanel'], N)
    elif backend[0] == 'FFX':
        numofchan = np.array(obsinfo['ichanel'])[flag]
    else:
        raise fm.utils.FMFlowError('invalid logging type')

    columns = []
    columns.append(fits.Column('ARRAYID', 'A4', '', array=arrayid))
    columns.append(fits.Column('SIDEBAND', 'A4', '', array=sideband))
    columns.append(fits.Column('RESTFREQ', 'D', 'Hz', array=restfreq))
    columns.append(fits.Column('INTMFREQ', 'D', 'Hz', array=intmfreq))
    columns.append(fits.Column('BEAMSIZE', 'D', 'deg', array=beamsize))
    columns.append(fits.Column('BANDWIDTH', 'D', 'Hz', array=bandwidth))
    columns.append(fits.Column('CHANWIDTH', 'D', 'Hz', array=chanwidth))
    columns.append(fits.Column('NUMOFCHAN', 'J', 'ch', array=numofchan))
    columns.append(fits.Column('INTERVAL', 'D', 's', array=interval))
    columns.append(fits.Column('INTEGTIME', 'D', 's', array=integtime))
    columns.append(fits.Column('FRONTEND', 'A10', '', array=frontend))
    columns.append(fits.Column('BACKEND', 'A8', '', array=backend))

    hdu = fits.BinTableHDU.from_columns(columns, header)
    return hdu

