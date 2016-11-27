# coding: utf-8

from __future__ import absolute_import as _absolute_import
from __future__ import division as _division
from __future__ import print_function as _print_function

# the standard library
import os

# dependent packages¬
import numpy as np
from astropy.io import fits
from fmflow import utils as ut

# imported items
__all__ = ['getarray']


def getarray(fitsname, arrayid, scantype):
    """Create a modulated fmarray from a FMFITS.
    
    Args:
    - fitsname (str): File name of a FMFITS.
    - arrayid (str): An array ID with which the output fmarray is created.
    - scantype (str): A scan type with which the output fmarray is created.
    
    Returns:
    - fmarray (FMArray): An output fmarray of the spacified array ID and scan type.
    """
    from fmflow import fm
    with fits.open(os.path.expanduser(fitsname)) as f:
        oi = f['OBSINFO']
        fl = f['FMLOLOG']
        be = f['BACKEND']

        # flags
        flag_oi = oi.data.ARRAYID == arrayid
        flag_fl = fl.data.SCANTYPE == scantype
        flag_be = (be.data.ARRAYID == arrayid) & (be.data.SCANTYPE == scantype)

        # array
        array = np.squeeze(be.data.ARRAYDATA[flag_be])

        if scantype != 'ON':
            return array

        # info
        keys = [name.lower() for name in oi.data[flag_oi].names]
        values = oi.data[flag_oi][0]

        info = dict(zip(keys, values))
        info['fitstype'] = oi.header['FITSTYPE']
        info['telescop'] = oi.header['TELESCOP']
        info['frontend'] = oi.header['FRONTEND']
        info['backend']  = oi.header['BACKEND']

        # fmch
        t_fl = fl.data.STARTTIME[flag_fl]
        t_be = be.data.STARTTIME[flag_be]
        fmfreq = fl.data.FMFREQ[flag_fl]

        if len(set(t_be)-set(t_fl)) > 0:
            raise ut.FMFlowError('time range of FMLOLOG does not cover that of BACKEND')

        fmfreq_matched = []
        for t in t_be:
            fmfreq_matched.append(fmfreq[t_fl==t][0])

        fmch = (np.asarray(fmfreq_matched) / info['chanwidth']).astype(int)
        
        fmarray = fm.array(array, fmch, info=info)
        return fmarray
