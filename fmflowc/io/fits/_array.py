# coding: utf-8

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

# dependent libraries¬
import numpy as np
from astropy.io import fits
from fmflowc.utils import exceptions as e


def fromfits(fitsname, arrayid, scantype):
    '''Make a FMArray from a FMFITS.
    
    Args:
    - fitsname (str):
    - arrayid (str):
    - scantype (str):
    
    Returns
    - array (FMArray):
    '''
    from fmflowc import fm
    
    with fits.open(fitsname) as f:
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
        info = dict(zip(oi.data[flag_oi].names, oi.data[flag_oi][0]))
        info['FITSTYPE'] = oi.header['FITSTYPE']
        #info['TELESCOP'] = oi.header['TELESCOP']
        info['FRONTEND'] = oi.header['FRONTEND']
        info['BACKEND']  = oi.header['BACKEND']

        # fmch
        t_fl = fl.data.STARTTIME[flag_fl]
        t_be = be.data.STARTTIME[flag_be]
        fmfreq = fl.data.FMFREQ[flag_fl]

        if len(set(t_be)-set(t_fl)) > 0:
            e.FMflowError('time range of FMLOLOG does not cover that of BACKEND')

        fmfreq_matched = []
        for t in t_be:
            fmfreq_matched.append(fmfreq[t_fl==t][0])

        fmch = (np.asarray(fmfreq_matched) / info['CHANWIDTH']).astype(int)
        
        array = fm.FMArray(array, fmch, info)
        return array