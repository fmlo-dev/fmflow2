# coding: utf-8

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

# dependent libraries¬
import numpy as np
from astropy.io import fits
from fmflowc import fm
from fmflowc.utils import exceptions as e


def getarray(fitsname, arrayid, scantype):
    '''Make a FMArray from a FMFITS.
    
    Args:
    - fitsname (str):
    - arrayid (str):
    - scantype (str):
    
    Returns
    - fmarray (FMArray):
    '''
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
            raise e.FMflowError('time range of FMLOLOG does not cover that of BACKEND')

        fmfreq_matched = []
        for t in t_be:
            fmfreq_matched.append(fmfreq[t_fl==t][0])

        fmch = (np.asarray(fmfreq_matched) / info['chanwidth']).astype(int)
        
        array = fm.array(array, fmch, info=info)
        return array
