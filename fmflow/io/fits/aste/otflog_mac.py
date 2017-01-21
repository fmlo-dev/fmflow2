# coding: utf-8

"""Module for defining structures of OTF logging (ASTE/MAC).

Structure is an ordered dictionary (collections.OrderedDict),
which defines each name of parameter and corresponding (format, shape)
and will be used with the fmflow.utils.readbinary function.
Format must be compatible with the Python's struct module.
For example, '>i' (int with big endian), or '10s' (10 chars).
For more information, see http://docs.python.jp/2/library/struct.html.

Available attributes:
- OBS: Definition of record header info.
- DAT: Definition of control info.
"""

from __future__ import absolute_import as __absolute_import
from __future__ import division as __division
from __future__ import print_function as __print_function

# the standard library
from collections import OrderedDict

# constants
ORDER = '<'    # little endian
ARYMAX = 8     # maximum number of arrays of spectrometer
CHMAX  = 1024  # maxinum number of channels of each array


# definition of observation info
OBS = OrderedDict()
OBS['cgroup']          = (ORDER+'8s',)
OBS['cproject']        = (ORDER+'16s',)
OBS['cobs_file']       = (ORDER+'24s',)
OBS['cobs_user']       = (ORDER+'40s',)
OBS['cst_time']        = (ORDER+'16s',)
OBS['clog_id']         = (ORDER+'24s',)
OBS['cobj_name']       = (ORDER+'16s',)
OBS['csrc_comment']    = (ORDER+'24s',)
OBS['cepoch']          = (ORDER+'8s',)
OBS['cseq_ptn']        = (ORDER+'120s',)
OBS['clagwin']         = (ORDER+'8s',)
OBS['ctrk_type']       = (ORDER+'8s',)
OBS['cscan_coord']     = (ORDER+'8s',)
OBS['cdmy1']           = (ORDER+'24s',) # reserved
OBS['icoord_type']     = (ORDER+'i',)
OBS['iscan_coord']     = (ORDER+'i',)
OBS['isw_mod']         = (ORDER+'i',)
OBS['icalib_int']      = (ORDER+'i',)
OBS['ichanel']         = (ORDER+'i',)
OBS['iary_num']        = (ORDER+'i',)
OBS['ivdef']           = (ORDER+'i',)
OBS['ivref']           = (ORDER+'i',)
OBS['imult_flag']      = (ORDER+'i',)
OBS['ich_bind']        = (ORDER+'i',)
OBS['ich_range']       = (ORDER+'i', 2)
OBS['iary_usefg']      = (ORDER+'i', ARYMAX)
OBS['iintensity_flag'] = (ORDER+'i',)
OBS['idmy1']           = (ORDER+'i', 9,) # reserved
OBS['dsrc_pos']        = (ORDER+'d', (2,3))
OBS['dsub_off']        = (ORDER+'d', 3)
OBS['dant_off']        = (ORDER+'d', 2)
OBS['dstrv']           = (ORDER+'d',)
OBS['dbsw_freq']       = (ORDER+'d',)
OBS['dbsw_ofd']        = (ORDER+'d',)
OBS['dmap_pos']        = (ORDER+'d',)
OBS['dalctime']        = (ORDER+'d',)
OBS['dmult_off']       = (ORDER+'d',)
OBS['diptim']          = (ORDER+'d',)
OBS['dsb_factor']      = (ORDER+'d', ARYMAX)
OBS['ccmt_t']          = (ORDER+'24s',)
OBS['dcmt_q']          = (ORDER+'d',)
OBS['dcmt_e']          = (ORDER+'d',)
OBS['dcmt_somega']     = (ORDER+'d',)
OBS['dcmt_node']       = (ORDER+'d',)
OBS['dcmt_i']          = (ORDER+'d',)
OBS['dant_vel']        = (ORDER+'d', ARYMAX)
OBS['dcent_freq']      = (ORDER+'d', ARYMAX)
OBS['dtrk_freq']       = (ORDER+'d', ARYMAX)
OBS['dflif']           = (ORDER+'d', ARYMAX)
OBS['dmult_scf']       = (ORDER+'d', ARYMAX)
OBS['dbebw']           = (ORDER+'d', ARYMAX)
OBS['dberes']          = (ORDER+'d', ARYMAX)
OBS['dbechwid']        = (ORDER+'d', ARYMAX)
OBS['dsw_freq']        = (ORDER+'d', ARYMAX)
OBS['dhpbw']           = (ORDER+'d', ARYMAX)
OBS['deffa']           = (ORDER+'d', ARYMAX)
OBS['deffb']           = (ORDER+'d', ARYMAX)
OBS['deffl']           = (ORDER+'d', ARYMAX)
OBS['defss']           = (ORDER+'d', ARYMAX)
OBS['dgain']           = (ORDER+'d', ARYMAX)
OBS['dpoldr']          = (ORDER+'d', ARYMAX)
OBS['ipordr']          = (ORDER+'i', ARYMAX)
OBS['iipint']          = (ORDER+'i', ARYMAX)
OBS['iary_refno']      = (ORDER+'i', ARYMAX)
OBS['imlt_no']         = (ORDER+'i', ARYMAX)
OBS['cfe_type']        = (ORDER+'10s', ARYMAX)
OBS['csid_type']       = (ORDER+'4s', ARYMAX)
OBS['chorn_typ']       = (ORDER+'2s', ARYMAX)
OBS['cpol_typ']        = (ORDER+'5s', ARYMAX)
OBS['cband_type']      = (ORDER+'2s', ARYMAX)
OBS['cposition']       = (ORDER+'32s',)
OBS['ctelescope']      = (ORDER+'16s',)
OBS['cdmy2']           = ('16s',) # reserved
OBS['ifqdat_nf']       = (ORDER+'i', ARYMAX)
OBS['dfqdat_f0']       = (ORDER+'d', ARYMAX)
OBS['dfqdat_fq']       = (ORDER+'d', (ARYMAX,10))
OBS['dfqdat_ch']       = (ORDER+'d', (ARYMAX,10))
OBS['dfqdat_cw']       = (ORDER+'d', (ARYMAX,10))

# definition of data info
DAT = OrderedDict()
DAT['cint_sttm']   = (ORDER+'24s',)
DAT['cscan_type']  = (ORDER+'8s',)
DAT['cary_name']   = (ORDER+'4s',)
DAT['cdmy4']       = (ORDER+'4s',) # reserved
DAT['idmy_flag']   = (ORDER+'i',)
DAT['iscan_no']    = (ORDER+'i',)
DAT['dweather']    = (ORDER+'d', 5)
DAT['dtsys']       = (ORDER+'d',)
DAT['dant_vel']    = (ORDER+'d',)
DAT['dalc_v']      = (ORDER+'d',)
DAT['dary_scf']    = (ORDER+'d',)
DAT['dary_offset'] = (ORDER+'d',)
DAT['iline_no']    = (ORDER+'i',)
DAT['iscan_line']  = (ORDER+'i',)
DAT['iinteg_time'] = (ORDER+'i',)
DAT['idmy1']       = (ORDER+'i',) # reserved
DAT['dalpha']      = (ORDER+'d',)
DAT['dbeta']       = (ORDER+'d',)
DAT['dganma']      = (ORDER+'d',)
DAT['cdmy1']       = (ORDER+'216s',) # reserved
DAT['iary_data']   = (ORDER+'i', CHMAX)
