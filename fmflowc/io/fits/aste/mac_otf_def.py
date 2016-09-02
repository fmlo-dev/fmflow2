# coding: utf-8

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

# standard libraries
from collections import OrderedDict


ARYMAX = 8
CHMAX  = 1024

# record header info
HEAD = OrderedDict()
HEAD['crec_type'] = ('<4s',)
HEAD['irec_len']  = ('<i',)

# control info
CTL = OrderedDict()
CTL['cversion']  = ('<8s',)
CTL['clog_name'] = ('<60s',)
CTL['cbe_type']  = ('<8s',)

# observation info
OBS = OrderedDict()
OBS['cgroup']          = ('<8s',)
OBS['cproject']        = ('<16s',)
OBS['cobs_file']       = ('<24s',)
OBS['cobs_user']       = ('<40s',)
OBS['cst_time']        = ('<16s',)
OBS['clog_id']         = ('<24s',)
OBS['cobj_name']       = ('<16s',)
OBS['csrc_comment']    = ('<24s',)
OBS['cepoch']          = ('<8s',)
OBS['cseq_ptn']        = ('<120s',)
OBS['clagwin']         = ('<8s',)
OBS['ctrk_type']       = ('<8s',)
OBS['cscan_coord']     = ('<8s',)
OBS['cdmy1']           = ('<24s',) # reserved
OBS['icoord_type']     = ('<i',)
OBS['iscan_coord']     = ('<i',)
OBS['isw_mod']         = ('<i',)
OBS['icalib_int']      = ('<i',)
OBS['ichanel']         = ('<i',)
OBS['iary_num']        = ('<i',)
OBS['ivdef']           = ('<i',)
OBS['ivref']           = ('<i',)
OBS['imult_flag']      = ('<i',)
OBS['ich_bind']        = ('<i',)
OBS['ich_range']       = ('<i', 2)
OBS['iary_usefg']      = ('<i', ARYMAX)
OBS['iintensity_flag'] = ('<i',)
OBS['idmy1']           = ('<i', 9,) # reserved
OBS['dsrc_pos']        = ('<d', (2,3))
OBS['dsub_off']        = ('<d', 3)
OBS['dant_off']        = ('<d', 2)
OBS['dstrv']           = ('<d',)
OBS['dbsw_freq']       = ('<d',)
OBS['dbsw_ofd']        = ('<d',)
OBS['dmap_pos']        = ('<d',)
OBS['dalctime']        = ('<d',)
OBS['dmult_off']       = ('<d',)
OBS['diptim']          = ('<d',)
OBS['dsb_factor']      = ('<d', ARYMAX)
OBS['ccmt_t']          = ('<24s',)
OBS['dcmt_q']          = ('<d',)
OBS['dcmt_e']          = ('<d',)
OBS['dcmt_somega']     = ('<d',)
OBS['dcmt_node']       = ('<d',)
OBS['dcmt_i']          = ('<d',)
OBS['dant_vel']        = ('<d', ARYMAX)
OBS['dcent_freq']      = ('<d', ARYMAX)
OBS['dtrk_freq']       = ('<d', ARYMAX)
OBS['dflif']           = ('<d', ARYMAX)
OBS['dmult_scf']       = ('<d', ARYMAX)
OBS['dbebw']           = ('<d', ARYMAX)
OBS['dberes']          = ('<d', ARYMAX)
OBS['dbechwid']        = ('<d', ARYMAX)
OBS['dsw_freq']        = ('<d', ARYMAX)
OBS['dhpbw']           = ('<d', ARYMAX)
OBS['deffa']           = ('<d', ARYMAX)
OBS['deffb']           = ('<d', ARYMAX)
OBS['deffl']           = ('<d', ARYMAX)
OBS['defss']           = ('<d', ARYMAX)
OBS['dgain']           = ('<d', ARYMAX)
OBS['dpoldr']          = ('<d', ARYMAX)
OBS['ipordr']          = ('<i', ARYMAX)
OBS['iipint']          = ('<i', ARYMAX)
OBS['iary_refno']      = ('<i', ARYMAX)
OBS['imlt_no']         = ('<i', ARYMAX)
OBS['cfe_type']        = ('<10s', ARYMAX)
OBS['csid_type']       = ('<4s', ARYMAX)
OBS['chorn_typ']       = ('<2s', ARYMAX)
OBS['cpol_typ']        = ('<5s', ARYMAX)
OBS['cband_type']      = ('<2s', ARYMAX)
OBS['cposition']       = ('<32s',)
OBS['ctelescope']      = ('<16s',)
OBS['cdmy2']           = ('16s',) # reserved
OBS['ifqdat_nf']       = ('<i', ARYMAX)
OBS['dfqdat_f0']       = ('<d', ARYMAX)
OBS['dfqdat_fq']       = ('<d', (ARYMAX,10))
OBS['dfqdat_ch']       = ('<d', (ARYMAX,10))
OBS['dfqdat_cw']       = ('<d', (ARYMAX,10))

# data info
DAT = OrderedDict()
DAT['cint_sttm']   = ('<24s',)
DAT['cscan_type']  = ('<8s',)
DAT['cary_name']   = ('<4s',)
DAT['cdmy4']       = ('<4s',) # reserved
DAT['idmy_flag']   = ('<i',)
DAT['iscan_no']    = ('<i',)
DAT['dweather']    = ('<d', 5)
DAT['dtsys']       = ('<d',)
DAT['dant_vel']    = ('<d',)
DAT['dalc_v']      = ('<d',)
DAT['dary_scf']    = ('<d',)
DAT['dary_offset'] = ('<d',)
DAT['iline_no']    = ('<i',)
DAT['iscan_line']  = ('<i',)
DAT['iinteg_time'] = ('<i',)
DAT['idmy1']       = ('<i',) # reserved
DAT['dalpha']      = ('<d',)
DAT['dbeta']       = ('<d',)
DAT['dganma']      = ('<d',)
DAT['cdmy1']       = ('<216s',) # reserved
DAT['iary_data']   = ('<i', CHMAX)
