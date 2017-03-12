# coding: utf-8

# Python 3.x compatibility
from __future__ import absolute_import as _absolute_import
from __future__ import division as _division
from __future__ import print_function as _print_function

# the Python Package Index
import numpy as np
import fmflow as fm
import astropy.units as u
from numpy.linalg import norm
from scipy.interpolate import interp1d
from scipy.ndimage.interpolation import map_coordinates
from scipy.special import j1

# imported items
__all__ = ['GridConverter']


class GridConverter(object):
    def __init__(self, gridsize, gdmax=None, gcf=None):
        self.gsize = gridsize
        self.gdmax = gdmax or 3.0
        self.gcf = getattr(self, (gcf or 'bessel_gauss'))
    
    def fit(self, fmarray, **kwargs):
        if fmarray.ismodulated:
            fmarray = fm.demodulate(fmarray)
        
        info = fmarray.info.copy()
        info.update(kwargs)
        
        x, y = fmarray.coord.T        
        
        try:
            x0, y0 = info['ra'], info['dec']
        except:
            x0, y0 = np.mean(x), np.mean(y)

        gs = self.gsize
        gx_rel_min = gs * np.floor(np.min(x-x0)/gs)
        gx_rel_max = gs * np.ceil(np.max(x-x0)/gs)
        gy_rel_min = gs * np.floor(np.min(y-y0)/gs)
        gy_rel_max = gs * np.ceil(np.max(y-y0)/gs)
        gx = x0 + np.arange(gx_rel_min, gx_rel_max+gs, gs)
        gy = y0 + np.arange(gy_rel_min, gy_rel_max+gs, gs)
        mgx, mgy = np.meshgrid(gx, gy)
        
        coord = fmarray.coord.copy()
        gcoord = np.vstack([mgx.flatten(), mgy.flatten()]).T
        distance = norm(coord[np.newaxis]-gcoord[:,np.newaxis], axis=2) / gs
        fmtemplate = fm.zeros_like(fmarray)
        
        self.coord, self.gcoord = coord, gcoord
        self.distance = distance
        self.fmtemplate = fmtemplate

    def transform(self, fmarray):
        if fmarray.ismodulated:
            fmarray = fm.demodulate(fmarray)        
        
        if np.any(fmarray.coord != self.coord):
            fm.utils.FMFlowError('inconsistent coordinate')
        
        array = fm.toarray(fmarray)
        
        @fm.numchunk
        def grid(distance):
            flag = (distance <= self.gdmax)
            garray = np.zeros([len(distance), array.shape[1]])

            for i in range(len(garray)):
                weight = self.gcf(distance[i][flag[i]])
                warray = weight[:,np.newaxis] * array[flag[i]]
                garray[i] = np.sum(warray, axis=0) / np.sum(weight)

            return garray
    
        garray = grid(self.distance)
        gfmarray = fm.array(garray, coord=self.gcoord, info=fmarray.info)
        return gfmarray
    
    def fit_transform(self, fmarray, **kwargs):
        self.fit(fmarray, **kwargs)
        return self.transform(fmarray)
    
    def inverse_transform(self, gfmarray):
        x, y = self.coord.T
        gx = np.unique(self.gcoord[:,0])
        gy = np.unique(self.gcoord[:,1])

        garray = fm.toarray(gfmarray)
        shape = (len(gy), len(gx), garray.shape[1])
        gcube = np.rollaxis(np.reshape(garray, shape), 2)
        
        idxs_x = interp1d(gx, np.arange(len(gx)))(x)
        idxs_y = interp1d(gy, np.arange(len(gy)))(y)
        
        @fm.numchunk
        def inverse_grid(gcube):
            array = np.zeros([len(gcube), len(self.coord)])
            
            for i in range(len(gcube)):
                array[i] = map_coordinates(gcube[i], (idxs_y, idxs_x))
            
            return array
        
        array = inverse_grid(gcube)
        fmarray = self.fmtemplate.copy()
        fmarray[:] = array.T
        
        return fmarray

    @staticmethod
    def bessel_gauss(r, a=1.55, b=2.52):
        condlist = [r != 0.0]
        choicelist = [2.0*j1(np.pi*r/a)/(np.pi*r/a) * np.exp(-(r/b)**2)]
        return np.select(condlist, choicelist, 1.0)

    @staticmethod
    def sinc_gauss(r, a=1.55, b=2.52):
        return np.sinc(r/a) * np.exp(-(r/b)**2)

    @staticmethod
    def gauss(r, a=1.00):
        return np.exp(-(r/a)**2)

