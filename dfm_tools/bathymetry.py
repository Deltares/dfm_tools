# -*- coding: utf-8 -*-
"""
Created on Mon May 31 17:17:09 2021

@author: veenstra
"""

import numpy as np


def write_bathy_toasc(filename_asc,lon_sel_ext,lat_sel_ext,elev_sel_ext,asc_fmt='%9.2f',nodata_val=32767):
    """
    empty docstring
    """

    print('writing to asc file')
    if elev_sel_ext.shape != (lat_sel_ext.shape[0], lon_sel_ext.shape[0]):
        raise ValueError('resulting shape of elev_sel_ext %s is not consistent with lat_sel_ext/lon_sel_ext vars %s'%(elev_sel_ext.shape,(lat_sel_ext.shape[0], lon_sel_ext.shape[0])))
    if isinstance(elev_sel_ext,np.ma.core.MaskedArray): #masked has to be filled in order for the nans to be visible
        elev_sel_ext = elev_sel_ext.filled(np.nan)
    if np.isnan(elev_sel_ext).sum()>0:
        elev_sel_ext = elev_sel_ext.copy()
        elev_sel_ext[np.isnan(elev_sel_ext)] = nodata_val
        print('replaced nan values with %d'%(nodata_val))
    resinv_lonlat = np.round(1/np.diff(lon_sel_ext[:2]),2)
    resinv_lat = np.round(1/np.diff(lat_sel_ext[:2]),2)
    if resinv_lonlat!=resinv_lat:
        raise ValueError('inconsistent resolution over lat/lon')
    
    with open(filename_asc,'w') as file_asc:
        file_asc.write('ncols         %d\n'%(elev_sel_ext.shape[1]))
        file_asc.write('nrows         %d\n'%(elev_sel_ext.shape[0]))
        file_asc.write('xllcenter     %13.8f\n'%(lon_sel_ext[0]))
        file_asc.write('yllcenter     %13.8f\n'%(lat_sel_ext[0]))
        file_asc.write('cellsize      %16.11f\n'%(1/resinv_lonlat))
        #file_asc.write('NODATA_value  %d\n'%(nodata_val))
        file_asc.write('NODATA_value  '+asc_fmt%(nodata_val)+'\n')
    with open(filename_asc,'a') as file_asc:
        np.savetxt(file_asc,np.flip(elev_sel_ext,axis=0),fmt=asc_fmt)
    print('...finished')

