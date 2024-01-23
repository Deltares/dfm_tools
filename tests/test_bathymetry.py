# -*- coding: utf-8 -*-
"""
Created on Wed Oct 11 10:48:06 2023

@author: veenstra
"""

import os
import pytest
import numpy as np
import dfm_tools as dfmt


@pytest.mark.unittest
def test_read_write_asc(tmp_path):
    lat_range = np.arange(-8,2,0.25)
    lon_range = np.arange(-10,1,0.25)
    data = np.cos(lon_range) * lat_range[np.newaxis].T
    file_asc = os.path.join(tmp_path, 'dummy.asc')
    dfmt.write_bathy_toasc(file_asc, lon_sel_ext=lon_range, lat_sel_ext=lat_range, elev_sel_ext=data, asc_fmt='%14.9f',nodata_val=-999)
    
    ds_asc = dfmt.read_asc(file_asc)
    
    assert (np.abs(ds_asc['lon'].to_numpy() - lon_range) < 1e-9).all()
    assert (np.abs(ds_asc['lat'].to_numpy() - lat_range) < 1e-9).all()
    assert (np.abs(ds_asc['data'].to_numpy() - data) < 1e-9).all()
