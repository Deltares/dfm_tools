# -*- coding: utf-8 -*-
"""
Created on Tue Jan 23 17:12:48 2024

@author: veenstra
"""


import os
import pytest
import xarray as xr
import numpy as np
from dfm_tools.xarray_helpers import prevent_dtype_int
# from .test_download import file_nc_era5_pattern

#TODO: many xarray_helpers tests are still in test_dfm_tools.py


@pytest.mark.unittest
def test_prevent_dtype_int(tmp_path, file_nc_era5_pattern):
    # file_nc_era5_pattern comes from conftest.py
    varn = "v10n"
    
    #optional encoding
    for zlib, size_expected in zip([False, True], [100e6, 50e6]):
        data_xr = xr.open_mfdataset(file_nc_era5_pattern)
        prevent_dtype_int(data_xr, zlib=zlib)
        
        #write to netcdf file
        file_out = os.path.join(tmp_path, 'era5_mslp_prevent_dtype_int.nc')
        data_xr.to_netcdf(file_out)
        data_xr_check = xr.open_dataset(file_out)
        print(data_xr_check[varn].encoding)
        
        absdiff = (data_xr_check - data_xr).apply(np.fabs)
        absdiff_max = absdiff[varn].max(dim=['longitude','latitude'])
        
        assert np.allclose(absdiff_max, 0)
    
        del data_xr
        del data_xr_check
        
        file_size = os.path.getsize(file_out)
        print(file_size)
        assert file_size < size_expected
