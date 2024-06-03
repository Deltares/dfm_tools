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
import dfm_tools as dfmt
import pandas as pd

#TODO: many xarray_helpers tests are still in test_dfm_tools.py


@pytest.mark.requiressecrets
@pytest.mark.unittest
@pytest.mark.timeout(60) # useful since CDS downloads are terribly slow sometimes, so skip in that case
def test_prevent_dtype_int(tmp_path, file_nc_era5_pattern):
    # file_nc_era5_pattern comes from file_nc_era5_pattern() in conftest.py
    varn = "msl"
    file_nc = os.path.join(tmp_path,f'era5_{varn}_*.nc')
    
    #optional encoding
    for zlib, size_expected in zip([False, True], [480000, 250000]):
        data_xr = xr.open_mfdataset(file_nc)
        prevent_dtype_int(data_xr, zlib=zlib)
        
        #write to netcdf file
        file_out = os.path.join(tmp_path, f'era5_prevent_dtype_int_{varn}_zlib_{zlib}.nc')
        data_xr.to_netcdf(file_out)
        data_xr_check = xr.open_dataset(file_out)
        print(data_xr_check[varn].encoding)
        
        assert np.allclose(data_xr[varn], data_xr_check[varn])
        
        del data_xr
        del data_xr_check
        
        file_size = os.path.getsize(file_out)
        print(file_size)
        assert file_size < size_expected


@pytest.mark.requiressecrets
@pytest.mark.unittest
def test_merge_meteofiles(file_nc_era5_pattern):
    # file_nc_era5_pattern comes from file_nc_era5_pattern() in conftest.py
    ds = dfmt.merge_meteofiles(file_nc=file_nc_era5_pattern, 
                               preprocess=dfmt.preprocess_ERA5, 
                               time_slice=slice("2010-01-30","2010-02-01")
                               )
    assert ds.sizes["time"] == 72
    assert ds.time.to_pandas().iloc[0] == pd.Timestamp('2010-01-30')
    assert ds.time.to_pandas().iloc[-1] == pd.Timestamp('2010-02-01 23:00')
    assert "msl" in ds.data_vars
