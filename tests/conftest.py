# -*- coding: utf-8 -*-
"""
Created on Thu May 16 15:51:44 2024

@author: veenstra
"""
import pytest
import dfm_tools as dfmt
import os


@pytest.mark.requiressecrets
@pytest.fixture
def file_nc_era5_pattern(tmp_path):
    date_min = '2010-01-31'
    date_max = '2010-02-01'
    longitude_min, longitude_max, latitude_min, latitude_max =    2,   3,  51, 52 #test domain
    variables_era5 = ['msl']#'v10n'] # check variables_dict in dfmt.download_ERA5() for valid names
    for varkey in variables_era5:
        dfmt.download_ERA5(varkey, 
                           longitude_min=longitude_min, longitude_max=longitude_max, latitude_min=latitude_min, latitude_max=latitude_max,
                           date_min=date_min, date_max=date_max,
                           dir_output=tmp_path, overwrite=True)
    
    # assert downloaded files
    file_nc_era5_pattern = os.path.join(tmp_path, "*.nc")
    return file_nc_era5_pattern
