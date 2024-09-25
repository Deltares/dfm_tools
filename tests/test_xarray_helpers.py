# -*- coding: utf-8 -*-
"""
Created on Tue Jan 23 17:12:48 2024

@author: veenstra
"""

import pytest
import dfm_tools as dfmt
import pandas as pd

#TODO: many xarray_helpers tests are still in test_dfm_tools.py


@pytest.mark.unittest
@pytest.mark.requiressecrets
@pytest.mark.era5slow # temporarily skip these on github
@pytest.mark.timeout(60) # useful since CDS downloads are terribly slow sometimes, so skip in that case
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
