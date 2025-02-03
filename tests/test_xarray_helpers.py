# -*- coding: utf-8 -*-
"""
Created on Tue Jan 23 17:12:48 2024

@author: veenstra
"""

import pytest
import dfm_tools as dfmt
import pandas as pd
import xarray as xr


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


@pytest.mark.unittest
def test_preprocess_hisnc():
    """
    not too much added value, but good to check dropping of duplicated labels
    in this case it happens for source_sinks, not sure if this is useful here
    """
    file_nc = dfmt.data.fm_grevelingen_his(return_filepath=True)
    ds1 = xr.open_dataset(file_nc)
    ds2 = xr.open_mfdataset(file_nc, preprocess=dfmt.preprocess_hisnc)
    assert ds1.sizes['source_sink'] == 46
    assert ds2.sizes['source_sink'] == 1
