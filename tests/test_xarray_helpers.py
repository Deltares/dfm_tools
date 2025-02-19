# -*- coding: utf-8 -*-
"""
Created on Tue Jan 23 17:12:48 2024

@author: veenstra
"""

import pytest
import dfm_tools as dfmt
from dfm_tools.xarray_helpers import file_to_list
import pandas as pd
import xarray as xr
from pathlib import Path


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
def test_file_to_list_pathlib_path():
    file_nc = dfmt.data.fm_curvedbend_his(return_filepath=True)
    file_nc_list = file_to_list(Path(file_nc))
    assert file_nc_list == [file_nc]


@pytest.mark.unittest
def test_file_to_list_empty_path():
    file_nc = "path/to/dummy/dir"
    with pytest.raises(FileNotFoundError) as e:
        _ = file_to_list(file_nc)
    assert "file(s) not found, empty file_nc_list" in str(e.value)


@pytest.mark.unittest
def test_file_to_list_already_list():
    file_nc = dfmt.data.fm_curvedbend_his(return_filepath=True)
    file_nc_list = file_to_list([file_nc])
    assert file_nc_list == [file_nc]


@pytest.mark.unittest
def test_preprocess_hisnc():
    """
    not too much added value, but good to check dropping of duplicated labels
    in this case it happens for source_sinks, not sure if this is useful here.
    More useful would be a hisfile with duplicated station names like
    p:\\archivedprojects\\11206813-006-kpp2021_rmm-2d\\C_Work\\31_RMM_FMmodel
    \\computations\\model_setup\\run_206\\results\\RMM_dflowfm_0000_his.nc
    """
    file_nc = dfmt.data.fm_grevelingen_his(return_filepath=True)
    ds1 = xr.open_dataset(file_nc)
    ds2 = xr.open_mfdataset(file_nc, preprocess=dfmt.preprocess_hisnc)
    assert ds1.sizes['source_sink'] == 46
    assert ds2.sizes['source_sink'] == 1
