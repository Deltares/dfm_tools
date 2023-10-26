# -*- coding: utf-8 -*-
"""
Created on Thu Oct 26 16:08:26 2023

@author: veenstra
"""

import os
import pytest
import dfm_tools as dfmt
import hydrolib.core.dflowfm as hcdfm
import xarray as xr


@pytest.mark.systemtest
def test_cmems_nc_to_ini():
    """
    not so covering test, preferrably include
    - conversion of quantity names
    - sal/tem combined in file, fail if no tem
    - depends on file created by other test
    """
    # TODO: depends on other test to create this
    dir_pattern = "../temp_cmems_2day_*.nc"
    
    ext_old = hcdfm.ExtOldModel()
    
    ext_old = dfmt.cmems_nc_to_ini(ext_old=ext_old,
                                   dir_output='.',
                                   list_quantities=["salinitybnd"],
                                   tstart="2020-01-01",
                                   dir_pattern=dir_pattern)

    file_expected = "./nudge_salinity_temperature_2020-01-01_00-00-00.nc"
    
    times_expected =  ['2019-12-31 12:00:00', '2020-01-01 12:00:00']
    
    assert os.path.exists(file_expected)
    ds_out = xr.open_dataset(file_expected)
    times_actual = ds_out.time.to_pandas().dt.strftime("%Y-%m-%d %H:%M:%S").tolist()
    assert "so" in ds_out.data_vars
    assert times_expected == times_actual
    ds_out.close()
    os.remove(file_expected)
    
    