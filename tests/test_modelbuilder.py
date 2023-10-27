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
    """
    # TODO: create fixture
    from tests.test_interpolate_grid2bnd import cmems_dataset_4times
    ds1 = cmems_dataset_4times().isel(time=slice(None,2))
    ds2 = cmems_dataset_4times().isel(time=slice(2,None))
    
    dir_pattern = "./temp_cmems_2day_*.nc"
    file_nc1 = dir_pattern.replace("*","sal_p1")
    file_nc2 = dir_pattern.replace("*","sal_p2")
    file_nc3 = dir_pattern.replace("*","tem_p1")
    file_nc4 = dir_pattern.replace("*","tem_p2")
    ds1.to_netcdf(file_nc1)
    ds2.to_netcdf(file_nc2)
    ds1.rename({"so":"thetao"}).to_netcdf(file_nc3)
    ds2.rename({"so":"thetao"}).to_netcdf(file_nc4)
    
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
    assert times_expected == times_actual
    assert "so" in ds_out.data_vars
    assert ds_out.so.isnull().sum().load() == 0
    
    # cleanup
    del ds_out
    del ds1
    del ds2
    os.remove(file_expected)
    os.remove(file_nc1)
    os.remove(file_nc2)
    os.remove(file_nc3)
    os.remove(file_nc4)
    