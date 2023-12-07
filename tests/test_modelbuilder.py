# -*- coding: utf-8 -*-
"""
Created on Thu Oct 26 16:08:26 2023

@author: veenstra
"""

import os
import pytest
import dfm_tools as dfmt
import pandas as pd
import hydrolib.core.dflowfm as hcdfm
import xarray as xr
from dfm_tools.hydrolib_helpers import get_ncbnd_construct


@pytest.mark.systemtest
def test_cmems_nc_to_ini_midday_centered():
    
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
    try:
        ext_old = dfmt.cmems_nc_to_ini(ext_old=ext_old,
                                   dir_output='.',
                                   list_quantities=["salinitybnd"],
                                   tstart="2020-01-01",
                                   dir_pattern=dir_pattern)
    except ValueError as e:
        assert "less than two timesteps" in str(e)
    
    # times_expected were '2019-12-31 12:00:00' and '2020-01-01 12:00:00'


@pytest.mark.systemtest
def test_cmems_nc_to_ini_midnight_centered():
    
    # TODO: create fixture
    from tests.test_interpolate_grid2bnd import cmems_dataset_4times
    ds1 = cmems_dataset_4times().isel(time=slice(None,2))
    ds1["time"] = ds1["time"] + pd.Timedelta(hours=12)
    ds2 = cmems_dataset_4times().isel(time=slice(2,None))
    ds2["time"] = ds2["time"] + pd.Timedelta(hours=12)
    
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
    
    times_expected =  ['2020-01-01 00:00:00', '2020-01-02 00:00:00']
    
    assert os.path.exists(file_expected)
    ds_out = xr.open_dataset(file_expected)
    
    times_actual = ds_out.time.to_pandas().dt.strftime("%Y-%m-%d %H:%M:%S").tolist()
    assert times_expected == times_actual
    
    assert "so" in ds_out.data_vars
    assert ds_out.so.isnull().sum().load() == 0
    
    ncbnd_construct = get_ncbnd_construct()
    varn_depth = ncbnd_construct['varn_depth']
    assert varn_depth in ds_out.coords
    # the below is inconsistent since depth is actually defined positive up, but FM requires this for inifields: https://issuetracker.deltares.nl/browse/UNST-7455
    assert ds_out[varn_depth].attrs['positive'] == 'down'
    
    # cleanup
    del ds_out
    del ds1
    del ds2
    os.remove(file_expected)
    os.remove(file_nc1)
    os.remove(file_nc2)
    os.remove(file_nc3)
    os.remove(file_nc4)


@pytest.mark.unittest
def test_create_model_exec_files():
    mdu_file = "./temp_test.mdu"
    file_dimr = "./dimr_config.xml"
    
    nproc = 1 # number of processes
    dimrset_folder = None
    mdu = hcdfm.FMModel()
    mdu.save(mdu_file)
    dfmt.create_model_exec_files(file_mdu=mdu_file, nproc=nproc, dimrset_folder=dimrset_folder)
    
    assert os.path.isfile(file_dimr)
    
    os.remove(mdu_file)
    os.remove(file_dimr)

