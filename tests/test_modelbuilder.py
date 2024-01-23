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
def test_cmems_nc_to_ini_midnight_centered(tmp_path):
    
    # TODO: create fixture
    from tests.test_interpolate_grid2bnd import cmems_dataset_4times
    ds1 = cmems_dataset_4times().isel(time=slice(None,2))
    ds1["time"] = ds1["time"] + pd.Timedelta(hours=12)
    ds2 = cmems_dataset_4times().isel(time=slice(2,None))
    ds2["time"] = ds2["time"] + pd.Timedelta(hours=12)
    
    dir_pattern = os.path.join(tmp_path, "temp_cmems_2day_*.nc")
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
                                   dir_output=tmp_path,
                                   list_quantities=["salinitybnd"],
                                   tstart="2020-01-01",
                                   dir_pattern=dir_pattern)
    
    file_expected = tmp_path / "nudge_salinity_temperature_2020-01-01_00-00-00.nc"
    
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


@pytest.mark.unittest
def test_create_model_exec_files_none(tmp_path):
    mdu_file = tmp_path / "temp_test.mdu"
    file_dimr = tmp_path / "dimr_config.xml"
    
    nproc = 1 # number of processes
    dimrset_folder = None
    mdu = hcdfm.FMModel()
    mdu.save(mdu_file)
    dfmt.create_model_exec_files(file_mdu=mdu_file, nproc=nproc, dimrset_folder=dimrset_folder)
    
    assert os.path.isfile(file_dimr)


@pytest.mark.unittest
def test_create_model_exec_files_docker(tmp_path):
    mdu_file = tmp_path / "temp_test.mdu"
    file_dimr = tmp_path / "dimr_config.xml"
    file_docker = tmp_path / "run_docker.sh"
    
    nproc = 1 # number of processes
    dimrset_folder = "docker"
    mdu = hcdfm.FMModel()
    mdu.save(mdu_file)
    dfmt.create_model_exec_files(file_mdu=mdu_file, nproc=nproc, dimrset_folder=dimrset_folder)
    
    assert os.path.isfile(file_dimr)
    assert os.path.isfile(file_docker)
    
    # check for unix line endings, this is required to be found by docker somehow
    with open(file_docker, 'rb') as f:
        data = f.readline()
    assert data == b'#!/bin/bash\n'


@pytest.mark.unittest
def test_make_paths_relative(tmp_path):
    file_pli = os.path.join(tmp_path,'test_model.pli')
    with open(file_pli,'w') as f:
        f.write("""name
                    2    2
                    1.0    2.0
                    3.0    4.0
                    """)
    
    # new ext
    ext_file_new = os.path.join(tmp_path, 'test_new.ext')
    ext_new = hcdfm.ExtModel()
    ForcingModel_object = hcdfm.ForcingModel()
    boundary_object = hcdfm.Boundary(quantity='waterlevelbnd', #the FM quantity for tide is also waterlevelbnd
                                     locationfile=file_pli,
                                     forcingfile=ForcingModel_object)
    ext_new.boundary.append(boundary_object)
    ext_new.save(filepath=ext_file_new)
        
    # old ext
    ext_file_old = os.path.join(tmp_path, 'test_old.ext')
    ext_old = hcdfm.ExtOldModel()
    forcing_old = hcdfm.ExtOldForcing(quantity='airpressure_windx_windy_charnock',
                                      filename=file_pli,
                                      filetype=hcdfm.ExtOldFileType.NetCDFGridData, #11
                                      method=hcdfm.ExtOldMethod.InterpolateTimeAndSpaceSaveWeights, #3
                                      operand=hcdfm.Operand.override) #O
    ext_old.forcing.append(forcing_old)
    ext_old.save(filepath=ext_file_old)
    
    # initialize mdu file, update settings and save
    mdu_file = os.path.join(tmp_path, 'test_model.mdu')
    mdu = hcdfm.FMModel()
    mdu.geometry.fixedweirfile = file_pli
    mdu.external_forcing.extforcefile = ext_old
    mdu.external_forcing.extforcefilenew = ext_new
    mdu.save(mdu_file)
    
    # make relative and test
    dfmt.make_paths_relative(mdu_file)
    
    # assertions: including equal sign is essential to check for relative paths
    with open(mdu_file, 'r') as file:
        filedata = file.read()
    assert "= test_model.pli" in filedata
    assert "= test_old.ext" in filedata
    assert "= test_new.ext" in filedata
    
    with open(ext_file_old, 'r') as file:
        filedata = file.read()
    assert "FILENAME=test_model.pli" in filedata
    
    with open(ext_file_new, 'r') as file:
        filedata = file.read()
    assert "locationFile = test_model.pli" in filedata
