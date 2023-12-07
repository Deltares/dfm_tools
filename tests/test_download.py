# -*- coding: utf-8 -*-
"""
Created on Thu Jun 22 09:41:49 2023

@author: veenstra
"""

import os
import pytest
import pandas as pd
from dfm_tools.download import cds_credentials, copernicusmarine_credentials
import dfm_tools as dfmt


@pytest.mark.requireslocaldata
@pytest.mark.unittest
def test_cds_credentials():
    # check if the credentials are present on this system
    val = cds_credentials()
    assert not val
    
    # some platforms depend on environ for url/apikey, check if the default cdsapi_url is set after running cds_credentials()
    assert "CDSAPI_URL" in os.environ.keys()


#TODO: properly set environment variables in github would prevent localness
@pytest.mark.requireslocaldata
@pytest.mark.unittest
def test_copernicusmarine_credentials():
    copernicusmarine_credentials()


#TODO: properly set environment variables in github would prevent localness
@pytest.mark.requireslocaldata
@pytest.mark.unittest
def test_download_era5():
    date_min = '2010-01-01'
    date_max = '2010-01-02'
    longitude_min, longitude_max, latitude_min, latitude_max =    2,   3,  51, 52 #test domain
    variables_era5 = ['msl']#'v10n'] # check variables_dict in dfmt.download_ERA5() for valid names
    dir_output = './tests/tests_output/era5_temp'
    for varkey in variables_era5:
        os.makedirs(dir_output, exist_ok=True)
        
        dfmt.download_ERA5(varkey, 
                           longitude_min=longitude_min, longitude_max=longitude_max, latitude_min=latitude_min, latitude_max=latitude_max,
                           date_min=date_min, date_max=date_max,
                           dir_output=dir_output, overwrite=True)
    # os.rmdir(dir_output)


#TODO: properly set environment variables in github would prevent localness
@pytest.mark.requireslocaldata
@pytest.mark.unittest
def test_download_cmems_my():
    date_min = '2010-01-01'
    date_max = '2010-01-02'
    longitude_min, longitude_max, latitude_min, latitude_max =    2,   3,  51, 52 #test domain
    varlist_cmems = ['bottomT','no3'] # avaliable variables differ per product, examples are ['bottomT','mlotst','siconc','sithick','so','thetao','uo','vo','usi','vsi','zos','no3']. More info on https://data.marine.copernicus.eu/products
    dir_output = './tests/tests_output/cmems_temp_my'
    for varkey in varlist_cmems:
        file_prefix = 'cmems_'
        dfmt.download_CMEMS(varkey=varkey,
                            longitude_min=longitude_min, longitude_max=longitude_max, latitude_min=latitude_min, latitude_max=latitude_max,
                            date_min=date_min, date_max=date_max,
                            dir_output=dir_output, file_prefix=file_prefix, overwrite=True)
    # os.rmdir(dir_output)


#TODO: properly set environment variables in github would prevent localness
@pytest.mark.requireslocaldata
@pytest.mark.unittest
def test_download_cmems_forecast():
    date_min = pd.Timestamp.today()
    date_max = pd.Timestamp.today() + pd.Timedelta(days=1)
    longitude_min, longitude_max, latitude_min, latitude_max =    2,   3,  51, 52 #test domain
    varlist_cmems = ['tob','no3'] # avaliable variables differ per product, examples are ['bottomT','mlotst','siconc','sithick','so','thetao','uo','vo','usi','vsi','zos','no3']. More info on https://data.marine.copernicus.eu/products
    dir_output = './tests/tests_output/cmems_temp_forecast'
    for varkey in varlist_cmems:
        file_prefix = 'cmems_'
        dfmt.download_CMEMS(varkey=varkey,
                            longitude_min=longitude_min, longitude_max=longitude_max, latitude_min=latitude_min, latitude_max=latitude_max,
                            date_min=date_min, date_max=date_max,
                            dir_output=dir_output, file_prefix=file_prefix, overwrite=True)
    # os.rmdir(dir_output)


@pytest.mark.unittest
def test_download_hycom():
    # domain
    longitude_min, longitude_max, latitude_min, latitude_max =    2,   3,  51, 52 #test domain
    date_min = '2010-01-01'
    date_max = '2010-01-02'
    varlist_hycom = ['surf_el']#'water_temp'] #['tau','water_u','water_v','water_temp','salinity','surf_el']
    
    dir_output = './tests/tests_output/hycom_temp'
    os.makedirs(dir_output, exist_ok=True)
    for varkey in varlist_hycom:
        # Path(dir_output).mkdir(parents=True, exist_ok=True)
        period_range_years = pd.period_range(date_min,date_max,freq='Y')
        dataset_url = [f'https://tds.hycom.org/thredds/dodsC/GLBu0.08/expt_19.1/{year}' for year in period_range_years] #list is possible with hycom, since it uses xr.open_mfdataset()
        file_prefix = 'hycom_'
        dfmt.download_OPeNDAP(dataset_url=dataset_url,
                              varkey=varkey,
                              longitude_min=longitude_min, longitude_max=longitude_max, latitude_min=latitude_min, latitude_max=latitude_max,
                              date_min=date_min, date_max=date_max,
                              dir_output=dir_output, file_prefix=file_prefix, overwrite=True)
