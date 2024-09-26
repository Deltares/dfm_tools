# -*- coding: utf-8 -*-
"""
Created on Thu Jun 22 09:41:49 2023

@author: veenstra
"""

import os
import pytest
import pandas as pd
import cdsapi
from dfm_tools.download import (cds_credentials,
                                cds_set_credentials,
                                cds_set_credentials_rcfile,
                                cds_remove_credentials_raise,
                                copernicusmarine_credentials,
                                copernicusmarine_get_buffer,
                                )
import dfm_tools as dfmt
import xarray as xr
import glob
import numpy as np


def get_cds_url_key():
    try:
        cds_url, cds_apikey, _ = cdsapi.api.get_url_key_verify(url=None, key=None, verify=None)
    except Exception as e:
        if "Missing/incomplete configuration file" in str(e):
            cds_url = None
            cds_apikey = None
        else:
            raise e
    
    return cds_url, cds_apikey


def set_cds_credentials_ifnot_none(cds_url, cds_apikey):
    if None not in [cds_url, cds_apikey]:
        cds_set_credentials(cds_url, cds_apikey)
   

@pytest.mark.requiressecrets
@pytest.mark.unittest
def test_cds_credentials():
    cds_credentials()


@pytest.mark.requiressecrets
@pytest.mark.unittest
def test_cds_credentials_onlykey_envvars():
    cds_url, cds_apikey = get_cds_url_key()
    
    # remove credentials envvars and file
    with pytest.raises(ValueError):
        cds_remove_credentials_raise()
    
    assert "CDSAPI_URL" not in os.environ.keys()
    os.environ["CDSAPI_KEY"] = cds_apikey
    
    cds_credentials()
    set_cds_credentials_ifnot_none(cds_url, cds_apikey)


@pytest.mark.unittest
def test_cds_credentials_newurl_incorrectkey_rcfile():
    cds_url, cds_apikey = get_cds_url_key()
    
    # remove credentials envvars and file
    with pytest.raises(ValueError):
        cds_remove_credentials_raise()
    
    cds_url_temp = "https://cds-beta.climate.copernicus.eu/api"
    cds_apikey_temp = "INCORRECT-APIKEY"
    cds_set_credentials_rcfile(cds_url_temp, cds_apikey_temp)
    
    with pytest.raises(ValueError) as e:
        cds_credentials()
    set_cds_credentials_ifnot_none(cds_url, cds_apikey)
    assert "Authentication failed" in str(e.value)
    assert "The CDS/ECMWF apikey file (~/.cdsapirc) was deleted" in str(e.value)    


@pytest.mark.unittest
def test_cds_credentials_newurl_incorrectkey_envvars():
    cds_url, cds_apikey = get_cds_url_key()
    
    os.environ["CDSAPI_URL"] = "https://cds-beta.climate.copernicus.eu/api"
    os.environ["CDSAPI_KEY"] = "INCORRECT-APIKEY"
    
    with pytest.raises(ValueError) as e:
        cds_credentials()
    set_cds_credentials_ifnot_none(cds_url, cds_apikey)
    assert "Authentication failed" in str(e.value)
    assert "The CDS/ECMWF apikey file (~/.cdsapirc) was deleted" in str(e.value)


@pytest.mark.unittest
def test_cds_credentials_oldurl_incorrectkey_rcfile():
    cds_url, cds_apikey = get_cds_url_key()
    
    # remove credentials envvars and file
    with pytest.raises(ValueError):
        cds_remove_credentials_raise()
    
    cds_url_temp = "https://cds.climate.copernicus.eu/api/v2"
    cds_apikey_temp = "INCORRECT-APIKEY"
    cds_set_credentials_rcfile(cds_url_temp, cds_apikey_temp)
    
    with pytest.raises(ValueError) as e:
        cds_credentials()
    set_cds_credentials_ifnot_none(cds_url, cds_apikey)
    assert "Authentication failed" in str(e.value) # should actually be "Old CDS URL found", but the url from the file is ignored, which is acceptable
    assert "The CDS/ECMWF apikey file (~/.cdsapirc) was deleted" in str(e.value)


@pytest.mark.unittest
def test_cds_credentials_oldurl_incorrectkey_envvars():
    cds_url, cds_apikey = get_cds_url_key()
    
    os.environ["CDSAPI_URL"] = "https://cds.climate.copernicus.eu/api/v2"
    os.environ["CDSAPI_KEY"] = "INCORRECT-APIKEY"
    
    with pytest.raises(ValueError) as e:
        cds_credentials()
    set_cds_credentials_ifnot_none(cds_url, cds_apikey)
    assert "Old CDS URL found" in str(e.value)
    assert "The CDS/ECMWF apikey file (~/.cdsapirc) was deleted" in str(e.value)


@pytest.mark.unittest
def test_cds_credentials_newurl_oldkey_rcfile():
    cds_url, cds_apikey = get_cds_url_key()
    
    # remove credentials envvars and file
    with pytest.raises(ValueError):
        cds_remove_credentials_raise()
    
    cds_url_temp = "https://cds-beta.climate.copernicus.eu/api"
    cds_apikey_temp = "olduid:old-api-key"
    cds_set_credentials_rcfile(cds_url_temp, cds_apikey_temp)
    
    with pytest.raises(ValueError) as e:
        cds_credentials()
    set_cds_credentials_ifnot_none(cds_url, cds_apikey)
    assert "Old CDS API-key found (with :)" in str(e.value)
    assert "The CDS/ECMWF apikey file (~/.cdsapirc) was deleted" in str(e.value)


@pytest.mark.unittest
def test_cds_credentials_newurl_oldkey_envvars():
    cds_url, cds_apikey = get_cds_url_key()
    
    os.environ["CDSAPI_URL"] = "https://cds-beta.climate.copernicus.eu/api"
    os.environ["CDSAPI_KEY"] = "olduid:old-api-key"
    
    with pytest.raises(ValueError) as e:
        cds_credentials()
    set_cds_credentials_ifnot_none(cds_url, cds_apikey)
    assert "Old CDS API-key found (with :)" in str(e.value)
    assert "The CDS/ECMWF apikey file (~/.cdsapirc) was deleted" in str(e.value)


@pytest.mark.requiressecrets
@pytest.mark.unittest
def test_copernicusmarine_credentials():
    copernicusmarine_credentials()


@pytest.mark.requiressecrets
@pytest.mark.unittest
def test_copernicusmarine_get_buffer_notfound():
    dataset_id = 'cmems_obs-oc_glo_bgc-transp_my_l3-multi-4km_P1D'
    buffer = copernicusmarine_get_buffer(dataset_id)
    assert np.isclose(buffer, 0.5)
    

@pytest.mark.requiressecrets
@pytest.mark.unittest
def test_download_cmems_my(tmp_path):
    date_min = '2010-01-01'
    date_max = '2010-01-02'
    longitude_min, longitude_max, latitude_min, latitude_max =    2,   3,  51, 52 #test domain
    varlist_cmems = ['bottomT','no3'] # avaliable variables differ per product, examples are ['bottomT','mlotst','siconc','sithick','so','thetao','uo','vo','usi','vsi','zos','no3']. More info on https://data.marine.copernicus.eu/products
    for varkey in varlist_cmems:
        file_prefix = 'cmems_'
        dfmt.download_CMEMS(varkey=varkey,
                            longitude_min=longitude_min, longitude_max=longitude_max, latitude_min=latitude_min, latitude_max=latitude_max,
                            date_min=date_min, date_max=date_max,
                            dir_output=tmp_path, file_prefix=file_prefix, overwrite=True)
    
    # assert downloaded files
    file_nc_pat = os.path.join(tmp_path, "*.nc")
    ds = xr.open_mfdataset(file_nc_pat)
    assert ds.sizes["time"] == 2
    assert ds.time.to_pandas().iloc[0] == pd.Timestamp('2010-01-01')
    assert ds.time.to_pandas().iloc[-1] == pd.Timestamp('2010-01-02')


@pytest.mark.requiressecrets
@pytest.mark.unittest
def test_download_cmems_forecast(tmp_path):
    date_min = pd.Timestamp.today()
    date_max = pd.Timestamp.today() + pd.Timedelta(days=1)
    longitude_min, longitude_max, latitude_min, latitude_max =    2,   3,  51, 52 #test domain
    varlist_cmems = ['tob','no3'] # avaliable variables differ per product, examples are ['bottomT','mlotst','siconc','sithick','so','thetao','uo','vo','usi','vsi','zos','no3']. More info on https://data.marine.copernicus.eu/products
    for varkey in varlist_cmems:
        file_prefix = 'cmems_'
        dfmt.download_CMEMS(varkey=varkey,
                            longitude_min=longitude_min, longitude_max=longitude_max, latitude_min=latitude_min, latitude_max=latitude_max,
                            date_min=date_min, date_max=date_max,
                            dir_output=tmp_path, file_prefix=file_prefix, overwrite=True)

    # assert downloaded files
    file_nc_pat = os.path.join(tmp_path, "*.nc")
    ds = xr.open_mfdataset(file_nc_pat)
    assert ds.sizes["time"] == 3
    assert ds.time.to_pandas().iloc[0] == date_min.floor("D")
    assert ds.time.to_pandas().iloc[-1] == date_max.ceil("D")


@pytest.mark.unittest
def test_download_hycom(tmp_path):
    # domain
    longitude_min, longitude_max, latitude_min, latitude_max =    2,   3,  51, 52 #test domain
    date_min = '2010-01-01'
    date_max = '2010-01-02'
    varlist_hycom = ['surf_el']#'water_temp'] #['tau','water_u','water_v','water_temp','salinity','surf_el']
    
    for varkey in varlist_hycom:
        # Path(dir_output).mkdir(parents=True, exist_ok=True)
        period_range_years = pd.period_range(date_min,date_max,freq='Y')
        dataset_url = [f'https://tds.hycom.org/thredds/dodsC/GLBu0.08/expt_19.1/{year}' for year in period_range_years] #list is possible with hycom, since it uses xr.open_mfdataset()
        file_prefix = 'hycom_'
        dfmt.download_OPeNDAP(dataset_url=dataset_url,
                              varkey=varkey,
                              longitude_min=longitude_min, longitude_max=longitude_max, latitude_min=latitude_min, latitude_max=latitude_max,
                              date_min=date_min, date_max=date_max,
                              dir_output=tmp_path, file_prefix=file_prefix, overwrite=True)

    # assert downloaded files
    file_nc_pat = os.path.join(tmp_path, "*.nc")
    ds = xr.open_mfdataset(file_nc_pat)
    assert ds.sizes["time"] == 2
    assert ds.time.to_pandas().iloc[0] == pd.Timestamp('2010-01-01')
    assert ds.time.to_pandas().iloc[-1] == pd.Timestamp('2010-01-02')


@pytest.mark.requiressecrets
@pytest.mark.unittest
def test_download_era5_unsupported_varkey():
    date_min = '2010-01-31'
    date_max = '2010-02-01'
    longitude_min, longitude_max, latitude_min, latitude_max =    2,   3,  51, 52 #test domain
    varkey = 'unexisting'
    with pytest.raises(KeyError) as e:
        dfmt.download_ERA5(varkey, 
                            longitude_min=longitude_min, longitude_max=longitude_max, latitude_min=latitude_min, latitude_max=latitude_max,
                            date_min=date_min, date_max=date_max,
                            dir_output='.', overwrite=True)
    assert '"unexisting" not available' in str(e.value)


@pytest.mark.requiressecrets
@pytest.mark.unittest
@pytest.mark.era5slow # temporarily skip these on github
@pytest.mark.timeout(60) # useful since CDS downloads are terribly slow sometimes, so skip in that case
def test_download_era5(file_nc_era5_pattern):
    # file_nc_era5_pattern comes from conftest.py
    file_list = glob.glob(file_nc_era5_pattern)
    assert len(file_list) == 2
    
    ds = xr.open_mfdataset(file_nc_era5_pattern)
    
    assert 'valid_time' in ds.dims # TODO: if this fails, remove the exception below and in preprocess_ERA5
    
    timedim = 'time'
    # datasets retrieved with new cds-beta have valid_time instead of time dimn/varn
    # https://forum.ecmwf.int/t/new-time-format-in-era5-netcdf-files/3796/5?u=jelmer_veenstra
    # TODO: can be removed after https://github.com/Unidata/netcdf4-python/issues/1357 or https://forum.ecmwf.int/t/3796 is fixed
    if 'valid_time' in ds.dims:
        timedim = 'valid_time'
    
    assert ds.sizes[timedim] == 1416
    assert ds[timedim].to_pandas().iloc[0] == pd.Timestamp('2010-01-01')
    assert ds[timedim].to_pandas().iloc[-1] == pd.Timestamp('2010-02-28 23:00')
    
    # check if there are no integers in the dataset anymore
    # this was the case before CDS-beta in https://github.com/Deltares/dfm_tools/issues/239
    msl_encoding = ds['msl'].encoding
    assert str(msl_encoding['dtype']) == 'float32'
    assert 'scale_factor' not in msl_encoding.keys()
