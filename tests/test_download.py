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
                                copernicusmarine_get_dataset_id,
                                )
import dfm_tools as dfmt
import xarray as xr
import glob
import numpy as np
from unittest.mock import patch


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


@patch("getpass.getpass")
def test_cds_credentials_prompt(getpass):
    # backup credentials and remove credentials envvars and file
    cds_url, cds_apikey = get_cds_url_key()
    with pytest.raises(ValueError):
        cds_remove_credentials_raise()
    
    # provide apikey to cds_credentials() prompt
    getpass.return_value = cds_apikey
    cds_credentials()

    # restore credentials file/envvars
    set_cds_credentials_ifnot_none(cds_url, cds_apikey)


@pytest.mark.requiressecrets
@pytest.mark.unittest
def test_cds_credentials_onlykey_envvars():
    # backup credentials and remove credentials envvars and file
    cds_url, cds_apikey = get_cds_url_key()
    with pytest.raises(ValueError):
        cds_remove_credentials_raise()
    
    # test
    assert "CDSAPI_URL" not in os.environ.keys()
    os.environ["CDSAPI_KEY"] = cds_apikey
    cds_credentials()
    
    # restore credentials file/envvars
    set_cds_credentials_ifnot_none(cds_url, cds_apikey)


@pytest.mark.unittest
def test_cds_credentials_newurl_incorrectkey_rcfile():
    # backup credentials and remove credentials envvars and file
    cds_url, cds_apikey = get_cds_url_key()
    with pytest.raises(ValueError):
        cds_remove_credentials_raise()
    
    # test
    cds_url_temp = "https://cds.climate.copernicus.eu/api"
    cds_apikey_temp = "INCORRECT-APIKEY"
    cds_set_credentials_rcfile(cds_url_temp, cds_apikey_temp)
    with pytest.raises(ValueError) as e:
        cds_credentials()
    assert "Authentication failed" in str(e.value)
    assert "The CDS/ECMWF apikey file (~/.cdsapirc) was deleted" in str(e.value)    
    
    # restore credentials file/envvars
    set_cds_credentials_ifnot_none(cds_url, cds_apikey)


@pytest.mark.unittest
def test_cds_credentials_newurl_incorrectkey_envvars():
    # backup credentials and remove credentials envvars and file
    cds_url, cds_apikey = get_cds_url_key()
    with pytest.raises(ValueError):
        cds_remove_credentials_raise()
    
    # test
    os.environ["CDSAPI_URL"] = "https://cds.climate.copernicus.eu/api"
    os.environ["CDSAPI_KEY"] = "INCORRECT-APIKEY"
    with pytest.raises(ValueError) as e:
        cds_credentials()
    assert "Authentication failed" in str(e.value)
    assert "The CDS/ECMWF apikey file (~/.cdsapirc) was deleted" in str(e.value)
    
    # restore credentials file/envvars
    set_cds_credentials_ifnot_none(cds_url, cds_apikey)


@pytest.mark.unittest
def test_cds_credentials_oldurl_incorrectkey_rcfile():
    # backup credentials and remove credentials envvars and file
    cds_url, cds_apikey = get_cds_url_key()
    with pytest.raises(ValueError):
        cds_remove_credentials_raise()
    
    # test
    cds_url_temp = "https://cds.climate.copernicus.eu/api/v2"
    cds_apikey_temp = "INCORRECT-APIKEY"
    cds_set_credentials_rcfile(cds_url_temp, cds_apikey_temp)
    with pytest.raises(ValueError) as e:
        with pytest.warns(UserWarning) as w:
            cds_credentials()
    assert "Old CDS URL found" in str(e.value)
    assert "The CDS/ECMWF apikey file (~/.cdsapirc) was deleted" in str(e.value)
    assert "404 Client Error: Not Found for url:" in str(w[0].message)
    
    # test
    cds_url_temp = "https://cds-beta.climate.copernicus.eu/api"
    cds_apikey_temp = "INCORRECT-APIKEY"
    cds_set_credentials_rcfile(cds_url_temp, cds_apikey_temp)
    with pytest.raises(ValueError) as e:
        with pytest.warns(UserWarning) as w:
            cds_credentials()
    assert "Old CDS URL found" in str(e.value)
    assert "The CDS/ECMWF apikey file (~/.cdsapirc) was deleted" in str(e.value)
    assert "certificate verify failed" in str(w[0].message)
    
    # restore credentials file/envvars
    set_cds_credentials_ifnot_none(cds_url, cds_apikey)


@pytest.mark.unittest
def test_cds_credentials_oldurl_incorrectkey_envvars():
    # backup credentials and remove credentials envvars and file
    cds_url, cds_apikey = get_cds_url_key()
    with pytest.raises(ValueError):
        cds_remove_credentials_raise()
    
    # test
    os.environ["CDSAPI_URL"] = "https://cds.climate.copernicus.eu/api/v2"
    os.environ["CDSAPI_KEY"] = "INCORRECT-APIKEY"
    with pytest.raises(ValueError) as e:
        with pytest.warns(UserWarning) as w:
            cds_credentials()
    assert "Old CDS URL found" in str(e.value)
    assert "The CDS/ECMWF apikey file (~/.cdsapirc) was deleted" in str(e.value)
    assert "404 Client Error: Not Found for url:" in str(w[0].message)
    
    # test
    os.environ["CDSAPI_URL"] = "https://cds-beta.climate.copernicus.eu/api"
    os.environ["CDSAPI_KEY"] = "INCORRECT-APIKEY"
    with pytest.raises(ValueError) as e:
        with pytest.warns(UserWarning) as w:
            cds_credentials()
    assert "Old CDS URL found" in str(e.value)
    assert "The CDS/ECMWF apikey file (~/.cdsapirc) was deleted" in str(e.value)
    assert "certificate verify failed" in str(w[0].message)

    # restore credentials file/envvars
    set_cds_credentials_ifnot_none(cds_url, cds_apikey)


@pytest.mark.unittest
def test_cds_credentials_newurl_oldkey_rcfile():
    # backup credentials and remove credentials envvars and file
    cds_url, cds_apikey = get_cds_url_key()
    with pytest.raises(ValueError):
        cds_remove_credentials_raise()
    
    # test
    cds_url_temp = "https://cds.climate.copernicus.eu/api"
    cds_apikey_temp = "olduid:old-api-key"
    cds_set_credentials_rcfile(cds_url_temp, cds_apikey_temp)
    with pytest.raises(ValueError) as e:
        cds_credentials()
    assert "Old CDS API-key found (with :)" in str(e.value)
    assert "The CDS/ECMWF apikey file (~/.cdsapirc) was deleted" in str(e.value)
    
    # restore credentials file/envvars
    set_cds_credentials_ifnot_none(cds_url, cds_apikey)


@pytest.mark.unittest
def test_cds_credentials_newurl_oldkey_envvars():
    # backup credentials and remove credentials envvars and file
    cds_url, cds_apikey = get_cds_url_key()
    with pytest.raises(ValueError):
        cds_remove_credentials_raise()
    
    # test
    os.environ["CDSAPI_URL"] = "https://cds.climate.copernicus.eu/api"
    os.environ["CDSAPI_KEY"] = "olduid:old-api-key"
    with pytest.raises(ValueError) as e:
        cds_credentials()
    assert "Old CDS API-key found (with :)" in str(e.value)
    assert "The CDS/ECMWF apikey file (~/.cdsapirc) was deleted" in str(e.value)
    
    # restore credentials file/envvars
    set_cds_credentials_ifnot_none(cds_url, cds_apikey)


@pytest.mark.requiressecrets
@pytest.mark.unittest
def test_copernicusmarine_credentials():
    copernicusmarine_credentials()


@pytest.mark.requiressecrets
@pytest.mark.unittest
def test_copernicusmarine_get_dataset_id():
    # with string datetimes, these are parsed to pandas timestamps in
    # copernicusmarine_get_product()
    date_min = '2010-01-01'
    date_max = '2010-01-02'
    date_args = dict(date_min=date_min, date_max=date_max)
    dataset_id = copernicusmarine_get_dataset_id(varkey='bottomT', **date_args)
    assert dataset_id == 'cmems_mod_glo_phy_my_0.083deg_P1D-m'
    dataset_id = copernicusmarine_get_dataset_id(varkey='no3', **date_args)
    assert dataset_id == 'cmems_mod_glo_bgc_my_0.25deg_P1D-m'
    dataset_id = copernicusmarine_get_dataset_id(varkey='phyc', **date_args)
    assert dataset_id == 'cmems_mod_glo_bgc_my_0.25deg_P1M-m'
    
    # with pandas timestamps
    date_min = pd.Timestamp.today()
    date_max = pd.Timestamp.today() + pd.Timedelta(days=1)
    date_args = dict(date_min=date_min, date_max=date_max)
    dataset_id = copernicusmarine_get_dataset_id(varkey='tob', **date_args)
    assert dataset_id == 'cmems_mod_glo_phy_anfc_0.083deg_P1D-m'
    dataset_id = copernicusmarine_get_dataset_id(varkey='no3', **date_args)
    assert dataset_id == 'cmems_mod_glo_bgc-nut_anfc_0.25deg_P1D-m'


@pytest.mark.parametrize("varkey", [pytest.param(x, id=x) for x in ['bottomT','no3','so','tob']])
@pytest.mark.requiressecrets
@pytest.mark.unittest
def test_download_cmems(tmp_path, varkey):
    """
    Test whether downloading cmems data works properly.
    Especially whether time/spatial extents are correct (buffered)
    and whether daily means are corrected with a 12 hour offset.
    
    the variables retrieved are arbitrary, but are present in the respective datasets
    avaliable variables differ per product, examples are ['bottomT','mlotst','siconc','sithick','so','thetao','uo','vo','usi','vsi','zos','no3'].
    More info on https://data.marine.copernicus.eu/products
    
    In the variable/dataset selection it is ensured to cover both
    reanalysis and analysisforecasts, both daily and monthly means,
    both bio and phy, both 0.25deg and 0.083deg spatial resolutions
    and both variables including and excluding depth
    """
    
    # deliberately take inconvenient time/spatial subset to test if
    # coordinates_selection_method='outside' works properly
    longitude_min, longitude_max, latitude_min, latitude_max =    2.001,   3.001,  51.001, 52.001 #test domain
    dataset_id_dict = {'bottomT':'cmems_mod_glo_phy_my_0.083deg_P1D-m', # phy my daily mean no_depth
                       'no3':'cmems_mod_glo_bgc_my_0.25deg_P1D-m', # bio my daily mean wi_depth
                       'so':'cmems_mod_glo_phy_my_0.083deg_P1M-m', # phy my monthly mean wi_depth
                       'tob':'cmems_mod_glo_phy_anfc_0.083deg_P1D-m', # phy anfc daily mean no_depth
                       }
    file_prefix = 'cmems_'
    dataset_id = dataset_id_dict[varkey]
    # freq is converted from "D" to "M" automatically in for monthly datasets
    freq = "D"

    if "_my_" in dataset_id:
        date_min = '2010-01-01 01:00'
        date_max = '2010-01-01 23:00'
    elif "_anfc_"in dataset_id:
        date_today = pd.Timestamp.today().floor("1D")
        date_min = date_today + pd.Timedelta(hours=1)
        date_max = date_today + pd.Timedelta(hours=23)
    
    if varkey == 'so':
        # monthly mean, times are at midnight (no offset applied)
        times_expected = ['2010-01-01 00:00:00', '2010-02-01 00:00:00']
        fnames_expected = ['cmems_so_2010-01.nc', 'cmems_so_2010-02.nc']
    elif varkey == 'no3':
        # daily mean, so offset of 12 hours is applied
        times_expected = ['2009-12-31 12:00:00',
                          '2010-01-01 12:00:00',
                          '2010-01-02 12:00:00',
                          ]
        fnames_expected = ['cmems_no3_2009-12-31.nc',
                           'cmems_no3_2010-01-01.nc',
                           'cmems_no3_2010-01-02.nc',
                           ]
    elif varkey == 'bottomT':
        # daily mean, so offset of 12 hours is applied
        times_expected = ['2009-12-31 12:00:00',
                          '2010-01-01 12:00:00',
                          '2010-01-02 12:00:00']
        fnames_expected = ['cmems_bottomT_2009-12-31.nc',
                           'cmems_bottomT_2010-01-01.nc',
                           'cmems_bottomT_2010-01-02.nc',
                           ]
    elif varkey == 'tob':
        # daily mean, so offset of 12 hours is applied
        datetime_first = date_today - pd.Timedelta(hours=12)
        date_first = datetime_first.date()
        times_expected = [str(datetime_first + x*pd.Timedelta("1D")) for x in [0,1,2]]
        fnames_dates = [str(date_first + x*pd.Timedelta("1D")) for x in [0,1,2]]
        fnames_expected = [f"cmems_tob_{x}.nc" for x in fnames_dates]
    
    if '0.25deg' in dataset_id:
        lon_max_exp = 3.25
        lat_max_exp = 52.25
    elif '0.083deg' in dataset_id:
        lon_max_exp = 3.08333
        lat_max_exp = 52.083332
        
    dfmt.download_CMEMS(varkey=varkey,
                        longitude_min=longitude_min, longitude_max=longitude_max, latitude_min=latitude_min, latitude_max=latitude_max,
                        date_min=date_min, date_max=date_max, freq=freq,
                        # speed up tests by supplying datset_id
                        dataset_id=dataset_id,
                        dir_output=tmp_path, file_prefix=file_prefix, overwrite=True)

    # open downloaded files
    file_nc_pat = os.path.join(tmp_path, f"{file_prefix}{varkey}*.nc")
    fname_list = [os.path.basename(x) for x in glob.glob(file_nc_pat)]
    fname_list.sort()
    ds = xr.open_mfdataset(file_nc_pat)
    times_actual = ds.time.to_pandas().dt.strftime("%Y-%m-%d %H:%M:%S").tolist()
    
    # assertions
    assert varkey in set(ds.variables)
    assert ds.sizes["time"] == len(times_expected)
    assert fname_list == fnames_expected
    assert times_actual == times_expected
    assert np.isclose(ds.longitude.to_numpy().min(), 2)
    assert np.isclose(ds.longitude.to_numpy().max(), lon_max_exp)
    assert np.isclose(ds.latitude.to_numpy().min(), 51)
    assert np.isclose(ds.latitude.to_numpy().max(), lat_max_exp)


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
        # temporary fix to avoid RuntimeError: NetCDF: file not found
        # https://github.com/Deltares/dfm_tools/issues/1048
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
    longitude_min, longitude_max, latitude_min, latitude_max = 2, 3, 51, 52
    varkey = 'unexisting'
    with pytest.raises(KeyError) as e:
        dfmt.download_ERA5(
            varkey, 
            longitude_min=longitude_min, longitude_max=longitude_max,
            latitude_min=latitude_min, latitude_max=latitude_max,
            date_min=date_min, date_max=date_max,
            dir_output='.', overwrite=True)
    assert '"unexisting" not available' in str(e.value)


@pytest.mark.requiressecrets
@pytest.mark.unittest
def test_download_era5_incorrect_times():
    date_min = '2010-01-01'
    date_max = '2009-01-01'
    longitude_min, longitude_max, latitude_min, latitude_max = 2, 3, 51, 52
    varkey = 'msl'
    with pytest.raises(ValueError) as e:
        dfmt.download_ERA5(
            varkey, 
            longitude_min=longitude_min, longitude_max=longitude_max,
            latitude_min=latitude_min, latitude_max=latitude_max,
            date_min=date_min, date_max=date_max,
            dir_output='.', overwrite=True)
    assert 'resulted in empty period_range' in str(e.value)


@pytest.mark.requiressecrets
@pytest.mark.unittest
@pytest.mark.era5slow # temporarily skip these on github
@pytest.mark.timeout(60) # useful since CDS downloads are terribly slow sometimes, so skip in that case
def test_download_era5(file_nc_era5_pattern):
    # file_nc_era5_pattern comes from conftest.py
    file_list = glob.glob(file_nc_era5_pattern)
    assert len(file_list) == 2
    
    ds = xr.open_mfdataset(file_nc_era5_pattern)
    
    # datasets retrieved with intermediate CDS had expver dimension causing issues
    # if this fails, update the docstring of preprocess_ERA5
    assert 'expver' not in ds.dims
    # if this fails, update preprocess_ERA5 to not add coordinate varialbe
    assert 'expver' in ds.coords
    
    # datasets retrieved with new CDS have valid_time instead of time dim/var
    assert 'valid_time' in ds.dims # TODO: if this fails, remove the exception below and in preprocess_ERA5
    timedim = 'valid_time'
    
    assert ds.sizes[timedim] == 1416
    assert ds[timedim].to_pandas().iloc[0] == pd.Timestamp('2010-01-01')
    assert ds[timedim].to_pandas().iloc[-1] == pd.Timestamp('2010-02-28 23:00')
    
    # datasets retrieved with new CDS are float32 instead of scaled ints
    # ints raised problem in https://github.com/Deltares/dfm_tools/issues/239
    msl_encoding = ds['msl'].encoding
    assert str(msl_encoding['dtype']) == 'float32'
    assert 'scale_factor' not in msl_encoding.keys()
