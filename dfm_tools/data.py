# -*- coding: utf-8 -*-
"""
Created on Sat Jul  8 09:38:01 2023

@author: veenstra
"""

import os
import requests
import xarray as xr
import xugrid as xu
import pooch
import zipfile
from dfm_tools import open_partitioned_dataset, preprocess_hisnc, open_dataset_delft3d4


def get_dir_testdata():
    # create cache dir like %USERPROFILE%/AppData/Local/dfm_tools/dfm_tools/Cache
    # TODO: add SHA256 checking and more, like: https://github.com/Deltares/xugrid/blob/main/xugrid/data/sample_data.py
    dir_testdata = str(pooch.os_cache('dfm_tools'))
    os.makedirs(dir_testdata, exist_ok=True)
    return dir_testdata


def maybe_download_opendap_data(file_nc,dir_subfolder=None):
    if os.path.exists(file_nc): # only download if file does not exist already
        return
    
    #opendap catalog: https://opendap.deltares.nl/thredds/catalog/opendap/deltares/Delft3D/netcdf_example_files/catalog.html
    opendap_url = 'https://opendap.deltares.nl/thredds/fileServer/opendap/deltares/Delft3D/netcdf_example_files'
    if dir_subfolder is not None:
        opendap_url = f'{opendap_url}/{dir_subfolder}'
    fname = os.path.basename(file_nc)
    file_url = f'{opendap_url}/{fname}'

    print(f'downloading "{fname}" from opendap.deltares.nl to cachedir')
    r = requests.get(file_url, allow_redirects=True)
    r.raise_for_status() #raise HTTPError if url not exists
    with open(file_nc, 'wb') as f:
        f.write(r.content)


def fm_grevelingen_map(return_filepath:bool = False) -> xu.UgridDataset:
    
    dir_subfolder = 'DFM_grevelingen_3D'
    dir_testdata = get_dir_testdata()
    
    file_nc_pat = os.path.join(dir_testdata,'Grevelingen-FM_0*_map.nc')
    
    #download data if not present
    for part in range(8):
        file_nc = file_nc_pat.replace('0*',f'{part:04d}')
        maybe_download_opendap_data(file_nc,dir_subfolder)
    
    #potentially only return filepath of downloaded file(s)
    filepath = file_nc_pat
    if return_filepath:
        return filepath
    
    #open as xugrid.UgridDataset
    uds = open_partitioned_dataset(filepath)
    return uds


def fm_grevelingen_his(return_filepath:bool = False) -> xr.Dataset:
    
    dir_subfolder = 'DFM_grevelingen_3D'
    dir_testdata = get_dir_testdata()
    
    #download data if not present
    file_nc = os.path.join(dir_testdata,'Grevelingen-FM_0000_his.nc')
    maybe_download_opendap_data(file_nc,dir_subfolder)
    
    #potentially only return filepath of downloaded file(s)
    filepath = file_nc
    if return_filepath:
        return filepath
    
    #open as xarray.Dataset
    ds = xr.open_mfdataset(filepath,preprocess=preprocess_hisnc)
    return ds


def fm_grevelingen_net(return_filepath:bool = False) -> xu.UgridDataset:
    
    dir_subfolder = 'DFM_grevelingen_3D'
    dir_testdata = get_dir_testdata()
    
    #download data if not present
    file_nc = os.path.join(dir_testdata,'Grevelingen_FM_grid_20190603_net.nc')
    maybe_download_opendap_data(file_nc,dir_subfolder)
    
    #potentially only return filepath of downloaded file(s)
    filepath = file_nc
    if return_filepath:
        return filepath
    
    #open as xugrid.UgridDataset
    uds = open_partitioned_dataset(filepath)
    return uds


def fm_curvedbend_map(return_filepath:bool = False) -> xu.UgridDataset:
    
    dir_subfolder = 'DFM_curvedbend_3D'
    dir_testdata = get_dir_testdata()
    
    #download data if not present
    file_nc = os.path.join(dir_testdata,'cb_3d_map.nc')
    maybe_download_opendap_data(file_nc,dir_subfolder)
    
    #potentially only return filepath of downloaded file(s)
    filepath = file_nc
    if return_filepath:
        return filepath
    
    #open as xugrid.UgridDataset
    uds = open_partitioned_dataset(filepath)
    return uds


def fm_curvedbend_his(return_filepath:bool = False) -> xr.Dataset:
    
    dir_subfolder = 'DFM_curvedbend_3D'
    dir_testdata = get_dir_testdata()
    
    #download data if not present
    file_nc = os.path.join(dir_testdata,'cb_3d_his.nc')
    maybe_download_opendap_data(file_nc,dir_subfolder)
    
    #potentially only return filepath of downloaded file(s)
    filepath = file_nc
    if return_filepath:
        return filepath
    
    #open as xarray.Dataset
    ds = xr.open_mfdataset(filepath,preprocess=preprocess_hisnc)
    return ds


def fm_westernscheldt_map(return_filepath:bool = False) -> xu.UgridDataset:
    
    dir_testdata = get_dir_testdata()
    
    #download data if not present
    file_nc = os.path.join(dir_testdata,'westernscheldt_sph_map.nc')
    maybe_download_opendap_data(file_nc)
    
    #potentially only return filepath of downloaded file(s)
    filepath = file_nc
    if return_filepath:
        return filepath
    
    #open as UgridDataset
    uds = open_partitioned_dataset(filepath)
    return uds


def d3d_westernscheldt_trim(return_filepath:bool = False) -> xu.UgridDataset:
    
    dir_testdata = get_dir_testdata()
    
    #download data if not present
    file_nc = os.path.join(dir_testdata,'trim-westernscheldt_sph.nc')
    maybe_download_opendap_data(file_nc)
    
    #potentially only return filepath of downloaded file(s)
    filepath = file_nc
    if return_filepath:
        return filepath
    
    #open as UgridDataset
    uds = open_dataset_delft3d4(filepath)
    return uds


def d3d_curvedbend_trim(return_filepath:bool = False) -> xu.UgridDataset:
    
    dir_testdata = get_dir_testdata()
    
    #download data if not present
    file_nc = os.path.join(dir_testdata,'trim-cb2-sal-added-3d.nc')
    maybe_download_opendap_data(file_nc)
    
    #potentially only return filepath of downloaded file(s)
    filepath = file_nc
    if return_filepath:
        return filepath
    
    #open as UgridDataset
    uds = open_dataset_delft3d4(filepath)
    return uds


def d3d_curvedbend_trih(return_filepath:bool = False) -> xr.Dataset:
    
    dir_testdata = get_dir_testdata()
    
    #download data if not present
    file_nc = os.path.join(dir_testdata,'trih-cb2-sal-added-3d.nc')
    maybe_download_opendap_data(file_nc)
    
    #potentially only return filepath of downloaded file(s)
    filepath = file_nc
    if return_filepath:
        return filepath
    
    #open as UgridDataset
    ds = xr.open_mfdataset(filepath,preprocess=preprocess_hisnc)
    return ds


def gshhs_coastlines_shp() -> str:
    """
    Downloads and unzips GSHHS coastlines from NOAA if not present.
    Checks presence of files for all five resolutions.

    Returns
    -------
    str
        DESCRIPTION.

    """
    
    dir_testdata = get_dir_testdata()
    
    fname = 'gshhg-shp-2.3.7.zip'
    filepath_zip = os.path.join(dir_testdata,fname)
    dir_gshhs = os.path.join(dir_testdata,'gshhg-shp-2.3.7')
    
    #download zipfile if not present
    if not os.path.exists(filepath_zip) and not os.path.exists(dir_gshhs):
        file_url = f'https://www.ngdc.noaa.gov/mgg/shorelines/data/gshhg/latest/{fname}'
        print(f'downloading "{fname}" from www.ngdc.noaa.gov to cachedir')
        r = requests.get(file_url, allow_redirects=True)
        r.raise_for_status() #raise HTTPError if url not exists
        with open(filepath_zip, 'wb') as f:
            f.write(r.content)

    #unzip zipfile if unzipped folder not present
    if not os.path.exists(dir_gshhs):
        print(f'unzipping "{fname}"')
        with zipfile.ZipFile(filepath_zip, 'r') as zip_ref:
            zip_ref.extractall(dir_gshhs)
    
    #construct filepath list and check existence of shapefiles
    filepath_shp_list = [os.path.join(dir_gshhs,'GSHHS_shp',res,f'GSHHS_{res}_L1.shp') for res in ['f','h','i','l','c']]
    for filepath_shp in filepath_shp_list:
        assert os.path.exists(filepath_shp) #coastlines
        assert os.path.exists(filepath_shp.replace('L1.shp','L2.shp')) #lakes
        assert os.path.exists(filepath_shp.replace('L1.shp','L3.shp')) #islands-in-lakes
        assert os.path.exists(filepath_shp.replace('L1.shp','L6.shp')) #Antarctic grounding-line polygons
    
    return dir_gshhs

