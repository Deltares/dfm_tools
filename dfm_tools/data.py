# -*- coding: utf-8 -*-
"""
Created on Sat Jul  8 09:38:01 2023

@author: veenstra
"""

import os
import requests
import xarray as xr
import xugrid as xu
import zipfile
import geopandas as gpd
from dfm_tools import open_partitioned_dataset, preprocess_hisnc, open_dataset_delft3d4
# TODO: work with pooch instead, like: https://github.com/Deltares/xugrid/blob/main/xugrid/data/sample_data.py


def get_dir_testdata(dir_subfolder=''):
    dir_testdata = os.path.join(r'c:\DATA','dfm_tools_data', dir_subfolder) #on WCF
    if os.path.exists(dir_testdata):
        return dir_testdata
    dir_testdata = './dfm_tools_data' #for instance when running on github/binder
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

    print(f'downloading "{fname}" from OPeNDAP to "{os.path.dirname(file_nc)}"')
    r = requests.get(file_url, allow_redirects=True)
    r.raise_for_status() #raise HTTPError if url not exists
    with open(file_nc, 'wb') as f:
        f.write(r.content)


def fm_grevelingen_map(return_filepath:bool = False) -> xu.UgridDataset:
    
    dir_subfolder = 'DFM_grevelingen_3D'
    dir_testdata = get_dir_testdata(dir_subfolder)
    
    file_nc_pat = os.path.join(dir_testdata,'Grevelingen-FM_0*_map.nc')
    
    #download data if not present
    file_nc_list = [file_nc_pat.replace('0*',f'{i:04d}') for i in range(8)]
    for file_nc in file_nc_list:
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
    dir_testdata = get_dir_testdata(dir_subfolder)
    
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
    dir_testdata = get_dir_testdata(dir_subfolder)
    
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
    dir_testdata = get_dir_testdata(dir_subfolder)
    
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
    dir_testdata = get_dir_testdata(dir_subfolder)
    
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


def gshhs_coastlines_shp_download(return_filepath=True) -> str:
    
    dir_testdata = get_dir_testdata()
    
    fname = 'gshhg-shp-2.3.7.zip'
    filepath_zip = os.path.join(dir_testdata,fname)
    filepath_dir = os.path.join(dir_testdata,'gshhg-shp-2.3.7')
    
    #download zipfile if not present
    if not os.path.exists(filepath_zip) and not os.path.exists(filepath_dir):
        file_url = f'https://www.ngdc.noaa.gov/mgg/shorelines/data/gshhg/latest/{fname}'
        print(f'downloading "{fname}" from www.ngdc.noaa.gov to "{os.path.basename(filepath_zip)}"')
        r = requests.get(file_url, allow_redirects=True)
        r.raise_for_status() #raise HTTPError if url not exists
        with open(filepath_zip, 'wb') as f:
            f.write(r.content)

    #unzip zipfile if unzipped folder not present
    if not os.path.exists(filepath_dir):
        print(f'unzipping "{fname}"')
        with zipfile.ZipFile(filepath_zip, 'r') as zip_ref:
            zip_ref.extractall(filepath_dir)
    
    #construct filepath list and check existence of shapefiles
    filepaths_shp_dict = {res:os.path.join(filepath_dir,'GSHHS_shp',res,f'GSHHS_{res}_L1.shp') for res in ['f','h','i','l','c']}
    for filepath_shp in filepaths_shp_dict.values():
        assert os.path.exists(filepath_shp) #coastlines
        assert os.path.exists(filepath_shp.replace('L1.shp','L2.shp')) #lakes
        assert os.path.exists(filepath_shp.replace('L1.shp','L3.shp')) #islands-in-lakes
        assert os.path.exists(filepath_shp.replace('L1.shp','L6.shp')) #Antarctic grounding-line polygons
    
    return filepaths_shp_dict
    
    
    

if __name__ == "__main__":
    func_list = [gshhs_coastlines_shp_download,
                 #fm_grevelingen_map, fm_grevelingen_his, fm_grevelingen_net, 
                 #fm_curvedbend_map, fm_curvedbend_his, 
                 #fm_westernscheldt_map, 
                 #d3d_westernscheldt_trim
                 ]
    for func in func_list:
        file_nc = func(return_filepath=True)
        print(file_nc)
        uds = func()


