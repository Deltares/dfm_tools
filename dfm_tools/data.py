# -*- coding: utf-8 -*-
"""
Created on Sat Jul  8 09:38:01 2023

@author: veenstra
"""

import os
import requests
import xarray as xr
import xugrid as xu
from dfm_tools import open_partitioned_dataset, preprocess_hisnc, open_dataset_delft3d4
# TODO: work with pooch instead, like: https://github.com/Deltares/xugrid/blob/main/xugrid/data/sample_data.py


def get_dir_testdata(dir_subfolder=''):
    dir_testdata = os.path.join(r'c:\DATA','dfm_tools_testdata', dir_subfolder) #on WCF
    if os.path.exists(dir_testdata):
        return dir_testdata
    dir_testdata = './data' #for instance when running on github/binder
    os.makedirs(dir_testdata, exist_ok=True)
    return dir_testdata


def maybe_download_opendap_data(file_nc,dir_subfolder=None):
    if os.path.exists(file_nc): # only download if file does not exist already
        return
    
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


if __name__ == "__main__":
    func_list = [fm_grevelingen_map, fm_grevelingen_his, fm_grevelingen_net, 
                 fm_curvedbend_map, fm_curvedbend_his, 
                 #fm_westernscheldt_map, 
                 d3d_westernscheldt_trim]
    for func in func_list:
        file_nc = func(return_filepath=True)
        print(file_nc)
        uds = func()


