# -*- coding: utf-8 -*-
"""
Created on Mon Jun 19 22:18:46 2023

@author: veenstra
"""

import os
import requests
import pathlib


def maybe_download_testdata():
    #TODO: work with pooch instead, like: https://github.com/Deltares/xugrid/blob/main/xugrid/data/sample_data.py
    #TODO: make opendap folder structure the same as local testdata folder
    
    #find dir_testinput
    dir_testinput = os.path.join(r'c:\DATA','dfm_tools_testdata') #on WCF
    if os.path.exists(dir_testinput):
        print(f'local testdata found in {dir_testinput}')
        return dir_testinput
    dir_testinput = './dfm_tools_testdata' #for instance when running on github
    if os.path.exists(dir_testinput): #escape when data is already downloaded
        print(f'local testdata found in {dir_testinput}')
        return dir_testinput
    
    print(f'no local testdata found, files are downloaded to {dir_testinput}')
    
    fname_list = []
    fname_list += ['DFM_curvedbend_3D/cb_3d_map.nc']
    fname_list += ['DFM_curvedbend_3D/cb_3d_his.nc']
    fname_list += [f'DFM_grevelingen_3D/Grevelingen-FM_{i:04d}_map.nc' for i in range(8)]
    fname_list += ['DFM_grevelingen_3D/Grevelingen-FM_0000_his.nc']
    fname_list += ['DFM_grevelingen_3D/Grevelingen_FM_grid_20190603_net.nc']
    # fname_list += ['westernscheldt_sph_map.nc']

    for fname in fname_list:
        file_nc = os.path.join(dir_testinput,fname)
        file_dirname = os.path.dirname(file_nc)
        if os.path.exists(file_nc): #skip if file exists
            continue
        pathlib.Path(file_dirname).mkdir(parents=True, exist_ok=True) #make subdir if needed
        
        file_url = f'https://opendap.deltares.nl/thredds/fileServer/opendap/deltares/Delft3D/netcdf_example_files/{fname}'
        print(f'downloading {file_url} to {file_nc}')
        r = requests.get(file_url, allow_redirects=True)
        r.raise_for_status() #raise HTTPError if url not exists
        with open(file_nc, 'wb') as f:
            f.write(r.content)
    
    return dir_testinput
