# -*- coding: utf-8 -*-
"""
Created on Sun Jul  9 22:17:16 2023

@author: veenstra
"""

import pytest
import xarray as xr
import xugrid as xu
import dfm_tools as dfmt


@pytest.mark.unittest
def test_data_map():
    func_list = [dfmt.data.fm_grevelingen_map, 
                 dfmt.data.fm_grevelingen_net, 
                 dfmt.data.fm_curvedbend_map, 
                 dfmt.data.fm_westernscheldt_map, 
                 dfmt.data.fm_westernscheldt_fou, 
                 dfmt.data.fm_westernscheldt_rst, 
                 dfmt.data.d3d_westernscheldt_trim,
                 dfmt.data.d3d_curvedbend_trim]
    
    for func in func_list:
        file_nc = func(return_filepath=True)
        uds = func()
        assert isinstance(file_nc,str)
        assert isinstance(uds,xu.UgridDataset)


@pytest.mark.unittest
def test_data_his():
    func_list = [dfmt.data.fm_grevelingen_his, 
                 dfmt.data.fm_curvedbend_his,
                 dfmt.data.fm_westernscheldt_his,
                 dfmt.data.d3d_curvedbend_trih]
    
    for func in func_list:
        file_nc = func(return_filepath=True)
        uds = func()
        assert isinstance(file_nc,str)
        assert isinstance(uds,xr.Dataset)


@pytest.mark.unittest
def test_gshhs_coastlines_shp():
    gshhs_dir = dfmt.data.gshhs_coastlines_shp()
    assert isinstance(gshhs_dir,str)
