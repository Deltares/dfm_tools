# -*- coding: utf-8 -*-
"""
Created on Fri Jul 14 12:04:27 2023

@author: veenstra
"""

import pytest
import xarray as xr
import dfm_tools as dfmt

#TODO: many xugrid_helpers tests are still in test_dfm_tools.py


@pytest.mark.unittest
def test_remove_unassociated_edges():
    file_nc = dfmt.data.fm_westernscheldt_map(return_filepath=True) #sigmalayer
    ds = xr.open_dataset(file_nc)
    ds2 = dfmt.remove_unassociated_edges(ds)
    
    ds_edgedimsize = ds.dims['nmesh2d_edge']
    ds2_edgedimsize = ds2.dims['nmesh2d_edge']
    
    assert ds2_edgedimsize == ds_edgedimsize-1
