# -*- coding: utf-8 -*-
"""
Created on Fri Jul 14 12:04:27 2023

@author: veenstra
"""

import pytest
import xarray as xr
import xugrid as xu
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


@pytest.mark.requireslocaldata
@pytest.mark.unittest
def test_open_2Dnetwork_with_1Dtopology():
    file_nc = r'p:\1204257-dcsmzuno\2006-2012\3D-DCSM-FM\A18b_ntsu1\DCSM-FM_0_5nm_grid_20191202_depth_20181213_net.nc'
    uds = dfmt.open_partitioned_dataset(file_nc)
    assert isinstance(uds.grid, xu.ugrid.ugrid1d.Ugrid1d)
