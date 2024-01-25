# -*- coding: utf-8 -*-
"""
Created on Wed Dec 13 12:17:41 2023

@author: veenstra
"""


import pytest
import numpy as np
import dfm_tools as dfmt
import xarray as xr


@pytest.mark.unittest
def test_polyline_mapslice():
    uds = dfmt.data.fm_curvedbend_map()
    timestep = 72
    line_array = np.array([[ 104.15421399, 2042.7077107 ],
                            [2913.47878063, 2102.48057382]])
    
    uds_crs = dfmt.polyline_mapslice(uds.isel(time=timestep), line_array)
    assert len(uds_crs.grid.node_x) == 720
    assert uds_crs.grid.face_node_connectivity.shape == (180, 4)
    assert np.isclose(uds_crs.grid.node_x.min(), 484.02691589067194)
    assert np.isclose(uds_crs.grid.node_x.max(), 1737.0864257281437)
    assert np.isclose(uds_crs.grid.node_y.min(), -5)
    assert np.isclose(uds_crs.grid.node_y.max(), 0.9261683648147339)


def test_get_dataset_atdepths_hisfile():
    
    file_nc = dfmt.data.fm_grevelingen_his(return_filepath=True)
    ds = xr.open_dataset(file_nc)#, preprocess=dfmt.preprocess_hisnc)

    depths = [-1,-4,0,-6]
    data_fromhis_atdepths = dfmt.get_Dataset_atdepths(data_xr=ds, depths=depths, reference='z0')
    data_xr_selzt = data_fromhis_atdepths.isel(stations=2).isel(time=slice(40,100))

    tem_values = data_xr_selzt.temperature.isel(time=-1).to_numpy()
    exp_values = np.array([5.81065253, 5.43289777, 5.36916911, 5.36916911])

    assert data_xr_selzt.temperature.shape == (60,4)
    assert np.allclose(tem_values, exp_values)
    assert np.isclose(data_xr_selzt.temperature.sum(), 1295.56826688)


