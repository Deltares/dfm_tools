# -*- coding: utf-8 -*-
"""
Created on Wed Dec 13 12:17:41 2023

@author: veenstra
"""


import pytest
import numpy as np
import dfm_tools as dfmt


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
