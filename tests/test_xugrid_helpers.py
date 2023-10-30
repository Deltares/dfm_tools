# -*- coding: utf-8 -*-
"""
Created on Fri Jul 14 12:04:27 2023

@author: veenstra
"""

import pytest
import xarray as xr
import xugrid as xu
import dfm_tools as dfmt
import numpy as np

#TODO: many xugrid_helpers tests are still in test_dfm_tools.py


@pytest.mark.unittest
def test_remove_unassociated_edges():
    file_nc = dfmt.data.fm_westernscheldt_map(return_filepath=True) #sigmalayer
    ds = xr.open_dataset(file_nc)
    ds2 = dfmt.remove_unassociated_edges(ds)
    
    ds_edgedimsize = ds.dims['nmesh2d_edge']
    ds2_edgedimsize = ds2.dims['nmesh2d_edge']
    
    assert ds2_edgedimsize == ds_edgedimsize-1


@pytest.mark.unittest
def test_open_2Dnetwork_with_1Dtopology():
    file_nc = r'p:\1204257-dcsmzuno\2006-2012\3D-DCSM-FM\A18b_ntsu1\DCSM-FM_0_5nm_grid_20191202_depth_20181213_net.nc'
    uds = dfmt.open_partitioned_dataset(file_nc)
    assert isinstance(uds.grid, xu.ugrid.ugrid1d.Ugrid1d)
    assert not hasattr(uds.grid, "face_node_connectivity")
    
    uds_withcellinfo = dfmt.add_network_cellinfo(uds)
    assert isinstance(uds_withcellinfo.grid, xu.Ugrid2d)
    assert hasattr(uds_withcellinfo.grid, "face_node_connectivity")


@pytest.mark.systemtest
def test_uda_edges_to_faces_interfaces_to_centers():
    file_nc = dfmt.data.fm_grevelingen_map(return_filepath=True) #zlayer
    
    uds = xu.open_dataset(file_nc.replace('0*','0002')) #partition 0002 of grevelingen contains both triangles as squares
    dimn_faces = uds.grid.face_dimension
    dimn_layer, _ = dfmt.get_vertical_dimensions(uds)
    
    for varn_edge in ['mesh2d_vicwwu','mesh2d_edge_type']:
        #vicwwu includes interface to layer interpolation
        #edge_type array is completely different when masking is forgotten (when uda contains both triangles and squares)
        uda_edge = uds[varn_edge]
        
        uda_face_int = dfmt.uda_edges_to_faces(uda_edge)
        uda_face = dfmt.uda_interfaces_to_centers(uda_face_int)

        if varn_edge == 'mesh2d_vicwwu':
            uda_face_sel = uda_face.isel({'time':-1, dimn_faces:slice(None,20), dimn_layer:-1})
            uda_face_sel_expected = np.array([0.000385  , 0.00036918, 0.00037161, 0.00055217, 0.00048461,
                    0.00043027, 0.00042121, 0.00063324, 0.00044923, 0.00052215,
                    0.00042954, 0.0004004 , 0.00051421, 0.00049618, 0.00051311,
                    0.00076876, 0.00089225, 0.00065682, 0.00069086, 0.000523  ])
            uds_face_dims_expected = ('time', dimn_faces, dimn_layer)
        elif varn_edge == 'mesh2d_edge_type':
            uda_face_sel = uda_face.isel({dimn_faces:slice(None,20)})
            uda_face_sel_expected = np.array([1.66666667, 1.66666667, 1.66666667, 1.66666667, 1.66666667,
                   1.        , 1.        , 1.        , 1.66666667, 1.66666667,
                   1.66666667, 1.66666667, 1.        , 1.        , 1.        ,
                   1.        , 1.        , 1.66666667, 1.66666667, 1.        ])
            uds_face_dims_expected = (dimn_faces,)
        
        assert uda_face.dims == uds_face_dims_expected
    
        assert (np.abs(uda_face_sel.to_numpy()-uda_face_sel_expected)<1e-6).all()
