# -*- coding: utf-8 -*-
"""
Created on Fri Jul 14 12:04:27 2023

@author: veenstra
"""

import os
import pytest
import xarray as xr
import xugrid as xu
import dfm_tools as dfmt
import numpy as np
from dfm_tools.xugrid_helpers import (remove_unassociated_edges,
                                      get_vertical_dimensions
                                      )


@pytest.mark.unittest
def test_open_partitioned_dataset():
    """
    Checks whether ghost cells are properly taken care of by asserting shape
    """
    file_nc = dfmt.data.fm_grevelingen_map(return_filepath=True)
    data_xr_map = dfmt.open_partitioned_dataset(file_nc)
    data_varsel = data_xr_map['mesh2d_sa1'].isel(time=2)
    assert data_varsel.shape == (44796, 36)


@pytest.mark.parametrize("file_nc, expected_size", [pytest.param(dfmt.data.fm_grevelingen_map(return_filepath=True), (44796, 4, 2), id='partitioned mapfile Grevelingen'),
                                                    pytest.param(dfmt.data.fm_curvedbend_map(return_filepath=True), (550, 4, 2), id='mapfile curvedbend'),
                                                    pytest.param(dfmt.data.fm_grevelingen_net(return_filepath=True), (44804,4,2), id='network Grevelingen')])
@pytest.mark.unittest
def test_facenodecoordinates_shape(file_nc, expected_size):
    
    uds = dfmt.open_partitioned_dataset(file_nc)
    facenodecoordinates = uds.grid.face_node_coordinates
    
    assert facenodecoordinates.shape == expected_size


@pytest.mark.unittest
def test_remove_unassociated_edges():
    file_nc = dfmt.data.fm_westernscheldt_map(return_filepath=True) #sigmalayer
    ds = xr.open_dataset(file_nc)
    ds2 = remove_unassociated_edges(ds)
    
    ds_edgedimsize = ds.sizes['mesh2d_nEdges']
    ds2_edgedimsize = ds2.sizes['mesh2d_nEdges']
    
    assert ds2_edgedimsize == ds_edgedimsize-1


@pytest.mark.unittest
def test_remove_nan_fillvalue_attrs(tmp_path):
    """
    xarray writes {"_FillValue": np.nan} to encoding for variables without _FillValue attribute.
    This test checks if that is still the case and checks if dfmt.open_partitioned_dataset removes them.
    """
    file_nc = dfmt.data.fm_curvedbend_map(return_filepath=True)
    file_out = os.path.join(tmp_path, "temp_fillvals_map.nc")
    ds_org = xr.open_dataset(file_nc)
    ds_org.to_netcdf(file_out)
    
    ds_out_xr = xr.open_dataset(file_out)
    ds_out_dfmt = dfmt.open_partitioned_dataset(file_out, chunks="auto")
    
    print("nan fillvalue attrs in dataset written by xugrid/xarray")
    ds = ds_out_xr
    count_xr = 0
    for varn in ds.variables:
        if '_FillValue' in ds[varn].encoding:
            if np.isnan(ds[varn].encoding['_FillValue']):
                print(varn, ds[varn].encoding['_FillValue'])
                count_xr += 1
    
    print("nan fillvalue attrs in dataset written by xugrid/xarray, read with dfm_tools")
    ds = ds_out_dfmt
    count_dfmt = 0
    for varn in ds.variables:
        if '_FillValue' in ds[varn].encoding:
            if np.isnan(ds[varn].encoding['_FillValue']):
                print(varn, ds[varn].encoding['_FillValue'])
                count_dfmt += 1
    
    assert count_xr == 10
    assert count_dfmt == 0


@pytest.mark.unittest
def test_uds_auto_set_crs_cartesian():
    file_nc = dfmt.data.fm_grevelingen_map(return_filepath=True).replace('0*','0002')
    uds = dfmt.open_partitioned_dataset(file_nc)
    assert uds.ugrid.crs['mesh2d'] is not None
    assert uds.ugrid.crs['mesh2d'].to_epsg() == 28992
    assert uds.grid.is_geographic is False


@pytest.mark.unittest
def test_uds_auto_set_crs_none():
    file_nc = dfmt.data.fm_curvedbend_map(return_filepath=True)
    uds = dfmt.open_partitioned_dataset(file_nc)
    assert uds.ugrid.crs['mesh2d'] is None
    assert uds.grid.is_geographic is False


@pytest.mark.unittest
def test_uds_auto_set_crs_spherical():
    # using dummy dataset, was tested before with
    # 'p:\1204257-dcsmzuno\2006-2012\3D-DCSM-FM\A18b_ntsu1\DCSM-FM_0_5nm_grid_20191202_depth_20181213_net.nc'
    file_nc = dfmt.data.fm_curvedbend_map(return_filepath=True)
    def replace_crs(ds):
        ds = ds.drop_vars('projected_coordinate_system')
        spherical_attrs = {'name': 'WGS84',
         'epsg': np.int32(4326),
         'grid_mapping_name': 'latitude_longitude',
         'EPSG_code': 'EPSG:4326',
         }
        ds['wgs84'] = xr.DataArray(0).assign_attrs(spherical_attrs)
        return ds
    uds = dfmt.open_partitioned_dataset(file_nc, preprocess=replace_crs)
    assert uds.ugrid.crs['mesh2d'] is not None
    assert uds.ugrid.crs['mesh2d'].to_epsg() == 4326
    assert uds.grid.is_geographic is True


@pytest.mark.unittest
def test_open_2Dnetwork_with_1Dtopology(tmp_path):
    file_nc = dfmt.data.fm_grevelingen_map(return_filepath=True).replace('0*','0002')
    
    # save as new network with only the 1D information in the mesh_topology
    mesh1d_attrs = {'cf_role': 'mesh_topology',
     'long_name': 'Topology data of 1D network',
     'topology_dimension': 1,
     'node_coordinates': 'mesh2d_node_x mesh2d_node_y',
     'node_dimension': 'nmesh2d_node',
     'edge_node_connectivity': 'mesh2d_edge_nodes',
     'edge_dimension': 'nmesh2d_edge'}
    ds = xr.open_dataset(file_nc, decode_cf=False, decode_times=False, decode_coords=False)
    ds['mesh2d'] = ds['mesh2d'].assign_attrs(mesh1d_attrs)
    file_nc_out = os.path.join(tmp_path, "network_1d_net.nc")
    ds.to_netcdf(file_nc_out)
    
    uds = dfmt.open_partitioned_dataset(file_nc_out)
    assert isinstance(uds.grid, xu.ugrid.ugrid1d.Ugrid1d)
    assert not hasattr(uds.grid, "face_node_connectivity")
    
    uds_withcellinfo = dfmt.add_network_cellinfo(uds)
    assert isinstance(uds_withcellinfo.grid, xu.Ugrid2d)
    assert hasattr(uds_withcellinfo.grid, "face_node_connectivity")
    
    # test projected attr, xugrid automatically sets standard_name 
    # based on projected property
    node_x = uds_withcellinfo.grid.to_dataset().mesh2d_node_x
    assert node_x.attrs['standard_name'] == "projection_x_coordinate"


@pytest.mark.unittest
def test_uda_edges_to_faces():
    file_nc = dfmt.data.fm_grevelingen_map(return_filepath=True) #zlayer
    
    uds = xu.open_dataset(file_nc.replace('0*','0002')) #partition 0002 of grevelingen contains both triangles as squares
    
    uda_edge = uds['mesh2d_vicwwu'].isel(time=0, nmesh2d_interface=0)
    uda_face = dfmt.uda_to_faces(uda_edge)
    
    assert uda_face.dims == (uds.grid.face_dimension,)
    assert hasattr(uda_face,'grid')


@pytest.mark.unittest
def test_get_vertical_dimensions():
    # hisfile
    ds = dfmt.data.fm_curvedbend_his()
    dimn_layer, dimn_interface = get_vertical_dimensions(ds)
    assert dimn_layer is None
    assert dimn_interface is None
   
    # mapfile
    file_nc = dfmt.data.fm_grevelingen_map(return_filepath=True)
    uds = xu.open_dataset(file_nc.replace('0*','0002')) #partition 0002 of grevelingen contains both triangles as squares
    dimn_layer, dimn_interface = get_vertical_dimensions(uds)
    assert dimn_layer == "nmesh2d_layer"
    assert dimn_interface == "nmesh2d_interface"

    # mapfile
    file_nc = dfmt.data.fm_curvedbend_map(return_filepath=True)
    uds = xu.open_dataset(file_nc)
    dimn_layer, dimn_interface = get_vertical_dimensions(uds)
    assert dimn_layer == "mesh2d_nLayers"
    assert dimn_interface == "mesh2d_nInterfaces"

    # networkfile (2D)
    file_nc = dfmt.data.fm_grevelingen_net(return_filepath=True)
    uds = xu.open_dataset(file_nc)
    dimn_layer, dimn_interface = get_vertical_dimensions(uds)
    assert dimn_layer is None
    assert dimn_interface is None


@pytest.mark.systemtest
def test_uda_edges_to_faces_interfaces_to_centers():
    file_nc = dfmt.data.fm_grevelingen_map(return_filepath=True) #zlayer
    
    uds = xu.open_dataset(file_nc.replace('0*','0002')) #partition 0002 of grevelingen contains both triangles as squares
    dimn_faces = uds.grid.face_dimension
    dimn_layer, _ = get_vertical_dimensions(uds)
    
    for varn_edge in ['mesh2d_vicwwu','mesh2d_edge_type']:
        #vicwwu includes interface to layer interpolation
        #edge_type array is completely different when masking is forgotten (when uda contains both triangles and squares)
        uda_edge = uds[varn_edge]
        
        uda_face_int = dfmt.uda_to_faces(uda_edge)
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
        
        assert set(uda_face.dims) == set(uds_face_dims_expected)
        assert hasattr(uda_face,'grid')
        
        assert (np.abs(uda_face_sel.to_numpy()-uda_face_sel_expected)<1e-6).all()


@pytest.mark.unittest
def test_uda_nodes_to_faces():
    file_nc = dfmt.data.fm_grevelingen_net(return_filepath=True)
    uds = dfmt.open_partitioned_dataset(file_nc)
    
    uda_node = uds.mesh2d_node_z
    uda_face = dfmt.uda_to_faces(uda_node)
    
    assert uda_face.dims == (uds.grid.face_dimension,)
    assert hasattr(uda_face,'grid')


@pytest.mark.unittest
def test_enrich_rst_with_map():
    """
    this tests whether rst file is correctly enriched with map by trying to
    plot a ugrid (face node) variable was tested before with
    'p:\\archivedprojects\\11206811-002-kpp-veerse-meer\\model\\runs_2011-2012\\VM_WQ_3D_run9_c\\DFM_OUTPUT_VM_WQ_3D\\VM_WQ_3D_0000_20130101_000000_rst.nc'
    
    """
    # mf1_rstfile (without topology var)
    
    file_nc_rst = dfmt.data.fm_westernscheldt_rst(return_filepath=True)
    
    uds_rst = dfmt.open_partitioned_dataset(file_nc_rst, preprocess=dfmt.enrich_rst_with_map)
    
    uds_rst.s1.isel(time=0).ugrid.plot()
    
    uda_rst = uds_rst['s1']
    assert "mesh2d_face_x" in uda_rst.coords


@pytest.mark.unittest
def test_open_dataset_delft3d4():
    file_nc = dfmt.data.d3d_curvedbend_trim(return_filepath=True)

    uds = dfmt.open_dataset_delft3d4(file_nc)

    assert "mesh2d" in uds.grid.to_dataset().data_vars
    assert "grid" not in uds.data_vars

    # test if plotting works, this is a basic validation of whether it is a proper ugrid dataset
    uds.umag.isel(time=-1, KMAXOUT_RESTR=-1).ugrid.plot()


@pytest.mark.requireslocaldata
@pytest.mark.unittest
def test_open_dataset_curvilinear():
    file_nc = r'P:\archivedprojects\11206304-futuremares-rawdata-preps\data\CMIP6_BC\CMCC-ESM2\vo_Omon_CMCC-ESM2_historical_r1i1p1f1_gn_*.nc'
    
    with pytest.warns(UserWarning) as w:
        uds = dfmt.open_dataset_curvilinear(
            file_nc,
            x_dim='i',
            y_dim='j',
            x_bounds='vertices_longitude',
            y_bounds='vertices_latitude',
            convert_360to180=True,
            )
    assert "A UGRID2D face requires at least" in str(w[0].message)

    # check if vertices lat/lon and lat/lon variables are dropped: https://github.com/Deltares/dfm_tools/issues/930
    assert set(uds.coords) == set(['latitude', 'longitude', 'lev', 'mesh2d_nFaces', 'time'])
    assert set(uds.data_vars) == set(['lev_bnds', 'time_bnds', 'vo'])
    
    # check absence of time dimension where it should not be: https://github.com/Deltares/dfm_tools/issues/928
    assert uds.lev_bnds.dims == ('lev', 'bnds')
    
    # check if all zero-sized cells were dropped: https://github.com/Deltares/dfm_tools/issues/926
    assert (uds.grid.area > 0).all()
