# -*- coding: utf-8 -*-
"""
Created on Wed Jul  5 12:59:23 2023

@author: veenstra
"""

import os
import pytest
import xugrid as xu
import dfm_tools as dfmt
import meshkernel
import xarray as xr
import numpy as np


@pytest.mark.unittest
def test_add_crs_to_dataset_cartesian():
    uds = xu.data.adh_san_diego()
    crs='EPSG:26946' # this is not the correct crs for this model, but that does not matter
    dfmt.add_crs_to_dataset(uds,is_geographic=False,crs=crs)
    
    assert 'projected_coordinate_system' in uds.data_vars
    crs_attrs = uds.projected_coordinate_system.attrs
    assert crs_attrs['EPSG_code'] == 'EPSG:26946'
    assert crs_attrs['epsg'] == 26946
    assert crs_attrs['grid_mapping_name'] == 'Unknown projected'


@pytest.mark.unittest
def test_add_crs_to_dataset_spherical():
    uds = xu.data.adh_san_diego()
    crs='EPSG:4326' # this is not the correct crs for this model, but that does not matter
    dfmt.add_crs_to_dataset(uds,is_geographic=True,crs=crs)
    
    assert 'wgs84' in uds.data_vars
    crs_attrs = uds.wgs84.attrs
    assert crs_attrs['EPSG_code'] == 'EPSG:4326'
    assert crs_attrs['epsg'] == 4326
    assert crs_attrs['grid_mapping_name'] == 'latitude_longitude'
    

@pytest.mark.systemtest
def test_meshkernel_check_geographic():
    """
    to check whether is_geographic can be correctly derived from the mk object, for cartesian as well spherical objects
    """
    lon_min, lon_max, lat_min, lat_max = -68.55, -67.9, 11.8, 12.6
    dxy = 0.05
    
    mk_cartesian = dfmt.make_basegrid(lon_min, lon_max, lat_min, lat_max, dx=dxy, dy=dxy, is_geographic=False)
    mk_cartesian_geograph = dfmt.meshkernel_check_geographic(mk_cartesian)
    mk_spherical = dfmt.make_basegrid(lon_min, lon_max, lat_min, lat_max, dx=dxy, dy=dxy, is_geographic=True)
    mk_spherical_geograph = dfmt.meshkernel_check_geographic(mk_spherical)
    
    assert mk_cartesian_geograph==False
    assert mk_spherical_geograph==True


@pytest.mark.unittest
def test_meshkernel_delete_withcoastlines():
    #generate basegrid
    lon_min, lon_max, lat_min, lat_max = -68.45, -68.1, 12, 12.35
    dxy = 0.005
    mk = dfmt.make_basegrid(lon_min, lon_max, lat_min, lat_max, dx=dxy, dy=dxy)
    
    assert len(mk.mesh2d_get().face_nodes) == 20732
    
    # remove cells with GSHHS coastlines
    dfmt.meshkernel_delete_withcoastlines(mk=mk, res='h')
    
    assert len(mk.mesh2d_get().face_nodes) == 17364


@pytest.mark.unittest
def test_meshkernel_delete_withgdf():
    #generate basegrid
    lon_min, lon_max, lat_min, lat_max = -68.45, -68.1, 12, 12.35
    dxy = 0.005
    mk = dfmt.make_basegrid(lon_min, lon_max, lat_min, lat_max, dx=dxy, dy=dxy)
    
    assert len(mk.mesh2d_get().face_nodes) == 20732
    
    # remove cells with custom ldb_gdf (can also come from polyfile)
    bbox = (lon_min, lat_min, lon_max, lat_max)
    ldb_gdf = dfmt.get_coastlines_gdb(bbox=bbox, res='h')
    dfmt.meshkernel_delete_withgdf(mk=mk, coastlines_gdf=ldb_gdf)
    
    assert len(mk.mesh2d_get().face_nodes) == 17364


@pytest.mark.systemtest
def test_meshkernel_to_UgridDataset():
    """
    generate grid with meshkernel. Then convert with `dfmt.meshkernel_to_UgridDataset()` from 0-based to 1-based indexing to make FM-compatible network.
    assert if _FillValue, start_index, min and max are the expected values, this ensures FM-compatibility
    """
    is_geographic = False #TODO: polygon refinement does not work for spherical grids: https://github.com/Deltares/MeshKernelPy/issues/78
    crs = 'EPSG:28992' #arbitrary non-spherical epsg code
    
    # create basegrid
    lon_min, lon_max, lat_min, lat_max = -6, 2, 48.5, 51.2
    dxy = 0.5
    make_grid_parameters = meshkernel.MakeGridParameters(origin_x=lon_min,
                                                         origin_y=lat_min,
                                                         upper_right_x=lon_max,
                                                         upper_right_y=lat_max,
                                                         block_size_x=dxy,
                                                         block_size_y=dxy)
    mk = meshkernel.MeshKernel(is_geographic=is_geographic)
    mk.curvilinear_make_uniform_on_extension(make_grid_parameters)
    mk.curvilinear_convert_to_mesh2d() #convert to ugrid/mesh2d
    
    # refine with polygon
    pol_x = np.array([-5,-4,0,-5], dtype=np.double)
    pol_y = np.array([49,51,49.5,49], dtype=np.double)
    geometry_list = meshkernel.GeometryList(pol_x, pol_y)
    mrp = meshkernel.MeshRefinementParameters()
    mk.mesh2d_refine_based_on_polygon(polygon=geometry_list, mesh_refinement_params=mrp)
    
    #convert to xugrid and write to netcdf
    xu_grid_uds = dfmt.meshkernel_to_UgridDataset(mk=mk, crs=crs)
    netfile = 'test_startindex_net.nc'
    xu_grid_uds.ugrid.to_netcdf(netfile)
    
    # plot
    # fig,ax = plt.subplots()
    # xu_grid_uds.grid.plot(ax=ax)
    # ax.plot(pol_x,pol_y,'r-')
    
    #assert output grid
    ds_out = xr.open_dataset(netfile,decode_cf=False).load()
    ds_out.close()
    os.remove(netfile)
    assert ds_out.mesh2d_face_nodes.attrs['_FillValue'] == -1
    assert ds_out.mesh2d_face_nodes.attrs['start_index'] == 1
    assert 0 not in ds_out.mesh2d_face_nodes.to_numpy()
    assert ds_out.mesh2d_face_nodes.to_numpy().min() == -1
    assert ds_out.mesh2d_face_nodes.to_numpy().max() == 135

