# -*- coding: utf-8 -*-
"""
Created on Wed Jul  5 12:59:23 2023

@author: veenstra
"""

import pytest
import xugrid as xu
import dfm_tools as dfmt
import meshkernel
import xarray as xr
import numpy as np
import geopandas as gpd
from shapely.geometry import Polygon
from dfm_tools.meshkernel_helpers import (geographic_to_meshkernel_projection, 
                                          uds_add_crs_attrs,
                                          crs_to_isgeographic
                                          )


@pytest.mark.unittest
def test_crs_to_isgeographic():
    is_geographic = crs_to_isgeographic(4326)
    assert is_geographic is True
    is_geographic = crs_to_isgeographic("EPSG:4326")
    assert is_geographic is True
    is_geographic = crs_to_isgeographic("WGS84")
    assert is_geographic is True

    is_geographic = crs_to_isgeographic(28992)
    assert is_geographic is False
    is_geographic = crs_to_isgeographic("EPSG:28992")
    assert is_geographic is False
    is_geographic = crs_to_isgeographic("Amersfoort / RD New")
    assert is_geographic is False
    
    is_geographic = crs_to_isgeographic(None)
    assert is_geographic is False


@pytest.mark.unittest
def test_uds_add_crs_attrs_cartesian():
    uds = xu.data.adh_san_diego()
    crs='EPSG:28992' # this is not the correct crs for this model, but that does not matter
    uds.ugrid.set_crs(crs)
    uds_add_crs_attrs(uds)
    
    assert "projected_coordinate_system" in uds.data_vars
    crs_attrs = uds["projected_coordinate_system"].attrs
    assert crs_attrs["name"] == "Amersfoort / RD New"
    assert crs_attrs["epsg"] == 28992
    assert crs_attrs["EPSG_code"] == "EPSG:28992"
    assert "grid_mapping_name" not in crs_attrs.keys()


@pytest.mark.unittest
def test_uds_add_crs_attrs_spherical():
    uds = xu.data.adh_san_diego()
    crs='EPSG:4326' # this is not the correct crs for this model, but that does not matter
    uds.ugrid.set_crs(crs)
    uds_add_crs_attrs(uds)
    
    assert "wgs84" in uds.data_vars
    crs_attrs = uds["wgs84"].attrs
    assert crs_attrs["name"] == "WGS 84"
    assert crs_attrs["epsg"] == 4326
    assert crs_attrs["EPSG_code"] == "EPSG:4326"
    assert "grid_mapping_name" in crs_attrs.keys()
    assert crs_attrs["grid_mapping_name"] == "latitude_longitude"


@pytest.mark.unittest
def test_meshkernel_refine_basegrid():
    # domain and resolution
    lon_min, lon_max, lat_min, lat_max = -68.55, -67.9, 11.8, 12.6
    dxy = 0.05
    crs = 'EPSG:4326'

    # grid generation and refinement with GEBCO bathymetry
    file_nc_bathy = r'https://opendap.deltares.nl/thredds/dodsC/opendap/deltares/Delft3D/netcdf_example_files/GEBCO_2022/GEBCO_2022_coarsefac08.nc'
    data_bathy = xr.open_dataset(file_nc_bathy)
    
    # TODO: in the future meshkernel supports all these four datatypes, so add more datatypes to list
    for dtype in ["float32", "int16"]: # ["float64", "float32", "int32", "int16"]:
        data_bathy_sel = data_bathy.sel(lon=slice(lon_min,lon_max),lat=slice(lat_min,lat_max)).elevation
        data_bathy_sel = data_bathy_sel.astype(dtype=dtype)
        
        # basegrid
        mk_object = dfmt.make_basegrid(lon_min, lon_max, lat_min, lat_max, dx=dxy, dy=dxy, crs=crs)
    
        # refine (fails if wrong type)
        min_edge_size = 300 #in meters
        dfmt.refine_basegrid(mk=mk_object, data_bathy_sel=data_bathy_sel, min_edge_size=min_edge_size)

        assert len(mk_object.mesh2d_get().node_x) == 4074
        assert len(mk_object.mesh2d_get().face_nodes) == 17298


@pytest.mark.unittest
def test_meshkernel_delete_withcoastlines():
    #generate basegrid
    lon_min, lon_max, lat_min, lat_max = -68.45, -68.1, 12, 12.35
    dxy = 0.005
    mk = dfmt.make_basegrid(lon_min, lon_max, lat_min, lat_max, dx=dxy, dy=dxy, crs=4326)
    assert mk.get_projection() == meshkernel.ProjectionType.SPHERICAL
    assert len(mk.mesh2d_get().face_nodes) == 20732
    
    # remove cells with GSHHS coastlines
    dfmt.meshkernel_delete_withcoastlines(mk=mk, res='h')
    
    assert len(mk.mesh2d_get().face_nodes) == 17368


@pytest.mark.unittest
def test_meshkernel_delete_withshp(tmp_path):
    # write shapefile from coords
    file_shp = tmp_path / 'mk_delete_test.shp'
    points_x = np.array([-68.40156631636155, -68.36143523088661, -68.28392176442131, -68.26413109213229,  
                         -68.20915700244058, -68.1965129618115, -68.20860726154366, -68.199811407193,    
                         -68.23059689742034, -68.23389534280184, -68.26303161033846, -68.29436684146273, 
                         -68.27512591007063, -68.3064611411949, -68.37462901241263, -68.41146165250606, 
    					 -68.41311087519682, -68.40156631636155])
    points_y = np.array([12.301796772360385, 12.303445995051137, 12.243524237287176, 12.233628901142668, 
                         12.225382787688911, 12.206141856296814, 12.184152220420131, 12.115984349202414, 
                         12.080800931799722, 12.028025805695682, 12.033523214664854, 12.113785385614745, 
                         12.144021134945184, 12.20119418822456,  12.220984860513575, 12.226482269482746, 
                         12.286404027246707, 12.301796772360385])
    
    geom = Polygon(zip(points_x, points_y))
    gdf = gpd.GeoDataFrame(geometry=[geom], crs='EPSG:4326')
    gdf.to_file(file_shp)
    
    #generate basegrid
    lon_min, lon_max, lat_min, lat_max = -68.45, -68.1, 12, 12.35
    dxy = 0.005
    mk = dfmt.make_basegrid(lon_min, lon_max, lat_min, lat_max, dx=dxy, dy=dxy, crs=4326)
    assert mk.get_projection() == meshkernel.ProjectionType.SPHERICAL
    
    assert len(mk.mesh2d_get().face_nodes) == 20732
    
    # remove cells with a shapefile
    dfmt.meshkernel_delete_withshp(mk=mk, coastlines_shp=file_shp)
    
    assert len(mk.mesh2d_get().face_nodes) == 17272


@pytest.mark.unittest
def test_meshkernel_delete_withgdf():
    #generate basegrid
    lon_min, lon_max, lat_min, lat_max = -68.45, -68.1, 12, 12.35
    dxy = 0.005
    mk = dfmt.make_basegrid(lon_min, lon_max, lat_min, lat_max, dx=dxy, dy=dxy, crs=4326)
    assert mk.get_projection() == meshkernel.ProjectionType.SPHERICAL

    assert len(mk.mesh2d_get().face_nodes) == 20732
    
    # remove cells with custom ldb_gdf (can also come from polyfile)
    bbox = (lon_min, lat_min, lon_max, lat_max)
    ldb_gdf = dfmt.get_coastlines_gdb(bbox=bbox, res='h')
    dfmt.meshkernel_delete_withgdf(mk=mk, coastlines_gdf=ldb_gdf)
    
    assert len(mk.mesh2d_get().face_nodes) == 17368


@pytest.mark.unittest
def test_geographic_to_meshkernel_projection():
    
    spherical = geographic_to_meshkernel_projection(is_geographic=True)
    cartesian = geographic_to_meshkernel_projection(is_geographic=False)
    spherical_mk = meshkernel.ProjectionType.SPHERICAL
    cartesian_mk = meshkernel.ProjectionType.CARTESIAN
    
    assert spherical == spherical_mk
    assert cartesian == cartesian_mk


@pytest.mark.systemtest
def test_meshkernel_to_UgridDataset_geographic_mismatch():
    """
    """
    projection = meshkernel.ProjectionType.CARTESIAN
    crs = 'EPSG:4326'
    mk = meshkernel.MeshKernel(projection=projection)
    try:
        _ = dfmt.meshkernel_to_UgridDataset(mk=mk, crs=crs)
    except ValueError as e:
        assert "crs has is_geographic=True and grid has is_geographic=False" in str(e)


@pytest.mark.systemtest
def test_meshkernel_to_UgridDataset(tmp_path):
    """
    generate grid with meshkernel. Then convert with `dfmt.meshkernel_to_UgridDataset()` from 0-based to 1-based indexing to make FM-compatible network.
    assert if _FillValue, start_index, min and max are the expected values, this ensures FM-compatibility
    """
    projection = meshkernel.ProjectionType.SPHERICAL
    crs = 'EPSG:4326'
    
    # create basegrid
    lon_min, lon_max, lat_min, lat_max = -6, 2, 48.5, 51.2
    dxy = 0.5
    make_grid_parameters = meshkernel.MakeGridParameters(origin_x=lon_min,
                                                         origin_y=lat_min,
                                                         upper_right_x=lon_max,
                                                         upper_right_y=lat_max,
                                                         block_size_x=dxy,
                                                         block_size_y=dxy)
    mk = meshkernel.MeshKernel(projection=projection)
    mk.curvilinear_compute_rectangular_grid_on_extension(make_grid_parameters)
    mk.curvilinear_convert_to_mesh2d() #convert to ugrid/mesh2d
    
    # refine with polygon
    pol_x = np.array([-5,-4,0,-5], dtype=np.double)
    pol_y = np.array([49,51,49.5,49], dtype=np.double)
    geometry_list = meshkernel.GeometryList(pol_x, pol_y)
    mrp = meshkernel.MeshRefinementParameters(min_edge_size=3000)
    mk.mesh2d_refine_based_on_polygon(polygon=geometry_list, mesh_refinement_params=mrp)
    
    #convert to xugrid and write to netcdf
    xu_grid_uds = dfmt.meshkernel_to_UgridDataset(mk=mk, crs=crs)
    netfile = tmp_path / 'test_startindex_net.nc'
    xu_grid_uds.ugrid.to_netcdf(netfile)
    
    # plot
    import matplotlib.pyplot as plt
    plt.close("all")
    fig,ax = plt.subplots()
    xu_grid_uds.grid.plot(ax=ax)
    ax.plot(pol_x,pol_y,'r-')
    
    #assert output grid
    ds_out = xr.open_dataset(netfile,decode_cf=False).load()
    assert ds_out.mesh2d_face_nodes.attrs['_FillValue'] == -1
    assert ds_out.mesh2d_face_nodes.attrs['start_index'] == 1
    assert 0 not in ds_out.mesh2d_face_nodes.to_numpy()
    assert ds_out.mesh2d_face_nodes.to_numpy().min() == -1
    assert ds_out.mesh2d_face_nodes.to_numpy().max() == 626
    assert "wgs84" in ds_out.variables


@pytest.mark.unittest
def test_generate_bndpli_cutland():
    # domain, resolution and expected values
    params_all = [[-68.55, -68.05, 11.95, 12.4,    1, 99],
                  [-68.55, -68.25, 11.95, 12.4,    1, 68],
                  [-68.55, -68.31, 11.95, 12.4,    2, 67],
                  [-68.31, -68.27, 12.10, 12.21,   2,  5]]
    
    dxy = 0.02
    
    for params in params_all:
        lon_min, lon_max, lat_min, lat_max, len_gdf, len_linestr0 = params
    
        mk_object = dfmt.make_basegrid(lon_min, lon_max, lat_min, lat_max, dx=dxy, dy=dxy, crs=4326)
        bnd_gdf = dfmt.generate_bndpli_cutland(mk=mk_object, res='h', buffer=0.01)
        
        # fig, ax = plt.subplots()
        # mk_object.mesh2d_get().plot_edges(ax=ax)
        # dfmt.plot_coastlines(ax=ax, crs=crs)
        # bnd_gdf.plot(ax=ax, color='r',marker='o',markersize=3)
        
        assert len(bnd_gdf) == len_gdf
        assert len(bnd_gdf.geometry[0].xy[0]) == len_linestr0

