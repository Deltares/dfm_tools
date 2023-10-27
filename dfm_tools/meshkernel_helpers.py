# -*- coding: utf-8 -*-
"""
Created on Thu Apr  6 18:46:35 2023

@author: veenstra
"""

import warnings
import xugrid as xu
import meshkernel
import xarray as xr
import datetime as dt
from dfm_tools import __version__
import getpass
import numpy as np
from dfm_tools.coastlines import get_coastlines_gdb
from netCDF4 import default_fillvals
import geopandas as gpd
from shapely import MultiPolygon, LineString, MultiLineString
from shapely.ops import linemerge


def meshkernel_delete_withcoastlines(mk:meshkernel.meshkernel.MeshKernel, res:str='f', min_area:float = 0, crs:(int,str) = None):
    """
    Wrapper around meshkernel_delete_withgdf, which automatically gets the bbox from the meshkernel object and retrieves the coastlines_gdf.

    Parameters
    ----------
    mk : meshkernel.meshkernel.MeshKernel
        DESCRIPTION.
    res : str, optional
        DESCRIPTION. The default is 'f'.
    min_area : float, optional
        DESCRIPTION. The default is 0.
    crs : (int,str), optional
        DESCRIPTION. The default is None.

    Returns
    -------
    None.

    """
    
    mesh_bnds = mk.mesh2d_get_mesh_boundaries_as_polygons()
    
    bbox = (mesh_bnds.x_coordinates.min(), mesh_bnds.y_coordinates.min(), mesh_bnds.x_coordinates.max(), mesh_bnds.y_coordinates.max())
    
    coastlines_gdf = get_coastlines_gdb(bbox=bbox, res=res, min_area=min_area, crs=crs)
    
    meshkernel_delete_withgdf(mk, coastlines_gdf)


def meshkernel_delete_withshp(mk:meshkernel.meshkernel.MeshKernel, coastlines_shp:str, crs:(int,str) = None):
    """
    Delete parts of mesh that are inside the shapefile polygon.

    Parameters
    ----------
    mk : meshkernel.meshkernel.MeshKernel
        DESCRIPTION.
    coastlines_shp : str
        Path to the shp file.
    crs : (int,str), optional
        DESCRIPTION. The default is None.

    Returns
    -------
    None.
    
    """
        
    mesh_bnds = mk.mesh2d_get_mesh_boundaries_as_polygons()
    bbox = (mesh_bnds.x_coordinates.min(), mesh_bnds.y_coordinates.min(), mesh_bnds.x_coordinates.max(), mesh_bnds.y_coordinates.max())
    coastlines_gdb = gpd.read_file(coastlines_shp, bbox=bbox)
    
    if crs:
        coastlines_gdb = coastlines_gdb.to_crs(crs)
    
    meshkernel_delete_withgdf(mk, coastlines_gdb)


def meshkernel_delete_withgdf(mk:meshkernel.meshkernel.MeshKernel, coastlines_gdf:gpd.GeoDataFrame):
    """
    Delete parts of mesh that are inside the polygons/Linestrings in a GeoDataFrame.

    Parameters
    ----------
    mk : meshkernel.meshkernel.MeshKernel
        DESCRIPTION.
    coastlines_gdf : gpd.GeoDataFrame
        DESCRIPTION.

    Returns
    -------
    None.

    """
    
    for coastline_geom in coastlines_gdf['geometry']:
        xx, yy = coastline_geom.exterior.coords.xy
        xx = np.array(xx)
        yy = np.array(yy)
        
        delete_pol_geom = meshkernel.GeometryList(x_coordinates=xx, y_coordinates=yy) #TODO: .copy()/to_numpy() makes the array contiguous in memory, which is necessary for meshkernel.mesh2d_delete()
        mk.mesh2d_delete(geometry_list=delete_pol_geom, 
                         delete_option=meshkernel.DeleteMeshOption(2), #ALL_COMPLETE_FACES/2: Delete all faces of which the complete face is inside the polygon
                         invert_deletion=False)


def meshkernel_check_geographic(mk:meshkernel.meshkernel.MeshKernel) -> bool:
    """
    Get projection from meshkernel instance

    Parameters
    ----------
    mk : meshkernel.meshkernel.MeshKernel
        DESCRIPTION.

    Returns
    -------
    bool
        DESCRIPTION.

    """
    
    if mk.get_projection()==meshkernel.ProjectionType.SPHERICAL:
        is_geographic = True
    else:
        is_geographic = False
    return is_geographic


def meshkernel_to_UgridDataset(mk:meshkernel.MeshKernel, crs:(int,str) = None, remove_noncontiguous:bool = False) -> xu.UgridDataset:
    """
    Convert a meshkernel object to a UgridDataset, including a variable with the crs (used by dflowfm to distinguish spherical/cartesian networks).
    The UgridDataset enables bathymetry interpolation and writing to netfile.

    Parameters
    ----------
    mk : meshkernel.MeshKernel
        DESCRIPTION.
    crs : (int,str), optional
        DESCRIPTION. The default is None.
    remove_noncontiguous : bool, optional
        DESCRIPTION. The default is False.

    Returns
    -------
    TYPE
        DESCRIPTION.

    """
    
    is_geographic = meshkernel_check_geographic(mk)
    
    mesh2d_grid = mk.mesh2d_get()
    
    xu_grid = xu.Ugrid2d.from_meshkernel(mesh2d_grid)
    
    #remove non-contiguous grid parts
    def xugrid_remove_noncontiguous(grid):
        #based on https://deltares.github.io/xugrid/examples/connectivity.html#connected-components
        #uses https://docs.scipy.org/doc/scipy/reference/sparse.csgraph.html
        #TODO: maybe replace with meshkernel?
        uda = xu.UgridDataArray(
            xr.DataArray(np.ones(grid.node_face_connectivity.shape[0]), dims=["face"]), grid
        )
        labels = uda.ugrid.connected_components()
        counts = labels.groupby(labels).count()
        most_frequent_label = counts["group"][np.argmax(counts.data)].item() #find largest contiguous part
        labels = labels.where(labels == most_frequent_label, drop=True)
        grid = labels.grid
        return grid
    if remove_noncontiguous:
        xu_grid = xugrid_remove_noncontiguous(xu_grid)
    
    #convert 0-based to 1-based indices for connectivity variables like face_node_connectivity
    xu_grid_ds = xu_grid.to_dataset()
    xu_grid_ds = xr.decode_cf(xu_grid_ds) #decode_cf is essential since it replaces fillvalues with nans
    ds_idx = xu_grid_ds.filter_by_attrs(start_index=0)
    for varn_conn in ds_idx.data_vars:
        xu_grid_ds[varn_conn] += 1 #from startindex 0 to 1 (fillvalues are now nans)
        xu_grid_ds[varn_conn].attrs["start_index"] += 1
        xu_grid_ds[varn_conn].encoding["_FillValue"] = -1 #can be any value <=0, but not 0 is currently the most convenient for proper xugrid plots.
    
    # convert to uds and add attrs and crs
    xu_grid_uds = xu.UgridDataset(xu_grid_ds)
    
    xu_grid_uds = xu_grid_uds.assign_attrs({#'Conventions': 'CF-1.8 UGRID-1.0 Deltares-0.10', #TODO: conventions come from xugrid, so this line is probably not necessary
                                          'institution': 'Deltares',
                                          'references': 'https://www.deltares.nl',
                                          'source': f'Created with meshkernel {meshkernel.__version__}, xugrid {xu.__version__} and dfm_tools {__version__}',
                                          'history': 'Created on %s, %s'%(dt.datetime.now().strftime('%Y-%m-%dT%H:%M:%S%z'),getpass.getuser()), #TODO: add timezone
                                          })
    #TODO: xugrid overwrites these global attributes upon saving the network file: https://github.com/Deltares/xugrid/issues/111
    
    add_crs_to_dataset(uds=xu_grid_uds,is_geographic=is_geographic,crs=crs)
    
    return xu_grid_uds


def add_crs_to_dataset(uds:(xu.UgridDataset,xr.Dataset),is_geographic:bool,crs:(str,int)):
    """
    

    Parameters
    ----------
    uds : (xu.UgridDataset,xr.Dataset)
        DESCRIPTION.
    is_geographic : bool
        whether it is a spherical (True) or cartesian (False), property comes from meshkernel instance.
    crs : (str,int)
        epsg, e.g. 'EPSG:4326' or 4326.

    Raises
    ------
    ValueError
        DESCRIPTION.

    Returns
    -------
    None.

    """
    #get crs information (name/num)
    if crs is None:
        crs_num = 0
        crs_name = ''
    else:
        crs_info = gpd.GeoSeries(crs=crs).crs #also contains area-of-use (name/bounds), datum (ellipsoid/prime-meridian)
        crs_num = crs_info.to_epsg()
        crs_name = crs_info.name
    crs_str = f'EPSG:{crs_num}'
    
    #check if combination of is_geographic and crs makes sense
    if is_geographic and crs_num!=4326:
        raise ValueError(f'provided grid is sperical (is_geographic=True) but crs="{crs}" while only "EPSG:4326" (WGS84) is supported for spherical grids') #TODO: is this true?
    if not is_geographic and crs_num==4326:
        raise ValueError('provided grid is cartesian (is_geographic=False) but crs="EPSG:4326" (WGS84), this combination is not supported')
    
    if is_geographic:
        grid_mapping_name = 'latitude_longitude'
        crs_varn = 'wgs84'
    else:
        grid_mapping_name = 'Unknown projected'
        crs_varn = 'projected_coordinate_system'
    
    attribute_dict = {
        'name': crs_name, # not required, but convenient for the user
        'epsg': np.array(crs_num, dtype=int), # epsg or EPSG_code should be present for the interacter to load the grid and by QGIS to recognize the epsg.
        'EPSG_code': crs_str, # epsg or EPSG_code should be present for the interacter to load the grid and by QGIS to recognize the epsg.
        'grid_mapping_name': grid_mapping_name, # without grid_mapping_name='latitude_longitude', interacter sees the grid as cartesian
        }
    
    uds[crs_varn] = xr.DataArray(np.array(default_fillvals['i4'],dtype=int),dims=(),attrs=attribute_dict)


def make_basegrid(lon_min,lon_max,lat_min,lat_max,dx,dy,angle=0,is_geographic=True):
    """
    empty docstring
    """
    # create base grid
    make_grid_parameters = meshkernel.MakeGridParameters(angle=angle,
                                                         origin_x=lon_min,
                                                         origin_y=lat_min,
                                                         upper_right_x=lon_max,
                                                         upper_right_y=lat_max,
                                                         block_size_x=dx,
                                                         block_size_y=dy)
    
    mk = meshkernel.MeshKernel(is_geographic=is_geographic)
    mk.curvilinear_make_uniform_on_extension(make_grid_parameters)
    mk.curvilinear_convert_to_mesh2d() #convert to ugrid/mesh2d
    
    return mk


def refine_basegrid(mk, data_bathy_sel, min_edge_size):
    """
    empty docstring
    """
    
    lon_np = data_bathy_sel.lon.to_numpy()
    lat_np = data_bathy_sel.lat.to_numpy()
    values_np = data_bathy_sel.elevation.to_numpy().flatten().astype('float') #TODO: astype to avoid "TypeError: incompatible types, c_short_Array_74880 instance instead of LP_c_double instance"
    gridded_samples = meshkernel.GriddedSamples(x_coordinates=lon_np,y_coordinates=lat_np,values=values_np)


    #refinement
    mesh_refinement_parameters = meshkernel.MeshRefinementParameters(min_edge_size=min_edge_size, #always in meters
                                                                     refinement_type=meshkernel.RefinementType(1), #Wavecourant/1,
                                                                     connect_hanging_nodes=True, #set to False to do multiple refinement steps (e.g. for multiple regions)
                                                                     smoothing_iterations=2,
                                                                     max_courant_time=120)
    
    mk.mesh2d_refine_based_on_gridded_samples(gridded_samples=gridded_samples,
                                               mesh_refinement_params=mesh_refinement_parameters,
                                               use_nodal_refinement=True) #TODO: what does this do?
    
    return mk


def generate_bndpli_cutland(mk:meshkernel.meshkernel.MeshKernel, res:str='f', min_area:float = 0, crs:(int,str) = None, buffer:float = 0):
    """
    Generate a boundary polyline from the meshkernel object and cut away the landward part.
    Be sure to do this on the base/refined grid, not on the grid where the landward cells were already cut.
    
    Parameters
    ----------
    mk : meshkernel.meshkernel.MeshKernel
        DESCRIPTION.
    res : str, optional
        DESCRIPTION. The default is 'f'.
    min_area : float, optional
        DESCRIPTION. The default is 0.
    crs : (int,str), optional
        DESCRIPTION. The default is None.
    buffer : float, optional
        DESCRIPTION. The default is 0.

    Returns
    -------
    bnd_gdf : TYPE
        DESCRIPTION.

    """
    
    mesh_bnds = mk.mesh2d_get_mesh_boundaries_as_polygons()
    if mesh_bnds.geometry_separator in mesh_bnds.x_coordinates:
        raise Exception('use dfmt.generate_bndpli_cutland() on an uncut grid')
    mesh_bnds_xy = np.c_[mesh_bnds.x_coordinates,mesh_bnds.y_coordinates]
    
    bbox = (mesh_bnds.x_coordinates.min(), mesh_bnds.y_coordinates.min(), mesh_bnds.x_coordinates.max(), mesh_bnds.y_coordinates.max())
    coastlines_gdf = get_coastlines_gdb(bbox=bbox, res=res, min_area=min_area, crs=crs)
    
    meshbnd_ls = LineString(mesh_bnds_xy)
    coastlines_mp = MultiPolygon(coastlines_gdf.geometry.tolist())
    coastlines_mp = coastlines_mp.buffer(buffer)
    bnd_ls = meshbnd_ls.difference(coastlines_mp)
    
    #attempt to merge MultiLineString to single LineString
    if isinstance(bnd_ls,MultiLineString):
        print('attemting to merge lines in MultiLineString to single LineString (if connected)')
        bnd_ls = linemerge(bnd_ls)
    
    #convert MultiLineString/LineString to GeoDataFrame
    if isinstance(bnd_ls,MultiLineString):
        bnd_gdf = gpd.GeoDataFrame(geometry=list(bnd_ls.geoms))
    elif isinstance(bnd_ls,LineString):
        bnd_gdf = gpd.GeoDataFrame(geometry=[bnd_ls])
    
    #set crs from coastlines
    bnd_gdf.crs = coastlines_gdf.crs
    return bnd_gdf


def interpolate_bndpli(bnd_gdf,res):
    """
    interpolate bnd_gdf to a new resolution
    """
    #TODO: keep corners of grid, maybe use mk.polygon_refine()
    
    bnd_gdf_interp = bnd_gdf.copy()
    for irow,row in bnd_gdf_interp.iterrows():
        bnd_ls = row.geometry
        interp_range = np.arange(0,bnd_ls.length,res)
        bnd_ls_interp_points = bnd_ls.interpolate(interp_range)
        if len(bnd_ls_interp_points)==1: #no change if interp results in only one point
            continue
        bnd_ls_interp = LineString(bnd_ls_interp_points)
        bnd_gdf_interp['geometry'][irow] = bnd_ls_interp
    return bnd_gdf_interp

