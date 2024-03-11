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

__all__ = [
    "meshkernel_delete_withcoastlines",
    "meshkernel_delete_withshp",
    "meshkernel_delete_withgdf",
    "meshkernel_to_UgridDataset",
    "make_basegrid",
    "refine_basegrid",
    "generate_bndpli_cutland",
    "interpolate_bndpli",
    ]

def meshkernel_delete_withcoastlines(mk:meshkernel.MeshKernel, res:str='f', min_area:float = 0, crs:(int,str) = None):
    """
    Wrapper around meshkernel_delete_withgdf, which automatically gets the bbox from the meshkernel object and retrieves the coastlines_gdf.

    Parameters
    ----------
    mk : meshkernel.MeshKernel
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


def meshkernel_delete_withshp(mk:meshkernel.MeshKernel, coastlines_shp:str, crs:(int,str) = None):
    """
    Delete parts of mesh that are inside the shapefile polygon.

    Parameters
    ----------
    mk : meshkernel.MeshKernel
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


def meshkernel_delete_withgdf(mk:meshkernel.MeshKernel, coastlines_gdf:gpd.GeoDataFrame):
    """
    Delete parts of mesh that are inside the polygons/Linestrings in a GeoDataFrame.

    Parameters
    ----------
    mk : meshkernel.MeshKernel
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
                         delete_option=meshkernel.DeleteMeshOption.INSIDE_NOT_INTERSECTED,
                         invert_deletion=False)


def geographic_to_meshkernel_projection(is_geographic:bool) -> meshkernel.ProjectionType:
    """
    converts is_geographic boolean to meshkernel.ProjectionType (SPHERICAL OR CARTESIAN)

    Parameters
    ----------
    is_geographic : bool
        DESCRIPTION.

    Returns
    -------
    projection : TYPE
        DESCRIPTION.

    """
    if is_geographic:
        projection = meshkernel.ProjectionType.SPHERICAL
    else:
        projection = meshkernel.ProjectionType.CARTESIAN
    return projection


def meshkernel_is_geographic(mk):
    if mk.get_projection()==meshkernel.ProjectionType.CARTESIAN:
        is_geographic = False
    else:
        is_geographic = True
    return is_geographic


def meshkernel_to_UgridDataset(mk:meshkernel.MeshKernel, crs:(int,str) = None) -> xu.UgridDataset:
    """
    Convert a meshkernel object to a UgridDataset, including a variable with the crs (used by dflowfm to distinguish spherical/cartesian networks).
    The UgridDataset enables bathymetry interpolation and writing to netfile.

    Parameters
    ----------
    mk : meshkernel.MeshKernel
        DESCRIPTION.
    crs : (int,str), optional
        DESCRIPTION. The default is None.

    Returns
    -------
    TYPE
        DESCRIPTION.

    """
    
    crs_is_geographic = crs_to_isgeographic(crs)
    
    mesh2d_grid = mk.mesh2d_get()
    
    #check if both crs and grid are geograpic or not
    #TODO: do this in xugrid: https://github.com/Deltares/xugrid/issues/188
    grid_is_geographic = meshkernel_is_geographic(mk)
    if crs_is_geographic != grid_is_geographic:
        raise ValueError(f"crs has is_geographic={crs_is_geographic} and grid has is_geographic={grid_is_geographic}. This is not allowed.")
    
    # TODO: below is not correctly handled by xugrid yet, projected=False does not give is_geographic=True
    # related issue is https://github.com/Deltares/dfm_tools/issues/686
    xu_grid = xu.Ugrid2d.from_meshkernel(mesh2d_grid, projected= not crs_is_geographic, crs=crs)
    
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
    
    # add crs including attrs
    if crs is not None:
        xu_grid_uds.ugrid.set_crs(crs)
        uds_add_crs_attrs(xu_grid_uds)
    return xu_grid_uds


def uds_get_crs(uds):
    uds_crs_dict = uds.ugrid.crs
    if uds_crs_dict is None:
        return None
    else:
        keys = list(uds_crs_dict.keys())
        crs = uds_crs_dict[keys[0]]
        return crs


def uds_add_crs_attrs(uds:(xu.UgridDataset,xr.Dataset)):
    """
    

    Parameters
    ----------
    uds : (xu.UgridDataset,xr.Dataset)
        DESCRIPTION.
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
    # TODO: upon crs conversion and to_netcdf(), the netcdf file could get two crs vars, catch this
    
    grid_is_geographic = uds.grid.is_geographic
    
    crs = uds_get_crs(uds)
    
    #get crs information (name/num)
    if crs is None:
        crs_num = 0
        crs_str = 'EPSG:0'
        crs_name = ''
        crs_is_geographic = False
    else:
        # get crs info, should actually be `import pyproj; pyproj.CRS.from_user_input(crs)`
        crs_info = gpd.GeoSeries(crs=crs).crs #also contains area-of-use (name/bounds), datum (ellipsoid/prime-meridian)
        crs_num = crs_info.to_epsg()
        crs_str = crs_info.to_string()
        crs_name = crs_info.name
        crs_is_geographic = crs_info.is_geographic
        # TODO: standard names for lat/lon in crs_info.cs_to_cf()
    
    #check if combination of is_geographic and crs makes sense
    if grid_is_geographic != crs_is_geographic:
        raise ValueError(f"`grid_is_geographic` mismatch between provided grid (is_geographic={grid_is_geographic}) and provided crs ({crs}, is_geographic={crs_is_geographic})")
    
    # TODO: consider always using the same crs_varn, align with xugrid
    # QGIS also does not recognize epsg anymore when renaming variable to `crs` 
    # or something else (`wgs84` and `projected_coordinate_system` both do work)
    if grid_is_geographic:
        crs_varn = 'wgs84'
    else:
        crs_varn = 'projected_coordinate_system'
    
    attribute_dict = {
        'name': crs_name, # not required, but convenient for the user
        'epsg': np.array(crs_num, dtype=int), # epsg or EPSG_code should be present for the interacter to load the grid and by QGIS to recognize the epsg.
        'EPSG_code': crs_str, # epsg or EPSG_code should be present for the interacter to load the grid and by QGIS to recognize the epsg.
        }
    if grid_is_geographic:
        # without grid_mapping_name="latitude_longitude", interacter sees the grid as cartesian
        # grid_mapping_name is not available for projected crs's
        attribute_dict['grid_mapping_name'] = crs_info.to_cf()['grid_mapping_name']
    
    uds[crs_varn] = xr.DataArray(np.array(default_fillvals['i4'],dtype=int),dims=(),attrs=attribute_dict)


def crs_to_isgeographic(crs=None):
    if crs is None:
        is_geographic = False
    else:
        import pyproj
        crs = pyproj.CRS.from_user_input(crs)
        is_geographic = crs.is_geographic
    return is_geographic


def make_basegrid(lon_min,lon_max,lat_min,lat_max,dx,dy,angle=0,
                  crs=None, is_geographic=None):
    """
    empty docstring
    """
    
    if is_geographic is not None:
        raise ValueError("'is_geographic' was deprecated as argument for dfmt.make_basegrid(), it is now derived from 'crs' instead")
    
    
    is_geographic = crs_to_isgeographic(crs)
    projection = geographic_to_meshkernel_projection(is_geographic)
    
    # create base grid
    make_grid_parameters = meshkernel.MakeGridParameters(angle=angle,
                                                         origin_x=lon_min,
                                                         origin_y=lat_min,
                                                         upper_right_x=lon_max,
                                                         upper_right_y=lat_max,
                                                         block_size_x=dx,
                                                         block_size_y=dy)
    
    mk = meshkernel.MeshKernel(projection=projection)
    mk.curvilinear_compute_rectangular_grid_on_extension(make_grid_parameters)
    mk.curvilinear_convert_to_mesh2d() #convert to ugrid/mesh2d
    
    return mk


def refine_basegrid(mk, data_bathy_sel, min_edge_size):
    """
    empty docstring
    """
    
    lon_np = data_bathy_sel.lon.to_numpy()
    lat_np = data_bathy_sel.lat.to_numpy()
    values_np = data_bathy_sel.to_numpy().flatten()
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


def generate_bndpli_cutland(mk:meshkernel.MeshKernel, res:str='f', min_area:float = 0, crs:(int,str) = None, buffer:float = 0):
    """
    Generate a boundary polyline from the meshkernel object and cut away the landward part.
    Be sure to do this on the base/refined grid, not on the grid where the landward cells were already cut.
    
    Parameters
    ----------
    mk : meshkernel.MeshKernel
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

