import os
import numpy as np
import xugrid as xu
import xarray as xr
import datetime as dt
import pandas as pd
import meshkernel
from dfm_tools.xarray_helpers import file_to_list
from netCDF4 import default_fillvals
import warnings

__all__ = [
    "open_partitioned_dataset",
    "open_dataset_curvilinear",
    "open_dataset_delft3d4",
    "uda_to_faces",
    "uda_interfaces_to_centers",
    "add_network_cellinfo",
    "enrich_rst_with_map",
]


def get_vertical_dimensions(uds): #TODO: maybe add layer_dimension and interface_dimension properties to xugrid?
    """
    get vertical_dimensions from grid_info of ugrid mapfile (this will fail for hisfiles). The info is stored in the layer_dimension and interface_dimension attribute of the mesh2d variable of the dataset (stored in uds.grid after reading with xugrid)
    
    processing cb_3d_map.nc
        >> found layer/interface dimensions in file: mesh2d_nLayers mesh2d_nInterfaces
    processing Grevelingen-FM_0*_map.nc
        >> found layer/interface dimensions in file: nmesh2d_layer nmesh2d_interface (these are updated in open_partitioned_dataset)
    processing DCSM-FM_0_5nm_0*_map.nc
        >> found layer/interface dimensions in file: mesh2d_nLayers mesh2d_nInterfaces
    processing MB_02_0*_map.nc
        >> found layer/interface dimensions in file: mesh2d_nLayers mesh2d_nInterfaces
    """
        
    if not hasattr(uds,'grid'): #early return in case of e.g. hisfile
        return None, None
        
    gridname = uds.grid.name
    grid_info = uds.grid.to_dataset()[gridname]
    if hasattr(grid_info,'layer_dimension'):
        return grid_info.layer_dimension, grid_info.interface_dimension
    else: # 2D model or networkfile
        return None, None


def remove_ghostcells(uds, fname): #TODO: remove ghostcells from output or align values between (non)ghost cells: https://issuetracker.deltares.nl/browse/UNST-6701
    """
    Dropping ghostcells if there is a domainno variable present and there is a domainno in the filename.
    Not using most-occurring domainno in var, since this is not a valid assumption for merged datasets and might be invalid for a very small partition.
    
    """
    gridname = uds.grid.name
    varn_domain = f'{gridname}_flowelem_domain'
    
    #check if dataset has domainno variable, return uds if not present
    if varn_domain not in uds.data_vars:
        print('[nodomainvar] ',end='')
        return uds
    
    #derive domainno from filename, return uds if not present
    if '_' not in fname: #safety escape in case there is no _ in the filename
        print('[nodomainfname] ',end='')
        return uds
    fname_splitted = fname.split('_')
    part_domainno_fromfname = fname_splitted[-2] #this is not valid for rstfiles (date follows after partnumber), but they cannot be read with xugrid anyway since they lack topology and node_x/node_y variables: https://issuetracker.deltares.nl/browse/UNST-7176
    if not part_domainno_fromfname.isnumeric() or len(part_domainno_fromfname)!=4:
        print('[nodomainfname] ',end='')
        return uds
    
    #drop ghostcells
    part_domainno_fromfname = int(part_domainno_fromfname)
    da_domainno = uds.variables[varn_domain]
    idx = np.flatnonzero(da_domainno == part_domainno_fromfname)
    uds = uds.isel({uds.grid.face_dimension:idx})
    return uds


def remove_unassociated_edges(ds: xr.Dataset) -> xr.Dataset:
    """
    Removes edges that are not associated to any of the faces, usecase in https://github.com/Deltares/xugrid/issues/68

    Parameters
    ----------
    ds : xr.Dataset
        DESCRIPTION.
    topology : str, optional
        DESCRIPTION. The default is None.

    Returns
    -------
    ds : xr.Dataset
        DESCRIPTION.

    """
    
    uds = xu.core.wrap.UgridDataset(ds)
    
    # escape for 1D networks
    if isinstance(uds.grid,xu.Ugrid1d):
        return ds
    
    associated = uds.grid.validate_edge_node_connectivity()
    edge_dimension = uds.grid.edge_dimension
    
    if not associated.all():
         print(f"[{(~associated).sum()} unassociated edges removed] ",end='')
         ds = ds.sel({edge_dimension: associated})
    return ds


def decode_default_fillvals(ds):
    """
    xarray only supports explicitly set _FillValue attrs, and therefore ignores the default netCDF4 fillvalue
    This function adds the default fillvalue as _FillValue attribute and decodes the dataset again.
    """
    # TODO: this function can be removed when xarray does it automatically: https://github.com/Deltares/dfm_tools/issues/490
    
    nfillattrs_added = 0
    for varn in ds.variables.keys():
        # TODO: possible to get always_mask boolean with `netCDF4.Dataset(file_nc).variables[varn].always_mask`, but this seems to be always True for FM mapfiles
        if '_FillValue' in ds.variables[varn].encoding:
            continue
        dtype_str = ds.variables[varn].dtype.str[1:]
        if dtype_str not in default_fillvals.keys():
            continue
        varn_fillval = default_fillvals[dtype_str]
        ds[varn] = ds[varn].assign_attrs({'_FillValue':varn_fillval})
        nfillattrs_added += 1
    print(f'[default_fillvals decoded for {nfillattrs_added} variables] ',end='')
    
    #decode the dataset with newly added _FillValue attrs again
    ds = xr.decode_cf(ds)
    return ds


def remove_nan_fillvalue_attrs(ds : (xr.Dataset, xu.UgridDataset)):
    """
    xarray writes {"_FillValue": np.nan} to encoding for variables without _FillValue attribute.
    Remove these again upon reading to avoid issues.
    """
    if isinstance(ds,xu.UgridDataset):
        ds = ds.obj
    
    count = 0
    for varn in ds.variables.keys():
        if '_FillValue' in ds.variables[varn].encoding:
            if np.isnan(ds.variables[varn].encoding['_FillValue']):
                ds[varn].encoding.pop('_FillValue')
                count += 1
    if count > 0:
        print(f"[{count} nan fillvalue attrs removed]", end="")


def uds_auto_set_crs(uds : xu.UgridDataset):
    # FM-mapfiles contain wgs84/projected_coordinate_system variables with epsg attr, xugrid has .crs property
    # TODO: parse+set crs in xugrid instead: https://github.com/Deltares/xugrid/issues/42
    # also adjusting projected_coordinate_system/wgs84 when using set_crs/to_crs
    
    uds_epsg = uds.filter_by_attrs(epsg=lambda v: v is not None)
    if len(uds_epsg.data_vars) != 1:
        return
    
    crs_varn = list(uds_epsg.data_vars)[0]
    epsg = uds[crs_varn].attrs["epsg"]
    from pyproj.exceptions import CRSError
    try:
        uds.ugrid.set_crs(epsg)
    except CRSError:
        return


def open_partitioned_dataset(file_nc:str, decode_fillvals:bool = False, remove_edges:bool = False, remove_ghost:bool = True, **kwargs): 
    """
    using xugrid to read and merge partitions, including support for delft3dfm mapformat1 
    by renaming old layerdim. Furthermore some optional extensions like removal of hanging
    edges and removal of ghost cells.

    Parameters
    ----------
    file_nc : str
        DESCRIPTION.
    decode_fillvals : bool, optional
        DESCRIPTION. The default is False.
    remove_edges : bool, optional
        Remove hanging edges from the mapfile, necessary to generate contour and contourf plots 
        with xugrid. Can be enabled always, but takes some additional computation time. The default is False.
    remove_ghost : bool, optional
        Remove ghostcells from the partitions. This is also done by xugrid automatically 
        upon merging, but then the domain numbers are not taken into account so 
        the result will be different. The default is True.
    file_nc : TYPE
        DESCRIPTION.
    kwargs : TYPE, optional
        arguments that are passed to xr.open_dataset. The chunks argument is set if not provided
        chunks={'time':1} increases performance significantly upon reading, but causes memory overloads when performing sum/mean/etc actions over time dimension (in that case 100/200 is better). The default is {'time':1}.

    Raises
    ------
    Exception
        DESCRIPTION.

    Returns
    -------
    ds_merged_xu : TYPE
        DESCRIPTION.
    
    """
    # TODO: timings
    # file_nc = 'p:\\1204257-dcsmzuno\\2006-2012\\3D-DCSM-FM\\A18b_ntsu1\\DFM_OUTPUT_DCSM-FM_0_5nm\\DCSM-FM_0_5nm_0*_map.nc' #3D DCSM
    # file_nc = 'p:\\archivedprojects\\11206813-006-kpp2021_rmm-2d\\C_Work\\31_RMM_FMmodel\\computations\\model_setup\\run_207\\results\\RMM_dflowfm_0*_map.nc' #RMM 2D
    # file_nc = 'p:\\1230882-emodnet_hrsm\\GTSMv5.0\\runs\\reference_GTSMv4.1_wiCA_2.20.06_mapformat4\\output\\gtsm_model_0*_map.nc' #GTSM 2D
    # file_nc = 'p:\\11208053-005-kpp2022-rmm3d\\C_Work\\01_saltiMarlein\\RMM_2019_computations_02\\computations\\theo_03\\DFM_OUTPUT_RMM_dflowfm_2019\\RMM_dflowfm_2019_0*_map.nc' #RMM 3D
    # file_nc = 'p:\\archivedprojects\\11203379-005-mwra-updated-bem\\03_model\\02_final\\A72_ntsu0_kzlb2\\DFM_OUTPUT_MB_02\\MB_02_0*_map.nc'
    # Timings (xu.open_dataset/xu.merge_partitions):
    # - DCSM 3D 20 partitions  367 timesteps: 231.5/ 4.5 sec (decode_times=False: 229.0 sec)
    # - RMM  2D  8 partitions  421 timesteps:  55.4/ 4.4 sec (decode_times=False:  56.6 sec)
    # - GTSM 2D  8 partitions  746 timesteps:  71.8/30.0 sec (decode_times=False: 204.8 sec)
    # - RMM  3D 40 partitions  146 timesteps: 168.8/ 6.3 sec (decode_times=False: 158.4 sec)
    # - MWRA 3D 20 partitions 2551 timesteps:  74.4/ 3.4 sec (decode_times=False:  79.0 sec)

    #TODO: add support for multiple grids via keyword? https://github.com/Deltares/dfm_tools/issues/497
    #TODO: speed up open_dataset https://github.com/Deltares/dfm_tools/issues/225 (also remove_ghost)
    
    if 'chunks' not in kwargs:
        kwargs['chunks'] = {'time':1}
    if 'decode_timedelta' not in kwargs:
        # avoid futurewarning: https://github.com/Deltares/dfm_tools/issues/1100
        kwargs['decode_timedelta'] = False
    
    dtstart_all = dt.datetime.now()
    file_nc_list = file_to_list(file_nc)
    
    print(f'>> xu.open_dataset() with {len(file_nc_list)} partition(s): ',end='')
    dtstart = dt.datetime.now()
    partitions = []
    for iF, file_nc_one in enumerate(file_nc_list):
        print(iF+1,end=' ')
        # suppress chunking warning: https://github.com/Deltares/dfm_tools/issues/947
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=UserWarning)
            ds = xr.open_mfdataset(file_nc_one, **kwargs)
        if decode_fillvals:
            ds = decode_default_fillvals(ds)
        if remove_edges:
            ds = remove_unassociated_edges(ds)
        if 'nFlowElem' in ds.dims and 'nNetElem' in ds.dims:
            print('[mapformat1] ',end='')
            #for mapformat1 mapfiles: merge different face dimensions (rename nFlowElem to nNetElem) to make sure the dataset topology is correct
            ds = ds.rename({'nFlowElem':'nNetElem'})
        remove_nan_fillvalue_attrs(ds)
        uds = xu.core.wrap.UgridDataset(ds)
        if remove_ghost: #TODO: this makes it way slower (at least for GTSM, although merging seems faster), but is necessary since values on overlapping cells are not always identical (eg in case of Venice ucmag)
            uds = remove_ghostcells(uds, file_nc_one)
        uds_auto_set_crs(uds)
        partitions.append(uds)
    print(': ',end='')
    print(f'{(dt.datetime.now()-dtstart).total_seconds():.2f} sec')
    
    if len(partitions) == 1: #do not merge in case of 1 partition
        return partitions[0]
    
    print(f'>> xu.merge_partitions() with {len(file_nc_list)} partition(s): ',end='')
    dtstart = dt.datetime.now()
    ds_merged_xu = xu.merge_partitions(partitions)
    print(f'{(dt.datetime.now()-dtstart).total_seconds():.2f} sec')
    
    #print variables that are dropped in merging procedure. Often only ['mesh2d_face_x_bnd', 'mesh2d_face_y_bnd'], which can be derived by combining node_coordinates (mesh2d_node_x mesh2d_node_y) and face_node_connectivity (mesh2d_face_nodes). >> can be removed from FM-mapfiles (email of 16-1-2023)
    varlist_onepart = list(partitions[0].variables.keys())
    varlist_merged = list(ds_merged_xu.variables.keys())
    varlist_dropped_bool = ~pd.Series(varlist_onepart).isin(varlist_merged)
    varlist_dropped = pd.Series(varlist_onepart).loc[varlist_dropped_bool]
    if varlist_dropped_bool.any():
        print(f'>> some variables dropped with merging of partitions: {varlist_dropped.tolist()}')
    
    print(f'>> dfmt.open_partitioned_dataset() total: {(dt.datetime.now()-dtstart_all).total_seconds():.2f} sec')
    return ds_merged_xu


def open_dataset_curvilinear(file_nc,
                             x_dim:str,
                             y_dim:str,
                             x_bounds:str,
                             y_bounds:str,
                             convert_360to180:bool = False ,
                             **kwargs) -> xu.UgridDataset:
    """
    Construct a UgridDataset from a curvilinear grid with 2D lat/lon variables
    with i/j indexes/dims, including vertices variables. Works for curvilinear
    datasets like CMCC and also for WAQUA files that are converted with getdata

    Parameters
    ----------
    file_nc : str or path
        DESCRIPTION.
    x_dim : str
        The x-dimension, like lon, i or N.
    y_dim : str
        The y-dimension, like lat, j or M.
    x_bounds : str
        The variable with the x-bounds, like vertices_longitude or grid_x.
    y_bounds : str
        The variable with the y-bounds, like vertices_latitude or grid_y.
    convert_360to180 : bool, optional
        Whether to convert from a 0 to 360 degree global model to a -180 to 180
        degree global model. The default is False.
    **kwargs : TYPE
        additional arguments are passed on to xr.open_mfdataset().

    Returns
    -------
    uds : xu.UgridDataset
        The resulting ugrid dataset.

    """
    
    if 'chunks' not in kwargs:
        kwargs['chunks'] = {'time':1}
    
    # data_vars='minimal' to avoid time dimension on vertices_latitude and
    # others when opening multiple files at once
    print('>> open_mfdataset: ',end='')
    dtstart = dt.datetime.now()
    ds = xr.open_mfdataset(file_nc, data_vars="minimal", **kwargs)
    print(f'{(dt.datetime.now()-dtstart).total_seconds():.2f} sec')

    # convert from 0to360 to -180 to 180
    if convert_360to180:
        ds[x_bounds] = (ds[x_bounds]+180) % 360 - 180
    
    topology = {"mesh2d":{"x":x_dim,
                          "y":y_dim,
                          "x_bounds":x_bounds,
                          "y_bounds":y_bounds,
                          },
                }
    
    print('>> convert to xugrid.UgridDataset: ',end='')
    dtstart = dt.datetime.now()
    uds = xu.UgridDataset.from_structured2d(ds, topology=topology)
    print(f'{(dt.datetime.now()-dtstart).total_seconds():.2f} sec')
    return uds


def delft3d4_get_nanmask(x,y):
    # -999.999 in kivu and 0.0 in curvedbend, both in westernscheldt
    bool_0 = (x==0) & (y==0)
    bool_1 = (x==-999) & (y==-999)
    bool_2 = (x==-999.999) & (y==-999.999)
    bool_mask = bool_0 | bool_1 | bool_2
    return bool_mask


def delft3d4_stack_shifted_coords(da):
    shift = 1
    np_stacked = np.stack([
        da, #ll
        da.shift(MC=shift), #lr
        da.shift(MC=shift, NC=shift), #ur
        da.shift(NC=shift), #ul
        ],axis=-1)
    da_stacked = xr.DataArray(np_stacked, dims=("M","N","four"))
    return da_stacked


def delft3d4_convert_uv(ds):
    # replace invalid values not with nan but with zero
    # otherwise the spatial coverage is affected
    mask_u1 = (ds.U1==-999) | (ds.U1==-999.999)
    mask_v1 = (ds.V1==-999) | (ds.V1==-999.999)
    u1_mn = ds.U1.where(~mask_u1, 0)
    v1_mn = ds.V1.where(~mask_v1, 0)
    
    # minus 0.5 since padding=low so corner value is representative for
    # previous face according to that logic, method=nearest might be better
    # (or just rename the dims)
    u1_mn_cen = u1_mn.interp(MC=u1_mn.MC-0.5, method='linear')
    v1_mn_cen = v1_mn.interp(NC=v1_mn.NC-0.5, method='linear')
    # rename corner dims to center dims since we shifted them with 0.5
    u1_mn_cen = u1_mn_cen.rename({'MC':'M'})
    v1_mn_cen = v1_mn_cen.rename({'NC':'N'})
    # TODO: since padding=low, just renaming the dims might even be better?
    # u1_mn_cen = u1_mn.rename({'MC':'M'})
    # v1_mn_cen = v1_mn.rename({'NC':'N'})
    # >> could also be done with `ds = ds.swap_dims({"M":"MC","N":"NC"})`
    
    # create combined uv mask (have to rename dimensions)
    mask_u1_mn = mask_u1.rename({'MC':'M'})
    mask_v1_mn = mask_v1.rename({'NC':'N'})
    mask_uv1_mn = mask_u1_mn & mask_v1_mn
    # drop all actual missing cells
    u1_mn_cen = u1_mn_cen.where(~mask_uv1_mn)
    v1_mn_cen = v1_mn_cen.where(~mask_uv1_mn)
    
    # to avoid creating large chunks, alternative is to overwrite the vars
    # with the MN-averaged vars, but it requires passing and updating of attrs
    ds = ds.drop_vars(['U1','V1'])
    
    # compute ux/uy/umag/udir
    # TODO: add attrs to variables
    alfas_rad = np.deg2rad(ds.ALFAS)
    vel_x = u1_mn_cen*np.cos(alfas_rad) - v1_mn_cen*np.sin(alfas_rad)
    vel_y = u1_mn_cen*np.sin(alfas_rad) + v1_mn_cen*np.cos(alfas_rad)
    ds['ux'] = vel_x
    ds['uy'] = vel_y
    ds['umag'] = np.sqrt(vel_x**2 + vel_y**2)
    ds['udir'] = np.rad2deg(np.arctan2(vel_y, vel_x))%360
    return ds


def open_dataset_delft3d4(file_nc, **kwargs) -> xu.UgridDataset:
    """
    Reads in a Delft3D4 netcdf outputfile (curvilinear/staggered) as a ugrid
    (xugrid.UgridDataset) dataset. This is a Delft3D4 specific version of
    dfmt.open_dataset_curvilinear().
    
    To get Delft3D4 to write netCDF output instead of .dat files, add these
    lines to your model settings file (.mdf):
    
    - FlNcdf=#maphis#
    - ncFormat=4

    Parameters
    ----------
    file_nc : str or path
        DESCRIPTION.
    **kwargs : TYPE
        additional arguments are passed on to xr.open_mfdataset().

    Returns
    -------
    uds : xu.UgridDataset
        The resulting ugrid dataset.

    """
    
    if 'chunks' not in kwargs:
        kwargs['chunks'] = {'time':1}
    
    ds = xr.open_dataset(file_nc, **kwargs)
    
    # prevent grid variable that might be confused with uds.grid accessor
    if 'grid' in ds.data_vars:
        ds = ds.rename_vars({'grid':'grid_original'})
    
    xcor_stacked = delft3d4_stack_shifted_coords(ds.XCOR)
    ycor_stacked = delft3d4_stack_shifted_coords(ds.YCOR)
    mask_xy = delft3d4_get_nanmask(xcor_stacked,ycor_stacked)
    ds['xcor_stacked'] = xcor_stacked.where(~mask_xy)
    ds['ycor_stacked'] = ycor_stacked.where(~mask_xy)
    
    if ('U1' in ds.data_vars) and ('V1' in ds.data_vars):
        ds = delft3d4_convert_uv(ds)
    
    # TODO: consider using same dims for variables on cell corners and faces
    # ds = ds.swap_dims({"M":"MC","N":"NC"})
    topology = {"mesh2d":{"x":"M",
                          "y":"N",
                          "x_bounds":"xcor_stacked",
                          "y_bounds":"ycor_stacked",
                          }
                }
    uds = xu.UgridDataset.from_structured2d(ds, topology=topology)
    
    return uds


def uda_to_faces(uda : xu.UgridDataArray) -> xu.UgridDataArray:
    """
    Interpolates a ugrid variable (xu.DataArray) with a node or edge dimension
    to the faces by averaging the 3/4 nodes/edges around each face.
    
    Parameters
    ----------
    uda_node : xu.UgridDataArray
        DESCRIPTION.

    Raises
    ------
    KeyError
        DESCRIPTION.

    Returns
    -------
    uda_face : xu.UgridDataArray
        DESCRIPTION.

    """
    grid = uda.grid
    
    dimn_faces = grid.face_dimension
    reduce_dim = 'nMax_face_nodes' #arbitrary dimname that is reduced anyway
    dimn_nodes = grid.node_dimension
    dimn_edges = grid.edge_dimension
    
    # construct indexing array
    if dimn_nodes in uda.dims:
        dimn_notfaces_name = "node"
        dimn_notfaces = dimn_nodes
        indexer_np = grid.face_node_connectivity
    elif dimn_edges in uda.dims:
        dimn_notfaces_name = "edge"
        dimn_notfaces = dimn_edges
        indexer_np = grid.face_edge_connectivity
    else:
        print(f'provided uda/variable "{uda.name}" does not have an node or edge dimension, returning unchanged uda')
        return uda
    
    # rechunk to make sure the node/edge dimension is not chunked, otherwise we
    # will get "PerformanceWarning: Slicing with an out-of-order index is
    # generating 384539 times more chunks."
    chunks = {dimn_notfaces:-1}
    uda = uda.chunk(chunks)

    indexer = xr.DataArray(indexer_np,dims=(dimn_faces,reduce_dim))
    indexer_validbool = indexer!=-1
    indexer = indexer.where(indexer_validbool,-1)
    
    print(f'{dimn_notfaces_name}-to-face interpolation: ',end='')
    dtstart = dt.datetime.now()
    # for each face, select all corresponding node/edge values
    # we do this via stack and unstack since 2D indexing does not
    # properly work in dask yet: https://github.com/dask/dask/pull/10237
    # this process converts the xu.UgridDataArray to a xr.DataArray, so we convert it back
    indexer_stacked = indexer.stack(__tmp_dim__=(dimn_faces, reduce_dim))
    uda_face_allnodes_ds_stacked = uda.isel({dimn_notfaces: indexer_stacked})
    uda_face_allnodes_ds = uda_face_allnodes_ds_stacked.unstack("__tmp_dim__")
    uda_face_allnodes = xu.UgridDataArray(uda_face_allnodes_ds,grid=grid)
    
    # replace nonexistent nodes/edges with nan
    # replace all values for fillvalue nodes/edges (-1) with nan
    uda_face_allnodes = uda_face_allnodes.where(indexer_validbool)
    # average node/edge values per face
    uda_face = uda_face_allnodes.mean(dim=reduce_dim,keep_attrs=True)
    #update attrs from node/edge to face
    face_attrs = {'location': 'face', 'cell_methods': f'{dimn_faces}: mean'}
    uda_face = uda_face.assign_attrs(face_attrs)
    print(f'{(dt.datetime.now()-dtstart).total_seconds():.2f} sec')
    
    return uda_face


def uda_interfaces_to_centers(uda_int : xu.UgridDataArray) -> xu.UgridDataArray:
    dimn_layer, dimn_interface = get_vertical_dimensions(uda_int)
    
    if dimn_interface not in uda_int.dims:
        print('no interface dimension found, returning original array')
        return uda_int
        
    # define layers to be halfway inbetween the interfaces
    nlayers = uda_int.sizes[dimn_interface] - 1
    array_shift_half = xr.DataArray(np.arange(0.5,nlayers),dims=dimn_layer)
    # interpolate from interfaces to layers
    uda_cen = uda_int.interp({dimn_interface:array_shift_half},assume_sorted=True)
    # drop interface coordinate
    uda_cen = uda_cen.drop_vars(dimn_interface)
    
    return uda_cen


def add_network_cellinfo(uds:xu.UgridDataset):
    """
    Reads a UgridDataset with a minimal topology as it occurs in dflowfm netfiles,
    this contains only node_x, node_y and edge_node_connectivity. xugrid reads this as Ugrid1d
    It converts it to a Ugrid2d grid which includes the face_node_connectivity.
    It also couples the variables like NetNode_z again, but drops the edge_dims to avoid conflicts.
    """
    #check if indeed 1D grid object
    assert isinstance(uds.grid, xu.Ugrid1d)
    
    # derive meshkernel from grid with Mesh1d
    mk1 = uds.grid.meshkernel
    mk_mesh1d = mk1.mesh1d_get()
    crs = uds.grid.crs
    
    # use Mesh1d nodes/edgenodes info for generation of meshkernel with Mesh2d
    is_geographic = uds.grid.is_geographic
    from dfm_tools.meshkernel_helpers import geographic_to_meshkernel_projection
    projection = geographic_to_meshkernel_projection(is_geographic)
    mk_mesh2d = meshkernel.Mesh2d(mk_mesh1d.node_x, mk_mesh1d.node_y, mk_mesh1d.edge_nodes)
    mk2 = meshkernel.MeshKernel(projection=projection)
    mk2.mesh2d_set(mk_mesh2d)
    mesh2d_grid = mk2.mesh2d_get() #this is a more populated version of mk_mesh2d, needed for xugrid
    #TODO: we have to supply is_geographic twice, necessary?
    # also "projected" is opposite of "is_geographic" according to the docstring
    xu_grid = xu.Ugrid2d.from_meshkernel(mesh2d_grid, projected = not is_geographic, crs=crs)
    
    # convert uds.obj (non-grid vars from dataset) to new xugrid standards
    rename_dims_dict = {uds.grid.node_dimension:xu_grid.node_dimension,
                        # uds.grid.edge_dimension:xu_grid.edge_dimension,
                        }
    uds_obj = uds.obj.rename_dims(rename_dims_dict)
    # drop edge dim since dim size can be changed in case of hanging edges
    uds_obj = uds_obj.drop_dims(uds.grid.edge_dimension)
    # drop node coordinate var since the order might conflict with the new grid and it has no meaning
    uds_obj = uds_obj.drop_vars(uds.grid.node_dimension)
    
    # combine again with rest of dataset
    uds_withinfo = xu.UgridDataset(obj=uds_obj, grids=[xu_grid])
    
    return uds_withinfo


def enrich_rst_with_map(ds_rst:xr.Dataset):
    """
    enriches rst dataset with topology from mapfile in the same folder.
    in order to merge, the bnd dimension has to be dropped.
    Also, the domain variable is currently not merged, so ghostcells are not removed.

    Parameters
    ----------
    ds_rst : xr.Dataset
        DESCRIPTION.

    Returns
    -------
    ds_rst : TYPE
        DESCRIPTION.

    """
    
    # derive fname of mapfile and open dataset
    file_nc_rst = ds_rst.encoding["source"]
    dir_output = os.path.dirname(file_nc_rst)
    fname_rst = os.path.basename(file_nc_rst)
    fname_map = "_".join(fname_rst.split("_")[:-3]) + "_map.nc"
    file_nc_map = os.path.join(dir_output, fname_map)
    ds_map = xr.open_dataset(file_nc_map)
    
    # remove unassociated edges from mapfile to align with rst file
    ds_map = remove_unassociated_edges(ds_map)
    
    # enrich rst file with topology variables from mapfile
    topology_varn = ds_map.ugrid_roles.topology[0]
    topo_var = ds_map[topology_varn]
    ds_rst[topology_varn] = topo_var
    nodecoords = topo_var.attrs["node_coordinates"]
    for topovar in nodecoords.split():
        ds_rst[topovar] = ds_map[topovar]
        ds_rst = ds_rst.set_coords(topovar)
    edgecoords = topo_var.attrs["edge_coordinates"]
    for topovar in edgecoords.split():
        ds_rst[topovar] = ds_map[topovar]
        ds_rst = ds_rst.set_coords(topovar)
    fnc_varn = topo_var.attrs["face_node_connectivity"]
    ds_rst[fnc_varn] = ds_map[fnc_varn]
    
    # rename old dims
    if 'nFlowLinkPts' in ds_rst.dims and 'nNetLinkPts' in ds_rst.dims:
        ds_rst = ds_rst.rename({'nFlowLinkPts':'nNetLinkPts'})
    if 'nFlowElem' in ds_rst.dims:
        facedim = topo_var.attrs["face_dimension"]
        ds_rst = ds_rst.rename({'nFlowElem':facedim})
    if 'nNetElem' in ds_rst.dims:
        facedim = topo_var.attrs["face_dimension"]
        ds_rst = ds_rst.rename({'nNetElem':facedim})
    if "nNetLink" in ds_rst.dims:
        edgedim = topo_var.attrs["edge_dimension"]
        ds_rst = ds_rst.rename({"nNetLink":edgedim})
    if "nNetElemMaxNode" in ds_rst.dims:
        mfndim = topo_var.attrs["max_face_nodes_dimension"]
        ds_rst = ds_rst.rename({"nNetElemMaxNode":mfndim})
    
    # drop bnd dim since it cannot be merged by xugrid
    if "nFlowElemBnd" in ds_rst.dims:
        ds_rst = ds_rst.drop_dims("nFlowElemBnd")
    
    return ds_rst
