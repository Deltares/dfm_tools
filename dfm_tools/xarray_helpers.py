# -*- coding: utf-8 -*-
"""
Created on Fri Oct 14 19:58:36 2022

@author: veenstra
"""

from netCDF4 import Dataset
import xarray as xr
import xugrid as xu
import matplotlib.pyplot as plt
plt.close('all')
import datetime as dt
import glob
import pandas as pd
import warnings


def preprocess_hisnc(ds):
    """
    Look for dim/coord combination and use this for Dataset.set_index(), to enable station/gs/crs/laterals label based indexing. If duplicate labels are found (like duplicate stations), these are dropped to avoid indexing issues.
    
    Parameters
    ----------
    ds : xarray.Dataset
        DESCRIPTION.
    drop_duplicate_stations : TYPE, optional
        DESCRIPTION. The default is True.

    Returns
    -------
    ds : TYPE
        DESCRIPTION.

    """
    
    #generate dim_coord_dict to set indexes, this will be something like {'stations':'station_name','cross_section':'cross_section_name'} after loop
    dim_coord_dict = {}
    for ds_coord in ds.coords.keys():
        ds_coord_dtype = ds[ds_coord].dtype
        ds_coord_dim = ds[ds_coord].dims[0] #these vars always have only one dim
        if ds_coord_dtype.str.startswith('|S'): #these are station/crs/laterals/gs names/ids
            dim_coord_dict[ds_coord_dim] = ds_coord
    
    #loop over dimensions and set corresponding coordinates/variables from dim_coord_dict as their index
    for dim in dim_coord_dict.keys():
        coord = dim_coord_dict[dim]
        coord_str = f'{coord}'#_str' #avoid losing the original variable by creating a new name
        ds[coord_str] = ds[coord].load().str.decode('utf-8',errors='ignore').str.strip() #.load() is essential to convert not only first letter of string.
        #ds = ds.set_index({dim:[coord_str,'station_x_coordinate','station_y_coordinate']}) #nearest station: "ValueError: multi-index does not support ``method`` and ``tolerance``"  #slice x/y: "TypeError: float() argument must be a string or a number, not 'slice'"
        ds = ds.set_index({dim:coord_str})
        
        #drop duplicate indices (stations/crs/gs), this avoids "InvalidIndexError: Reindexing only valid with uniquely valued Index objects"
        duplicated_keepfirst = ds[dim].to_series().duplicated(keep='first')
        if duplicated_keepfirst.sum()>0:
            print(f'dropping {duplicated_keepfirst.sum()} duplicate "{coord}" labels to avoid InvalidIndexError')
            ds = ds[{dim:~duplicated_keepfirst}]

    
    if 'source' in ds.attrs.keys():
        source_attr = ds.attrs["source"]
    else:
        source_attr = None
    try:
        source_attr_version = source_attr.split(', ')[1]
        source_attr_date = source_attr.split(', ')[2]
        if pd.Timestamp(source_attr_date) < dt.datetime(2020,11,28):
            warnings.warn(UserWarning(f'Your model was run with a D-FlowFM version from before 28-10-2020 ({source_attr_version} from {source_attr_date}), the layers in the hisfile are incorrect. Check UNST-2920 and UNST-3024 for more information, it was fixed from OSS 67858.'))
    except:
        #print('No source attribute present in hisfile, cannot check version')
        pass

    return ds


def preprocess_hirlam(ds):
    """
    add xy variables as longitude/latitude to avoid duplicate var/dim names
    add xy as variables again with help of NetCDF4 
    #TODO: this part is hopefully temporary, necessary since variables cannot have the same name as dimensions in xarray
    # background and future solution: https://github.com/pydata/xarray/issues/6293
    """
    
    print('adding x/y variables again as lon/lat')
    file_nc_one = ds.encoding['source']
    with Dataset(file_nc_one) as data_nc:
        data_nc_x = data_nc['x']
        data_nc_y = data_nc['y']
        ds['longitude'] = xr.DataArray(data_nc_x,dims=data_nc_x.dimensions,attrs=data_nc_x.__dict__)
        ds['latitude'] = xr.DataArray(data_nc_y,dims=data_nc_y.dimensions,attrs=data_nc_y.__dict__)
    ds = ds.set_coords(['latitude','longitude'])
    for varkey in ds.data_vars:
        del ds[varkey].encoding['coordinates'] #remove {'coordinates':'y x'} from encoding (otherwise set twice)
    return ds


def preprocess_ERA5(ds):
    """
    Reduces the expver dimension in some of the ERA5 data (mtpr and other variables), which occurs in files with very recent data. The dimension contains the unvalidated data from the latest month in the second index in the expver dimension. The reduction is done with mean, but this is arbitrary, since there is only one valid value per timestep and the other one is nan.
    """
    if 'expver' in ds.dims:
        ds = ds.mean(dim='expver')
    return ds


def Dataset_varswithdim(ds,dimname): #TODO: dit zit ook in xugrid, wordt nu gebruikt in hisfile voorbeeldscript en kan handig zijn, maar misschien die uit xugrid gebruiken?
    if dimname not in ds.dims:
        raise Exception(f'dimension {dimname} not in dataset, available are: {list(ds.dims)}')
    
    varlist_keep = []
    for varname in ds.variables.keys():
        if dimname in ds[varname].dims:
            varlist_keep.append(varname)
    ds = ds[varlist_keep]
    
    return ds


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
    gridname = uds.grid.name
    grid_info = uds.grid.to_dataset()[gridname]
    if hasattr(grid_info,'layer_dimension'):
        print('>> found layer/interface dimensions in file: ',end='')
        print(grid_info.layer_dimension, grid_info.interface_dimension) #combined in attr vertical_dimensions
        return grid_info.layer_dimension, grid_info.interface_dimension
    else:
        return None, None


def open_partitioned_dataset(file_nc, chunks={'time':1}): 
    """
    using xugrid to read and merge partitions, with some additional features (remaning old layerdim, timings, set zcc/zw as data_vars)
    Opmerkingen HB:
        - Dit werkt nu ook alleen als je enkel grid in je dataset hebt. Bij meerdere grids zouden we een keyword moeten toevoegen dat je aangeeft welke je gemerged wilt zien.    

    Parameters
    ----------
    file_nc : TYPE
        DESCRIPTION.
    chunks : TYPE, optional
        chunks={'time':1} increases performance significantly upon reading, but causes memory overloads when performing sum/mean/etc actions over time dimension (in that case 100/200 is better). The default is {'time':1}.

    Raises
    ------
    Exception
        DESCRIPTION.

    Returns
    -------
    TYPE
        DESCRIPTION.
    
    file_nc = 'p:\\1204257-dcsmzuno\\2006-2012\\3D-DCSM-FM\\A18b_ntsu1\\DFM_OUTPUT_DCSM-FM_0_5nm\\DCSM-FM_0_5nm_0*_map.nc' #3D DCSM
    file_nc = 'p:\\11206813-006-kpp2021_rmm-2d\\C_Work\\31_RMM_FMmodel\\computations\\model_setup\\run_207\\results\\RMM_dflowfm_0*_map.nc' #RMM 2D
    file_nc = 'p:\\1230882-emodnet_hrsm\\GTSMv5.0\\runs\\reference_GTSMv4.1_wiCA_2.20.06_mapformat4\\output\\gtsm_model_0*_map.nc' #GTSM 2D
    file_nc = 'p:\\11208053-005-kpp2022-rmm3d\\C_Work\\01_saltiMarlein\\RMM_2019_computations_02\\computations\\theo_03\\DFM_OUTPUT_RMM_dflowfm_2019\\RMM_dflowfm_2019_0*_map.nc' #RMM 3D
    file_nc = 'p:\\archivedprojects\\11203379-005-mwra-updated-bem\\03_model\\02_final\\A72_ntsu0_kzlb2\\DFM_OUTPUT_MB_02\\MB_02_0*_map.nc'
    Timings (open_dataset/merge_partitions):
        - DCSM 3D 20 partitions  367 timesteps: 219.0/ 4.8 sec >> merge 13.4 sec
        - RMM  2D  8 partitions  421 timesteps:  60.6/ 5.3 sec >> merge  8.6 sec
        - GTSM 2D  8 partitions  746 timesteps:  73.8/31.0 sec >> merge 37.0 sec
        - RMM  3D 40 partitions  146 timesteps: 166.0/ 7.6 sec >> merge 15.2 sec
        - MWRA 3D 20 partitions 2551 timesteps: 826.2/ 3.9 sec >> merge 73.2 sec
    
    """
    #TODO: FM-mapfiles contain wgs84/projected_coordinate_system variables. xugrid has .crs property, projected_coordinate_system/wgs84 should be updated to be crs so it will be automatically handled? >> make dflowfm issue (and https://github.com/Deltares/xugrid/issues/42)
    
    dtstart_all = dt.datetime.now()
    if isinstance(file_nc,list):
        file_nc_list = file_nc
    else:
        file_nc_list = glob.glob(file_nc)
    if len(file_nc_list)==0:
        raise Exception('file(s) not found, empty file_nc_list')
    
    print(f'>> xu.open_dataset() with {len(file_nc_list)} partition(s): ',end='')
    dtstart = dt.datetime.now()
    partitions = []
    for iF, file_nc_one in enumerate(file_nc_list):
        print(iF+1,end=' ')
        #TODO: speed up, for instance by doing decode after merging? (or is second-read than not faster anymore?) >> https://github.com/Deltares/dfm_tools/issues/225 >> c:\DATA\dfm_tools\tests\examples_workinprogress\xarray_largemapfile_profiler.py (copied MBAY partition to d-drive to check whether network causes it)
        ds = xu.open_dataset(file_nc_one, chunks=chunks)
        partitions.append(ds)
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


#TODO: remove all code below this line if xugrid.merge_partitions is fast again

from collections import defaultdict
def group_grids_by_name(partitions):
    grouped = defaultdict(list)
    for partition in partitions:
        for grid in partition.grids:
            grouped[grid.name].append(grid)

    validate_partition_topology(grouped, len(partitions))
    return grouped

def group_vars_by_ugrid_dim(data_objects, ugrid_dims):
    validate_partition_objects(data_objects)

    # Group variables by UGRID dimension.
    ds = data_objects[0]
    grouped = defaultdict(list)
    for var, da in ds.data_vars.items():
        intersection = ugrid_dims.intersection(da.dims)
        if intersection:
            if len(intersection) > 1:
                raise ValueError(
                    f"{var} contains more than one UGRID dimension: {intersection}"
                )
            dim = intersection.pop()
            grouped[dim].append(var)

    return grouped

def validate_partition_topology(grouped, n_partition: int):
    n = n_partition
    if not all(len(v) == n for v in grouped.values()):
        raise ValueError(
            f"Expected {n} UGRID topologies for {n} partitions, received: " f"{grouped}"
        )

    for name, grids in grouped.items():
        types = set(type(grid) for grid in grids)
        if len(types) > 1:
            raise TypeError(
                f"All partition topologies with name {name} should be of the "
                f"same type, received: {types}"
            )

        griddims = set(tuple(grid.dimensions) for grid in grids)
        if len(griddims) > 1:
            raise ValueError(
                f"Dimension names on UGRID topology {name} do not match "
                f"across partitions: {griddims[0]} versus {griddims[1]}"
            )

    return None

def validate_partition_objects(data_objects):
    # Check presence of variables.
    allvars = set(tuple(sorted(ds.dims)) for ds in data_objects)
    if len(allvars) > 1:
        raise ValueError(
            "These variables are present in some partitions, but not in "
            f"others: {set(allvars[0]).symmetric_difference(allvars[1])}"
        )
    # Check dimensions
    for var in allvars.pop():
        vardims = set(ds[var].dims for ds in data_objects)
        if len(vardims) > 1:
            raise ValueError(
                f"Dimensions for {var} do not match across partitions: "
                f"{vardims[0]} versus {vardims[1]}"
            )


def merge_partitions_OLD(partitions):
    types = set(type(obj) for obj in partitions)
    msg = "Expected UgridDataArray or UgridDataset, received: {}"
    if len(types) > 1:
        type_names = [t.__name__ for t in types]
        raise TypeError(msg.format(type_names))
    obj_type = types.pop()
    if obj_type not in (xu.UgridDataArray, xu.UgridDataset):
        raise TypeError(msg.format(obj_type.__name__))

    # Convert to dataset for convenience
    data_objects = [partition.obj for partition in partitions]
    data_objects = [
        obj.to_dataset() if isinstance(obj, xr.DataArray) else obj
        for obj in data_objects
    ]
    # Collect grids
    grids = [grid for p in partitions for grid in p.grids]
    ugrid_dims = set(dim for grid in grids for dim in grid.dimensions)
    grids_by_name = group_grids_by_name(partitions)
    vars_by_dim = group_vars_by_ugrid_dim(data_objects, ugrid_dims)

    merged_grids = []
    objects = data_objects
    for grids in grids_by_name.values():
        grid = grids[0]
        merged_grid, indexes = grid.merge_partitions(grids)
        merged_grids.append(merged_grid)
        for dim, dim_indexes in indexes.items():
            objects = [
                obj.isel({dim: index}) for obj, index in zip(objects, dim_indexes)
            ]

    merged = xr.Dataset()
    for dim, vars in vars_by_dim.items():
        for var in vars:
            try:
                das = [obj[var] for obj in objects]
                merged[var] = xr.concat(das, dim=dim)
            # This should mean that another dimension doesn't align, e.g. a
            # variable that depends on the n_max_node_per_face dimension,
            # e.g. for triangles and quadrangles.
            except ValueError:
                pass

    return xu.UgridDataset(merged, merged_grids)

