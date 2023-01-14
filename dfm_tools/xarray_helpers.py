# -*- coding: utf-8 -*-
"""
Created on Fri Oct 14 19:58:36 2022

@author: veenstra
"""

from netCDF4 import Dataset
import xarray as xr
import numpy as np
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


def Dataset_varswithdim(ds,dimname):
    if dimname not in ds.dims:
        raise Exception(f'dimension {dimname} not in dataset, available are: {list(ds.dims)}')
    
    varlist_keep = []
    for varname in ds.variables.keys():
        if dimname in ds[varname].dims:
            varlist_keep.append(varname)
    ds = ds[varlist_keep]
    
    return ds


def get_vertical_dimensions(uds):
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
    Opmerkingen HB:
        - Voor data op de edges zou het ook werken, maar dan is nog een andere isel noodzakelijk, specifiek voor de edge data.
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
    
    file_nc = 'p:\\1204257-dcsmzuno\\2006-2012\\3D-DCSM-FM\\A18b_ntsu1\\DFM_OUTPUT_DCSM-FM_0_5nm\\DCSM-FM_0_5nm_0000_map.nc' #3D DCSM
    file_nc = 'p:\\11206813-006-kpp2021_rmm-2d\\C_Work\\31_RMM_FMmodel\\computations\\model_setup\\run_207\\results\\RMM_dflowfm_0000_map.nc' #RMM 2D
    file_nc = 'p:\\1230882-emodnet_hrsm\\GTSMv5.0\\runs\\reference_GTSMv4.1_wiCA_2.20.06_mapformat4\\output\\gtsm_model_0*_map.nc' #GTSM 2D
    file_nc = 'p:\\11208053-005-kpp2022-rmm3d\\C_Work\\01_saltiMarlein\\RMM_2019_computations_02\\computations\\theo_03\\DFM_OUTPUT_RMM_dflowfm_2019\\RMM_dflowfm_2019_0*_map.nc' #RMM 3D
    file_nc = 'p:\\archivedprojects\\11203379-005-mwra-updated-bem\\03_model\\02_final\\A72_ntsu0_kzlb2\\DFM_OUTPUT_MB_02\\MB_02_0000_map.nc'
    Timings (open_dataset/merge_partitions): (update timings after new xugrid code)
        - DCSM 3D 20 partitions  367 timesteps: 219.0/4.68 sec
        - RMM  2D  8 partitions  421 timesteps:  60.6/1.4+0.1 sec
        - GTSM 2D  8 partitions  746 timesteps:  73.8/6.4+0.1 sec
        - RMM  3D 40 partitions  146 timesteps: 166.0/3.6+0.5 sec
        - MWRA 3D 20 partitions 2551 timesteps: 826.2/3.4+1.2 sec
    """
    
    dtstart_all = dt.datetime.now()
    if isinstance(file_nc,list):
        file_nc_list = file_nc
    else:
        file_nc_list = glob.glob(file_nc)
    if len(file_nc_list)==0:
        raise Exception('file(s) not found, empty file_nc_list')
    
    print(f'>> xu.open_dataset() for {len(file_nc_list)} partition(s): ',end='')
    dtstart = dt.datetime.now()
    partitions = []
    for iF, file_nc_one in enumerate(file_nc_list):
        print(iF,end=' ')
        partitions.append(xu.open_dataset(file_nc_one, chunks=chunks)) #TODO: speed up, for instance by doing decode after merging?
    print(': ',end='')
    print(f'{(dt.datetime.now()-dtstart).total_seconds():.2f} sec')
    
    print(f'>> xu.merge_partitions() with {len(file_nc_list)} partition(s): ',end='')
    dtstart = dt.datetime.now()
    ds_merged_xu = xu.merge_partitions(partitions)
    print(f'{(dt.datetime.now()-dtstart).total_seconds():.2f} sec')
    
    #add non face/node/edge variables back to merged dataset
    #TODO: add this to xugrid? contains data_vars ['projected_coordinate_system', 'timestep', 'mesh2d_layer_z', 'mesh2d_interface_z']
    #TODO: xugrid has .crs property, projected_coordinate_system/wgs should be updated to be crs so it will be automatically handled? >> make dflowfm issue
    facedim = ds_merged_xu.grid.face_dimension
    nodedim = ds_merged_xu.grid.node_dimension
    edgedim = ds_merged_xu.grid.edge_dimension
    ds_rest = partitions[0].drop_dims([facedim,nodedim,edgedim])
    for data_var in ds_rest.data_vars:
        ds_merged_xu[data_var] = ds_rest[data_var]

    gridname = ds_merged_xu.grid.name
    
    #update important variable names #TODO: add laydim/interfacedim to xugrid dataset? (like face_dimension property)
    rename_dict = {}
    # layer_nlayers_opts = [f'{gridname}_nLayers','laydim'] # options for old layer dimension name #TODO: others from get_varname_fromnc: ['nmesh2d_layer_dlwq','LAYER','KMAXOUT_RESTR','depth'
    # for opt in layer_nlayers_opts:
    #     if opt in ds_merged_xu.dims:
    #         #print(f'hardcoded replacing {opt} with nmesh2d_layer. Auto would replace "{partitions[0].ugrid.grid.to_dataset().mesh2d.vertical_dimensions}"')
    #         rename_dict.update({opt:f'n{gridname}_layer'})
    # layer_ninterfaces_opts = [f'{gridname}_nInterfaces']
    # for opt in layer_ninterfaces_opts:
    #     if opt in ds_merged_xu.dims:
    #         rename_dict.update({opt:f'n{gridname}_interface'})
    
    #TODO: below works if xugrid handles arbitrary grid names
    # gridspecs = partitions[0].ugrid.grid.to_dataset()[gridname]
    # if hasattr(gridspecs,'vertical_dimensions'):
    #     layer_dimn = gridspecs.layer_dimension
    #     rename_dict.update({layer_dimn:f'n{gridname}_layer'}) #old varnames for mesh2d_layer: 'mesh2d_nLayers','laydim','nmesh2d_layer_dlwq'
    # else:
    #     print(f'no layer dimension found in gridspecs ds.{gridname}')
    
    ds_merged_xu = ds_merged_xu.rename(rename_dict)
    
    varlist_onepart = list(partitions[0].variables.keys())
    varlist_merged = list(ds_merged_xu.variables.keys())
    varlist_dropped_bool = ~pd.Series(varlist_onepart).isin(varlist_merged)
    varlist_dropped = pd.Series(varlist_onepart).loc[varlist_dropped_bool]
    if varlist_dropped_bool.any():
        print(f'WARNING: some variables dropped with merging of partitions:\n{varlist_dropped}')
    #TODO: with DCSM the var mesh2d_flowelem_zw is dropped for unknown reason, fix this in xugrid
    
    print(f'>> dfmt.open_partitioned_dataset() total: {(dt.datetime.now()-dtstart_all).total_seconds():.2f} sec')
    return ds_merged_xu

