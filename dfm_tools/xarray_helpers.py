# -*- coding: utf-8 -*-
"""
Created on Fri Oct 14 19:58:36 2022

@author: veenstra
"""

from netCDF4 import Dataset
import xarray as xr


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
        ds = ds.set_index({dim:coord_str})
        
        #drop duplicate indices (stations/crs/gs), this avoids "InvalidIndexError: Reindexing only valid with uniquely valued Index objects"
        duplicated_keepfirst = ds[dim].to_series().duplicated(keep='first')
        if duplicated_keepfirst.sum()>0:
            print(f'dropping {duplicated_keepfirst.sum()} duplicate "{coord}" labels to avoid InvalidIndexError')
            ds = ds[{dim:~duplicated_keepfirst}]
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


def Dataset_varswithdim(ds,dimname):
    if dimname not in ds.dims:
        raise Exception(f'dimension {dimname} not in dataset, available are: {list(ds.dims)}')
    
    varlist_keep = []
    for varname in ds.variables.keys():
        if dimname in ds[varname].dims:
            varlist_keep.append(varname)
    ds = ds[varlist_keep]
    
    return ds



