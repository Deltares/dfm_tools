# -*- coding: utf-8 -*-
"""
Created on Fri Oct 14 19:58:36 2022

@author: veenstra
"""

from netCDF4 import Dataset
import xarray as xr


def preprocess_hisnc(ds, drop_duplicate_stations=True, silent=True):
    """
    Loop over keys of dict (potential dimensions) and set as Dataset index if present, to enable label based indexing. If duplicate labels are found (like duplicate stations), these are dropped to avoid indexing issues.
    
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
    dim_name_dict = {'stations':'station_name',
                     'cross_section':'cross_section_name',
                     'general_structures':'general_structure_id',
                     'source_sink':'source_sink_name',
                     }
    
    for ds_dim in ds.dims.keys():
        if ds_dim in ['time']:#,'laydim','laydimw','source_sink_pts','nFlowElemWithBnd','nFlowElemContourPts','nNetLink','nNetLinkPts','nFlowLink','nFlowLinkPts','station_geom_nNodes','source_sink_geom_nNodes']:
            continue
        if not ds_dim in dim_name_dict.keys():
            if silent:
                continue
            print(f'WARNING: dimension "{ds_dim}" does not have an index yet, you could extend the dim_name_dict')
    
    for dim in dim_name_dict.keys():
        if not dim in ds.dims:
            continue
        name = dim_name_dict[dim]
        name_str = f'{name}_str' #avoid losing the original variable by creating a new name
        ds[name_str] = ds[name].load().str.decode('utf-8',errors='ignore').str.strip() #.load() is essential to convert not only first letter of string.
        ds = ds.set_index({dim:name_str})
        
        #drop duplicate indices (stations/crs/gs), this avoids issues later on
        duplicated_keepfirst = ds[dim].to_series().duplicated(keep='first')
        if duplicated_keepfirst.sum()>0 and drop_duplicate_stations:
            print(f'dropping {duplicated_keepfirst.sum()} duplicates in "{name}" to avoid indexing issues. Prevent this with preprocess_hisnc(drop_duplicate_stations=False).')
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