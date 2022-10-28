# -*- coding: utf-8 -*-
"""
Created on Fri Oct 14 19:58:36 2022

@author: veenstra
"""

from netCDF4 import Dataset
import xarray as xr


def preprocess_hisnc(ds, drop_duplicate_stations=True, silent=False): #TODO: these arguments cannot be changed by user somehow
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
    dim_coord_dict = {'stations':'station_name',
                      'cross_section':'cross_section_name',
                      'general_structures':'general_structure_id',
                      'source_sink':'source_sink_name',
                      'lateral':'lateral_id',
                      }
    
    #loop over variables with cf_role=timeseries_id as attribute, these can be potentially used as index
    ds_cfrole_timeseriesid = ds.filter_by_attrs(cf_role='timeseries_id')
    for ds_coord in ds_cfrole_timeseriesid.coords.keys():
        if ds_coord in ['station_x_coordinate','station_y_coordinate','station_lon','station_lat']: #chosen index for dimension stations is variable station_name
            continue
        if not ds_coord in dim_coord_dict.values():
            print(f'WARNING: variable  {ds_coord} {ds_cfrole_timeseriesid[ds_coord].dims}  could be used as index, you could extend the dimcoord_name_dict')
    
    #loop over dimensions and set corresponding coordinates/variables from dim_coord_dict as their index
    for dim in dim_coord_dict.keys():
        if not dim in ds.dims:
            continue
        name = dim_coord_dict[dim]
        name_str = f'{name}_str' #avoid losing the original variable by creating a new name
        ds[name_str] = ds[name].load().str.decode('utf-8',errors='ignore').str.strip() #.load() is essential to convert not only first letter of string.
        ds = ds.set_index({dim:name_str})
        
        #drop duplicate indices (stations/crs/gs), this avoids issues later on
        duplicated_keepfirst = ds[dim].to_series().duplicated(keep='first')
        if duplicated_keepfirst.sum()>0 and drop_duplicate_stations:
            print(f'dropping {duplicated_keepfirst.sum()} duplicate labels in "{name}" to avoid indexing issues.')
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



