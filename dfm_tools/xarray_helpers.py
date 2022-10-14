# -*- coding: utf-8 -*-
"""
Created on Fri Oct 14 19:58:36 2022

@author: veenstra
"""

from netCDF4 import Dataset
import xarray as xr


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