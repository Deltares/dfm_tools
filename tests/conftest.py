# -*- coding: utf-8 -*-
"""
Created on Thu May 16 15:51:44 2024

@author: veenstra
"""
import pytest
import dfm_tools as dfmt
import os
import xarray as xr
import pandas as pd
import numpy as np
import geopandas as gpd


@pytest.fixture
def data_dcsm_gdf():
    # dummy gdf
    points_x = [-9.25, -9.5, -9.75]
    points_y = [43, 43, 43]
    points_n = [f'DCSM-FM_OB_all_20181108_{i+1:04d}' for i in range(3)]
    geom = gpd.points_from_xy(x=points_x, y=points_y)
    gdf_points = gpd.GeoDataFrame(geometry=geom, crs='EPSG:4326')
    gdf_points['station_id'] = points_n
    data_dcsm_gdf = gdf_points
    return data_dcsm_gdf


@pytest.fixture
def cmems_dataset_notime():
    # use hardcoded depth varname/dimname to simulate CMEMS dataset
    ds = xr.Dataset()
    so_np = np.array([[[35.819576, 35.82568 , 35.82873 ],
                       [35.819576, 35.824154, 35.831783],
                       [35.822628, 35.824154, 35.82873 ]],
                      
                      [[35.802788, 35.80584 , 35.815   ],
                       [35.815   , 35.810417, 35.821102],
                       [35.824154, 35.813473, 35.81805 ]],
                      
                      [[35.786003, 35.789055, np.nan],
                       [35.807365, 35.796684, np.nan],
                       [35.824154, 35.80584 , np.nan]],
                      
                      [[35.776848, np.nan,    np.nan],
                       [35.792107, np.nan,    np.nan],
                       [35.822628, np.nan,    np.nan]],
                                              
                      [[np.nan, np.nan,    np.nan],
                       [np.nan, np.nan,    np.nan],
                       [np.nan, np.nan,    np.nan]]])
    ds['so'] = xr.DataArray(so_np,dims=('depth','latitude','longitude'))
    ds['so'] = ds['so'].assign_attrs({"units":"dummyunit"})
    lons = [-9.6, -9.5, -9.4]
    lats = [42.9, 43.0, 43.1]
    depths = [-0.494025, -1.541375, -2.645669, -3.819495, -5.078224]
    
    depth_attrs = {'positive': 'up'}
    
    ds['longitude'] = xr.DataArray(lons, dims=('longitude'))
    ds['latitude'] = xr.DataArray(lats, dims=('latitude'))
    ds['depth'] = xr.DataArray(depths, dims=('depth')).assign_attrs(depth_attrs)
    cmems_dataset_notime = ds
    return cmems_dataset_notime


@pytest.fixture
def cmems_dataset_4times(cmems_dataset_notime):
    ds_notime = cmems_dataset_notime.copy()
    ds = xr.concat(4*[ds_notime.expand_dims('time')],dim='time')
    ds['time'] = xr.DataArray([-12,12,36,60],dims='time').assign_attrs({'standard_name':'time','units':'hours since 2020-01-01'})
    ds = xr.decode_cf(ds)
    cmems_dataset_4times = ds
    return cmems_dataset_4times


@pytest.fixture
def file_nc_era5_pattern(tmp_path):
    """
    requires CDS credentials
    """
    date_min = '2010-01-31'
    date_max = '2010-02-01'
    longitude_min, longitude_max, latitude_min, latitude_max =    2,   3,  51, 52 #test domain
    variables_era5 = ['msl']#'v10n'] # check variables_dict in dfmt.download_ERA5() for valid names
    for varkey in variables_era5:
        dfmt.download_ERA5(varkey, 
                           longitude_min=longitude_min, longitude_max=longitude_max, latitude_min=latitude_min, latitude_max=latitude_max,
                           date_min=date_min, date_max=date_max,
                           dir_output=tmp_path, overwrite=True)
    
    # assert downloaded files
    file_nc_era5_pattern = os.path.join(tmp_path, "*.nc")
    return file_nc_era5_pattern


@pytest.fixture
def ds_era5_empty():
    # create dummy dataset
    ds_era5_empty = xr.Dataset()
    ds_era5_empty['longitude'] = xr.DataArray()
    ds_era5_empty['latitude'] = xr.DataArray()
    ds_era5_empty = ds_era5_empty.set_coords(["longitude","latitude"])
    time_data = pd.date_range('2010-01-31', '2010-02-01', freq="3h")
    ds_era5_empty['time'] = xr.DataArray(time_data, dims='time')
    return ds_era5_empty
