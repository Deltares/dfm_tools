# -*- coding: utf-8 -*-
"""
Created on Wed Oct 27 22:07:34 2021
@author: veenstra
script to merge netcdf file from several folders into one file (concatenate time dimension)
too slow? Change fn_match_pattern to match less files if possible, or provide a list of (less) files as file_nc
"""

import datetime as dt
import os
import xarray as xr
import matplotlib.pyplot as plt
plt.close('all')
import dfm_tools as dfmt

mode = 'ERA5_heat_model' # 'HIRLAM_meteo' 'HIRLAM_meteo-heatflux' 'HARMONIE' 'ERA5_wind_pressure' 'ERA5_heat_model' 'ERA5_radiation' 'ERA5_rainfall'

if 'HIRLAM' in mode:
    if mode == 'HIRLAM_meteo': #1year voor meteo crasht (HIRLAM72_*\\h72_*) door conflicting dimension sizes, sourcefolders opruimen? meteo_heatflux folders zijn schoner dus daar werkt het wel
        dir_data = 'P:\\1204257-dcsmzuno\\2014\\data\\meteo\\HIRLAM72_*' #files contain: ['air_pressure_fixed_height','northward_wind','eastward_wind']
    elif mode == 'HIRLAM_meteo-heatflux':
        dir_data = 'P:\\1204257-dcsmzuno\\2014\\data\\meteo-heatflux\\HIRLAM72_*' # files contain: ['dew_point_temperature','air_temperature','cloud_area_fraction']
    fn_match_pattern = 'h72_20131*.nc'
    file_out_prefix = 'h72_'
    kwargs = dict(preprocess=None, combine="nested", concat_dim="time") # additional keyword arguments for xarray.open_mfdataset. 
    time_slice = slice('2013-12-30','2014-01-01')
elif mode == 'HARMONIE':
    dir_data = 'p:\\1204257-dcsmzuno\\data\\meteo\\HARMONIE\\nc\\air_*' #many invalid files, so subsetting here
    fn_match_pattern = 'HARMONIE_*_2020_*.nc'
    file_out_prefix = 'HARMONIE_'
    kwargs = dict() # additional keyword arguments for xarray.open_mfdataset. 
    time_slice = slice('2020-01-01','2020-01-02')
elif 'ERA5' in mode:
    if mode=='ERA5_u10':
        varkey_list = ['u10'] #for testing
    elif mode=='ERA5_v10n':
        varkey_list = ['v10n'] #for testing
    elif mode=='ERA5_wind_pressure':
        varkey_list = ['chnk','mslp','u10n','v10n'] #charnock, mean_sea_level_pressure, 10m_u_component_of_neutral_wind, 10m_v_component_of_neutral_wind
    elif mode=='ERA5_heat_model':
        varkey_list = ['d2m','t2m','tcc'] # 2m_dewpoint_temperature, 2m_temperature, total_cloud_cover
    elif mode=='ERA5_radiation':
        varkey_list = ['ssr','strd'] # surface_net_solar_radiation, surface_thermal_radiation_downwards
    elif mode=='ERA5_rainfall':
        varkey_list = ['mer','mtpr'] # mean_evaporation_rate, mean_total_precipitation_rate
    fn_match_pattern = f'era5_.*({"|".join(varkey_list)})_.*.nc' #simpler but selects more files: 'era5_*.nc'
    file_out_prefix = f'era5_{"_".join(varkey_list)}_'
    dir_data = 'p:\\metocean-data\\open\\ERA5\\data\\Irish_North_Baltic_Sea\\*'
    kwargs = dict(preprocess=dfmt.preprocess_ERA5) # additional keyword arguments for xarray.open_mfdataset. dfmt.preprocess_ERA5: reduce expver dimension if present and prevent dtype int
    time_slice = slice('2013-12-30','2014-01-01')
    #time_slice = slice('2005-01-01','2022-01-01') #for performance checking, was 12 minutes (for which varkey_list?)
else:
    raise KeyError('ERROR: wrong mode %s'%(mode))

dir_output = '.'
os.makedirs(dir_output, exist_ok=True)

file_nc = os.path.join(dir_data,fn_match_pattern)

data_xr_tsel = dfmt.merge_meteofiles(file_nc=file_nc, time_slice=time_slice, **kwargs)

#write to netcdf file
print('>> writing file (can take a while): ',end='')
dtstart = dt.datetime.now()
times_pd = data_xr_tsel['time'].to_series()
time_start_str = times_pd.iloc[0].strftime("%Y%m%d")
time_stop_str = times_pd.iloc[-1].strftime("%Y%m%d")
file_out = os.path.join(dir_output, f'{file_out_prefix}{time_start_str}to{time_stop_str}_{mode}.nc')
data_xr_tsel.to_netcdf(file_out)
print(f'{(dt.datetime.now()-dtstart).total_seconds():.2f} sec')


#load outputfile
data_xr_check = xr.open_dataset(file_out)

for varkey in data_xr_check.data_vars:
    varsel = data_xr_check[varkey]
    if not set(['longitude','latitude']).issubset(set(varsel.coords)): #skipping vars without lat/lon coordinate
        continue
    print(f'plotting {varkey}')
    fig,ax1 = plt.subplots()
    if 'HIRLAM' in mode:
        varsel.isel(time=0).plot(ax=ax1,x='longitude',y='latitude') #x/y are necessary since coords are not 1D and dims
    elif 'depth' in data_xr_tsel[varkey].coords:
        varsel.isel(time=0).sel(depth=0).plot(ax=ax1)
    else:
        varsel.isel(time=0).plot(ax=ax1)
    fig.savefig(file_out.replace('.nc',f'_{varkey}'))
