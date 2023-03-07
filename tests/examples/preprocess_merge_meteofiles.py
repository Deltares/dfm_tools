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

mode = 'HIRLAM_meteo' #'ERA5_wind_pressure' # 'HARMONIE' 'HIRLAM_meteo' 'HIRLAM_meteo-heatflux' 'HYCOM' 'ERA5_wind_pressure' 'ERA5_heat_model' 'ERA5_radiation' 'ERA5_rainfall'

script_tstart = dt.datetime.now()

if 'HIRLAM' in mode:
    if mode == 'HIRLAM_meteo': #1year voor meteo crasht (HIRLAM72_*\\h72_*) door conflicting dimension sizes, sourcefolders opruimen? meteo_heatflux folders zijn schoner dus daar werkt het wel
        dir_data = 'P:\\1204257-dcsmzuno\\2014\\data\\meteo\\HIRLAM72_*' #files contain: ['air_pressure_fixed_height','northward_wind','eastward_wind']
    elif mode == 'HIRLAM_meteo-heatflux':
        dir_data = 'P:\\1204257-dcsmzuno\\2014\\data\\meteo-heatflux\\HIRLAM72_*' # files contain: ['dew_point_temperature','air_temperature','cloud_area_fraction']
    fn_match_pattern = 'h72_20131*.nc'
    file_out_prefix = 'h72_'
    preprocess = dfmt.preprocess_hirlam #temporary fix for dims/vars with same names
    time_slice = slice('2013-12-30','2014-01-01')
elif mode == 'HARMONIE':
    dir_data = 'p:\\1204257-dcsmzuno\\data\\meteo\\HARMONIE\\nc\\air_*' #many invalid files, so subsetting here
    fn_match_pattern = 'HARMONIE_*_2020_*.nc'
    file_out_prefix = 'HARMONIE_'
    preprocess = None
    time_slice = slice('2020-01-01','2020-02-01')
elif mode == 'HYCOM':
    dir_data = 'c:\\DATA\\dfm_tools_testdata\\GLBu0.08_expt_91.2'
    fn_match_pattern = 'HYCOM_ST_GoO_*.nc'
    file_out_prefix = 'HYCOM_ST_GoO_'
    preprocess = None
    time_slice = slice('2016-04-28','2016-05-03')
    #rename_variables = {'salinity':'so', 'water_temp':'thetao'}
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
    file_out_prefix = f'era5_{"_".join(varkey_list)}'
    dir_data = 'p:\\metocean-data\\open\\ERA5\\data\\Irish_North_Baltic_Sea\\*'
    time_slice = slice('2013-12-30','2014-01-01')
    #time_slice = slice('2005-01-01','2022-01-01') #for performance checking, was 12 minutes (for which case?)
    preprocess = dfmt.preprocess_ERA5 #reduce expver dimension if present
elif mode == 'WOA': #fails since 0000 is not a valid year.
    dir_data = r'p:\1204257-dcsmzuno\data\WOA13'
    fn_match_pattern = 'woa13_decav_s*.nc'
    file_out_prefix = 'woa13_decav_s_'
    preprocess = dfmt.preprocess_woa #add 360-day calendar unit to time attrs before decode_cf
    time_slice = slice('0000-01-16','0000-03-16')
else:
    raise Exception('ERROR: wrong mode %s'%(mode))

dir_output = '.'
if not os.path.exists(dir_output):
    os.makedirs(dir_output)

file_nc = os.path.join(dir_data,fn_match_pattern)

data_xr_tsel = dfmt.merge_meteofiles(file_nc=file_nc, time_slice=time_slice, 
                                     preprocess=preprocess,
                                     add_global_overlap=False, #GTSM specific: extend data beyond -180 to 180 longitude
                                     zerostart=False) #GTSM specific: extend data with 0-value fields 1 and 2 days before all_tstart

#write to netcdf file
print('writing file (can take a while)')
times_pd_str = data_xr_tsel['time'].to_series().dt.strftime("%Y%m%d")
file_out = os.path.join(dir_output, f'{file_out_prefix}_{times_pd_str[0]}to{times_pd_str[-1]}_{mode}.nc')
data_xr_tsel.to_netcdf(file_out)


print('loading outputfile')
with xr.open_dataset(file_out) as data_xr_check:
    for varkey in data_xr_check.data_vars:
        print(f'plotting {varkey}')
        #continue #uncomment to skip plotting
        fig,ax1 = plt.subplots()
        if 'HIRLAM' in mode:
            data_xr_check[varkey].isel(time=0).plot(ax=ax1,x='longitude',y='latitude') #x/y are necessary since coords are not 1D and dims
        elif 'depth' in data_xr_tsel[varkey].coords:
            data_xr_check[varkey].isel(time=0).sel(depth=0).plot(ax=ax1)
        else:
            data_xr_check[varkey].isel(time=0).plot(ax=ax1)
        fig.savefig(file_out.replace('.nc',f'_{varkey}'))

script_telapsed = (dt.datetime.now()-script_tstart)
print(f'elapsed time: {script_telapsed}')

