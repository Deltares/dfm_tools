# -*- coding: utf-8 -*-
"""
Created on Wed Oct 27 22:07:34 2021

@author: veenstra

script to merge netcdf file from several folders into one file (concatenate time dimension)
too slow? Change fn_match_pattern to match less files if possible, or provide a list of (less) files as file_nc

"""

from netCDF4 import Dataset
import glob
import re
import datetime as dt
import os
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
plt.close('all')

mode = 'ERA5_heat_model'# 'HIRLAM_meteo' 'HIRLAM_meteo-heatflux' 'HYCOM' 'ERA5_wind_pressure' ERA5_heat_model ERA5_radiation ERA5_rainfall
all_tstart = dt.datetime(2013,12,30) # HIRLAM and ERA5
all_tstop = dt.datetime(2014,1,1)
#all_tstart = dt.datetime(2016,4,28) # HYCOM
#all_tstop = dt.datetime(2016,5,3)

script_tstart = dt.datetime.now()

if 'HIRLAM' in mode:
    if mode == 'HIRLAM_meteo': #1year voor meteo crasht (HIRLAM72_*\\h72_*) door conflicting dimension sizes, sourcefolders opruimen? meteo_heatflux folders zijn schoner dus daar werkt het wel
        dir_data = 'P:\\1204257-dcsmzuno\\2014\\data\\meteo\\HIRLAM72_*' #files contain: ['air_pressure_fixed_height','northward_wind','eastward_wind']
    elif mode == 'HIRLAM_meteo-heatflux':
        dir_data = 'P:\\1204257-dcsmzuno\\2014\\data\\meteo-heatflux\\HIRLAM72_*' # files contain: ['dew_point_temperature','air_temperature','cloud_area_fraction']
    fn_match_pattern = 'h72_20131*.nc'
    file_out_prefix = 'h72_'
    drop_variables = ['x','y'] #will be added again as longitude/latitude, this is a workaround
    rename_variables = None
elif mode == 'HYCOM':
    dir_data = 'c:\\DATA\\dfm_tools_testdata\\GLBu0.08_expt_91.2'
    fn_match_pattern = 'HYCOM_ST_GoO_*.nc'
    file_out_prefix = fn_match_pattern.replace('*.nc','')
    drop_variables = None
    rename_variables = {'salinity':'so', 'water_temp':'thetao'}
elif 'ERA5' in mode: #TODO: generates "PerformanceWarning: Slicing is producing a large chunk.", probably because of lots of files but probably also solveable
    # TODO: add features from c:\DATA\hydro_tools\ERA5\ERA52DFM.py (except for varRhoair_alt)
    if mode=='ERA5_wind_pressure':
        varkey_list = ['chnk','mslp','u10n','v10n'] #charnock, mean_sea_level_pressure, 10m_u_component_of_neutral_wind, 10m_v_component_of_neutral_wind
    elif mode=='ERA5_heat_model':
        varkey_list = ['d2m','t2m','tcc'] # 2m_dewpoint_temperature, 2m_temperature, total_cloud_cover
    elif mode=='ERA5_radiation':
        varkey_list = ['ssr','strd'] # surface_net_solar_radiation, surface_thermal_radiation_downwards
    elif mode=='ERA5_rainfall':
        varkey_list = ['mer','mtpr'] # mean_evaporation_rate, mean_total_precipitation_rate
    dir_data = 'p:\\metocean-data\\open\\ERA5\\data\\Irish_North_Baltic_Sea\\*' #TODO: add other vars with * (support separate files)
    fn_match_pattern = f'era5_.*({"|".join(varkey_list)})_.*\.nc' #to filter with re
    file_out_prefix = f'era5_{"_".join(varkey_list)}'
    drop_variables = None
    rename_variables = None
else:
    raise Exception('ERROR: wrong mode %s'%(mode))

dir_output = '.' #os.path.join(dir_data,'merged_netcdf_files_JV')
if not os.path.exists(dir_output):
    os.makedirs(dir_output)


#generating file list #TODO: can be more efficient
if '.*' in fn_match_pattern: #regex complex stuff
    file_pd_allnc = pd.Series(glob.glob(os.path.join(dir_data,'*.nc'))) #all nc files
    file_pd_allnc_varnames = file_pd_allnc.str.extract(fn_match_pattern)[0].str.replace('_','')
    file_pd_allnc_df = pd.DataFrame({'file':file_pd_allnc,'varname':file_pd_allnc_varnames})
    file_pd_pattern = file_pd_allnc_df.loc[~file_pd_allnc_df['varname'].isnull()]
    file_nc = file_pd_pattern['file'].tolist()
    file_list = file_nc
else:
    file_nc = os.path.join(dir_data,fn_match_pattern)
    file_list = glob.glob(file_nc)

print(f'opening multifile dataset of {len(file_list)} files matching "{fn_match_pattern}" (can take a while with lots of files)',end='')
#open_mfdataset when dropping x/y since both varname as dimname, which is not possible in xarray
data_xr = xr.open_mfdataset(file_nc,
                            drop_variables=drop_variables, #necessary since dims/vars with equal names are not allowed by xarray, add again later and requested matroos to adjust netcdf format.
                            parallel=True, #speeds up the process
                            #concat_dim="time", combine="nested", data_vars='minimal', coords='minimal', compat='override', #optional vars to look into: https://docs.xarray.dev/en/stable/user-guide/io.html#reading-multi-file-datasets
                            )
print('...done')

if 'HIRLAM' in mode:
    print('adding x/y variables again as lon/lat')
    #add xy as variables again with help of NetCDF4 #TODO: this part is hopefully temporary, necessary since variables cannot have the same name as dimensions in xarray
    file_nc_one = glob.glob(os.path.join(dir_data,'h72_201312.nc'))[0]
    data_nc = Dataset(file_nc_one)
    data_nc_x = data_nc['x']
    data_nc_y = data_nc['y']
    data_xr['longitude'] = xr.DataArray(data_nc_x,dims=data_nc_x.dimensions,attrs=data_nc_x.__dict__)
    data_xr['latitude'] = xr.DataArray(data_nc_y,dims=data_nc_y.dimensions,attrs=data_nc_y.__dict__)
    data_xr = data_xr.set_coords(['longitude','latitude'])

#rename variables
data_xr = data_xr.rename(rename_variables)
varkeys = data_xr.variables.mapping.keys()
#data_xr.attrs['comment'] = 'merged with dfm_tools from https://github.com/openearth/dfm_tools' #TODO: add something like this


#select time and do checks
print('time slicing and drop duplicates')
data_xr_tsel = data_xr.sel(time=slice(all_tstart,all_tstop))
data_xr_tsel = data_xr_tsel.sel(time=~data_xr_tsel.get_index('time').duplicated()) #drop duplicate timesteps
times_pd = data_xr_tsel['time'].to_series()

#check if there are times selected
if len(times_pd)==0:
    raise Exception('ERROR: no times selected, check tstart/tstop and file_nc')

#check if there are no gaps (more than one timestep)
timesteps_uniq = times_pd.diff().iloc[1:].unique()
if len(timesteps_uniq)>1:
    raise Exception(f'ERROR: gaps found in selected dataset (are there sourcefiles missing?), unique timesteps (hour): {timesteps_uniq/1e9/3600}')

#check if requested times are available in selected files (in times_pd)
if not all_tstart in times_pd.index:
    raise Exception(f'ERROR: all_tstart="{all_tstart}" not in selected files, timerange: "{times_pd.index[0]}" to "{times_pd.index[-1]}"')
if not all_tstop in times_pd.index:
    raise Exception(f'ERROR: all_tstop="{all_tstop}" not in selected files, timerange: "{times_pd.index[0]}" to "{times_pd.index[-1]}"')

#get actual tstart/tstop strings from times_pd
tstart_str = times_pd.index[0].strftime("%Y%m%d")
tstop_str = times_pd.index[-1].strftime("%Y%m%d")


#do conversions: 'dew_point_temperature' (K to C) ,'air_temperature' (K to C),'cloud_area_fraction' (0-1 to %)
for var_temperature_KtoC in ['air_temperature','dew_point_temperature','d2m','t2m']:
    if var_temperature_KtoC in varkeys:
        print(f'converting temperature units (K to C) for {var_temperature_KtoC}')
        data_xr_tsel[var_temperature_KtoC].attrs['units'] = 'C'
        data_xr_tsel[var_temperature_KtoC] = data_xr_tsel[var_temperature_KtoC] - 273.15
for var_fractoperc in ['cloud_area_fraction']: #TODO: add 'tcc'
    if var_fractoperc in varkeys:
        print(f'converting fraction to percentage for {var_fractoperc}')
        #data_xr_tsel[var_fractoperc].attrs['units'] = '%' #unit is al %
        data_xr_tsel[var_fractoperc] = data_xr_tsel[var_fractoperc] * 100
    


#write to netcdf file
print('writing file (can take a while)')
file_out = os.path.join(dir_output, f'{file_out_prefix}_{tstart_str}to{tstop_str}_{mode}.nc')
data_xr_tsel.to_netcdf(file_out) #TODO: maybe add different reftime?

print('loading outputfile')
with xr.open_dataset(file_out) as data_xr_check:
    #data_xr_check.close()
    if 1: #optionally plot to check
        print('plotting')
        for varkey in varkeys:
            if 'time' not in data_xr_tsel[varkey].coords:
                continue
            if varkey=='time':
                continue
            fig,ax1 = plt.subplots()
            if 'HIRLAM' in mode:
                data_xr_check[varkey].isel(time=0).plot(ax=ax1,x='longitude',y='latitude') #x/y are necessary since coords are not 1D and dims
            elif 'depth' in data_xr_tsel[varkey].coords:
                data_xr_check[varkey].isel(time=0).sel(depth=0).plot(ax=ax1)
            else:
                data_xr_check[varkey].isel(time=0).plot(ax=ax1)

script_telapsed = (dt.datetime.now()-script_tstart)
print(f'elapsed time: {script_telapsed}')





