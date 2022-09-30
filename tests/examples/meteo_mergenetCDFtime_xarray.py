# -*- coding: utf-8 -*-
"""
Created on Wed Oct 27 22:07:34 2021

@author: veenstra

this test tests merging a (eg meteo) netcdf over time

"""

from netCDF4 import Dataset
import glob
import datetime as dt
import os
import xarray as xr
#import matplotlib.pyplot as plt
#plt.close('all')
#from dfm_tools.get_nc_helpers import get_ncvardimlist#, get_ncfilelist
#from dfm_tools.io.netCDF_utils import merge_netCDF_time


mode = 'ERA5'# 'HIRLAM_meteo' 'HIRLAM_meteo_heatflux' 'HYCOM'
all_tstart = dt.datetime(2013,12,30) # HIRLAM
all_tstop = dt.datetime(2014,1,1)
#all_tstart = dt.datetime(2016,4,28) # HYCOM
#all_tstop = dt.datetime(2016,5,3)

script_tstart = dt.datetime.now()

if 'HIRLAM' in mode:
    if mode == 'HIRLAM_meteo': #1year voor meteo crasht (HIRLAM72_*\\h72_*) door conflicting dimension sizes, sourcefolders opruimen? meteo_heatflux folders zijn schoner dus daar werkt het wel
        dir_data = 'P:\\1204257-dcsmzuno\\2014\\data\\meteo\\HIRLAM72_*' #files contain: ['air_pressure_fixed_height','northward_wind','eastward_wind']
    elif mode == 'HIRLAM_meteo_heatflux':
        dir_data = 'P:\\1204257-dcsmzuno\\2014\\data\\meteo-heatflux\\HIRLAM72_*' # files contain: ['dew_point_temperature','air_temperature','cloud_area_fraction']
    fn_match_pattern = 'h72_20131*.nc'
    drop_variables = ['x','y'] #will be added again as longitude/latitude, this is a workaround
    rename_variables = None
elif mode == 'HYCOM':
    dir_data = 'c:\\DATA\\dfm_tools_testdata\\GLBu0.08_expt_91.2'
    fn_match_pattern = 'HYCOM_ST_GoO_*.nc'
    drop_variables = None
    rename_variables = {'salinity':'so', 'water_temp':'thetao'}
elif mode == 'ERA5':
    data_dir = 'p:\\metocean-data\\open\\ERA5\\data\\Irish_North_Baltic_Sea'
    
else:
    raise Exception('ERROR: wrong mode %s'%(mode))

dir_output = '.' #os.path.join(dir_data,'merged_netcdf_files_JV')
if not os.path.exists(dir_output):
    os.makedirs(dir_output)


#open_mfdataset when dropping x/y since both varname as dimname, which is not possible in xarray
file_nc = os.path.join(dir_data,fn_match_pattern)
print(f'opening multifile dataset of {len(glob.glob(file_nc))} files (can take a while with lots of files)',end='')
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


print('optionally converting units')
#do conversions: 'dew_point_temperature' (K to C) ,'air_temperature' (K to C),'cloud_area_fraction' (0-1 to %)
if 'air_temperature' in varkeys:
    data_xr_tsel['air_temperature'].attrs['units'] = 'C'
    data_xr_tsel['air_temperature'] = data_xr_tsel['air_temperature'] - 273.15
if 'dew_point_temperature' in varkeys:
    data_xr_tsel['dew_point_temperature'].attrs['units'] = 'C'
    data_xr_tsel['dew_point_temperature'] = data_xr_tsel['dew_point_temperature'] - 273.15
if 'cloud_area_fraction' in varkeys:
    #data_xr_tsel['cloud_area_fraction'].attrs['units'] = '%' #unit is al %
    data_xr_tsel['cloud_area_fraction'] = data_xr_tsel['cloud_area_fraction'] * 100


#write to netcdf file
print('writing file (can take a while)')
file_prefix = fn_match_pattern.replace('_*.nc','').replace('*.nc','')

file_out = os.path.join(dir_output, f'{file_prefix}_{tstart_str}to{tstop_str}_{mode}.nc')
data_xr_tsel.to_netcdf(file_out) #TODO: maybe add different reftime?

print('loading outputfile')
with xr.open_dataset(file_out) as data_xr_check:
    #data_xr_check.close()
    if 1: #optionally plot to check
        print('plotting')
        import matplotlib.pyplot as plt
        plt.close('all')
        fig,ax1 = plt.subplots()
        if 'northward_wind' in varkeys: # HIRLAM
            sel_var = 'northward_wind'
        elif 'air_temperature' in varkeys: # HIRLAM
            sel_var = 'air_temperature'
        elif 'so' in varkeys: # HYCOM
            sel_var = 'so'
        if 'HIRLAM' in mode:
            data_xr_check[sel_var].isel(time=0).plot(ax=ax1,x='longitude',y='latitude') #x/y are necessary since coords are not 1D and dims
        else:
            data_xr_check[sel_var].isel(time=0).sel(depth=0).plot(ax=ax1)

script_telapsed = (dt.datetime.now()-script_tstart).total_seconds()/60
print(f'elapsed time: {script_telapsed/60:.2f} hours or {script_telapsed:.2f} minutes')





