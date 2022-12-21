# -*- coding: utf-8 -*-
"""
Created on Wed Oct 27 22:07:34 2021

@author: veenstra

script to merge netcdf file from several folders into one file (concatenate time dimension)
too slow? Change fn_match_pattern to match less files if possible, or provide a list of (less) files as file_nc

"""

import glob
#import re
import datetime as dt
import os
import xarray as xr
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
plt.close('all')
import dfm_tools as dfmt

#TODO: crashes in pytest for some reason: "OSError: [Errno -51] NetCDF: Unknown file format"
#TODO: add ERA5 conversions and features from hydro_tools\ERA5\ERA52DFM.py (except for varRhoair_alt, request FM support of varying airpressure)
#TODO: request FM support for charnock (etc) separate meteo forcing (now merged files are required)
#TODO: add standard_name to all variables? >> better to rename varname itself to dfm quantity (not possible yet) or alternative dfm varname (or generate ext file with quantity/varname translation)
#TODO: add renamevars of add attrs['standard_name'] with ncvarnames/ncstdnames conversion table in https://svn.oss.deltares.nl/repos/delft3d/trunk/src/utils_lgpl/ec_module/packages/ec_module/src/ec_provider.f90 (line 2479) (alternatively set varname in extfile)
#varname_dict = {'u10': '',
#                'v10': '',
#                'msl': ''}
#stdname_dict = {'u10': 'eastward_wind',
#                'v10': 'northward_wind',
#                'msl': 'air_pressure'}
#TODO: add coordinate conversion (maybe only for models with multidimensional lat/lon variables like HARMONIE and HIRLAM). This should work: ds_reproj = ds.set_crs(4326).to_crs(28992)
#TODO: add CMCC etc from gtsmip repos (mainly calendar conversion)
#TODO: add convert_360to180 with "ds.coords['lon'] = (ds.coords['lon'] + 180) % 360 - 180; ds = ds.sortby('lon')" (without hardcoded longitude/lon names)
#TODO: move to function

add_global_overlap = False #GTSM specific: extend data beyond -180 to 180 longitude
zerostart = False #GTSM specific: extend data with 0-value fields 1 and 2 days before all_tstart

mode = 'ERA5_wind_pressure' # 'HARMONIE' 'HIRLAM_meteo' 'HIRLAM_meteo-heatflux' 'HYCOM' 'ERA5_wind_pressure' 'ERA5_heat_model' 'ERA5_radiation' 'ERA5_rainfall'
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
    preprocess = dfmt.preprocess_hirlam
    rename_variables = None
elif mode == 'HARMONIE':
    dir_data = 'p:\\1204257-dcsmzuno\\data\\meteo\\HARMONIE\\nc\\air_*' #many invalid files, so subsetting here
    fn_match_pattern = 'HARMONIE_*_2020_*.nc'
    file_out_prefix = fn_match_pattern.replace('*_2020_*.nc','')
    drop_variables = None
    preprocess = None
    rename_variables = {'x':'longitude', 'y':'latitude'}
    #all_tstart = dt.datetime(2020,1,1)
    #all_tstop = dt.datetime(2020,6,1)
elif mode == 'HYCOM':
    dir_data = 'c:\\DATA\\dfm_tools_testdata\\GLBu0.08_expt_91.2'
    fn_match_pattern = 'HYCOM_ST_GoO_*.nc'
    file_out_prefix = fn_match_pattern.replace('*.nc','')
    drop_variables = None
    preprocess = None
    rename_variables = {'salinity':'so', 'water_temp':'thetao'}
elif 'ERA5' in mode:
    #all_tstart = dt.datetime(2005,1,1) #for performance checking, was 12 minutes
    #all_tstop = dt.datetime(2022,1,1)
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
    if 1:
        dir_data = 'p:\\metocean-data\\open\\ERA5\\data\\Irish_North_Baltic_Sea\\*' #TODO: add other vars with * (support separate files)
        fn_match_pattern = f'era5_.*({"|".join(varkey_list)})_.*\.nc'
        file_out_prefix = f'era5_{"_".join(varkey_list)}'
    else: #global test (first use download_ERA5.py)
        dir_data = r'c:\DATA\dfm_tools\tests\examples\v10n'
        fn_match_pattern = 'era5_*.nc'
        file_out_prefix = fn_match_pattern.replace('*.nc','')
        all_tstart = dt.datetime(2021,1,1)
        all_tstop = dt.datetime(2021,2,28,23,0)
    drop_variables = None
    preprocess = None
    rename_variables = None
elif mode == 'WOA': #TODO: does not work since time units is 'months since 0000-01-01 00:00:00' and calendar is not set (360_day is the only one that supports that unit)
    def preprocess_woa(ds):
        ds.time.attrs['calendar'] = '360_day'
        return ds
    dir_data = r'p:\1204257-dcsmzuno\data\WOA13'
    fn_match_pattern = 'woa13_decav_s*.nc'
    file_out_prefix = fn_match_pattern.replace('*.nc','')
    drop_variables = None
    preprocess = preprocess_woa
    rename_variables = None
else:
    raise Exception('ERROR: wrong mode %s'%(mode))

dir_output = '.' #os.path.join(dir_data,'merged_netcdf_files_JV')
if not os.path.exists(dir_output):
    os.makedirs(dir_output)


#generating file list
if '.*' in fn_match_pattern: #regex complex stuff #TODO: make cleaner, but watch v10/v10n matching
    file_pd_allnc = pd.Series(glob.glob(os.path.join(dir_data,'*.nc'))) #all nc files
    file_pd_allnc_varnames = file_pd_allnc.str.extract(fn_match_pattern)[0].str.replace('_','')
    file_pd_allnc_df = pd.DataFrame({'file':file_pd_allnc,'varname':file_pd_allnc_varnames})
    file_pd_pattern = file_pd_allnc_df.loc[~file_pd_allnc_df['varname'].isnull()]
    file_nc = file_pd_pattern['file'].tolist()
    file_list = file_nc
else:
    file_nc = os.path.join(dir_data,fn_match_pattern)
    file_list = glob.glob(file_nc)

print(f'opening multifile dataset of {len(file_list)} files matching "{fn_match_pattern}" (can take a while with lots of files)')
#open_mfdataset when dropping x/y since both varname as dimname, which is not possible in xarray
data_xr = xr.open_mfdataset(file_nc,
                            drop_variables=drop_variables, #necessary since dims/vars with equal names are not allowed by xarray, add again later and requested matroos to adjust netcdf format.
                            parallel=True, #speeds up the process
                            preprocess=preprocess,
                            chunks={'time':1}, #this prevents large chunks and memory issues
                            #concat_dim="time", combine="nested", data_vars='minimal', coords='minimal', compat='override', #TODO: optional vars to look into: https://docs.xarray.dev/en/stable/user-guide/io.html#reading-multi-file-datasets. might also resolve large chunks warning with ERA5 (which has dissapeared somehow)
                            )
print('...done')


#rename variables
data_xr = data_xr.rename(rename_variables)
varkeys = data_xr.variables.mapping.keys()
#data_xr.attrs['comment'] = 'merged with dfm_tools from https://github.com/openearth/dfm_tools' #TODO: add something like this or other attributes? (some might also be dropped now)

#select time and do checks #TODO: check if calendar is standard/gregorian
print('time selection')
data_xr_tsel = data_xr.sel(time=slice(all_tstart,all_tstop))
if data_xr_tsel.get_index('time').duplicated().any():
    print('dropping duplicate timesteps')
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

#TODO: check conversion implementation with hydro_tools\ERA5\ERA52DFM.py
def get_unit(data_xr_var):
    if 'units' in data_xr_var.attrs.keys():
        unit = data_xr_var.attrs["units"]
    else:
        unit = '-'
    return unit
#convert Kelvin to Celcius
for varkey_sel in ['air_temperature','dew_point_temperature','d2m','t2m']: # 2 meter dewpoint temparature / 2 meter temperature
    if varkey_sel in varkeys:
        current_unit = get_unit(data_xr_tsel[varkey_sel])
        new_unit = 'C'
        print(f'converting {varkey_sel} unit from Kelvin to Celcius: [{current_unit}] to [{new_unit}]')
        data_xr_tsel[varkey_sel].attrs['units'] = new_unit
        data_xr_tsel[varkey_sel] = data_xr_tsel[varkey_sel] - 273.15
#convert fraction to percentage
for varkey_sel in ['cloud_area_fraction','tcc']: #total cloud cover
    if varkey_sel in varkeys:
        current_unit = get_unit(data_xr_tsel[varkey_sel])
        new_unit = '%' #unit is soms al %
        print(f'converting {varkey_sel} unit from fraction to percentage: [{current_unit}] to [{new_unit}]')
        data_xr_tsel[varkey_sel].attrs['units'] = new_unit
        data_xr_tsel[varkey_sel] = data_xr_tsel[varkey_sel] * 100
#convert kg/m2/s to mm/day
for varkey_sel in ['mer','mtpr']: #mean evaporation rate / mean total precipitation rate
    if varkey_sel in varkeys:
        current_unit = get_unit(data_xr_tsel[varkey_sel])
        new_unit = 'mm/day'
        print(f'converting {varkey_sel} unit from kg/m2/s to mm/day: [{current_unit}] to [{new_unit}]')
        data_xr_tsel[varkey_sel].attrs['units'] = new_unit
        data_xr_tsel[varkey_sel] = data_xr_tsel[varkey_sel] * 86400 # kg/m2/s to mm/day (assuming rho_water=1000)
#convert J/m2 to W/m2
for varkey_sel in ['ssr','strd']: #solar influx (surface_net_solar_radiation) / surface_thermal_radiation_downwards #TODO: 
    if varkey_sel in varkeys:
        current_unit = get_unit(data_xr_tsel[varkey_sel])
        new_unit = 'W m**-2'
        print(f'converting {varkey_sel} unit from J/m2 to W/m2: [{current_unit}] to [{new_unit}]')
        data_xr_tsel[varkey_sel].attrs['units'] = new_unit
        data_xr_tsel[varkey_sel] = data_xr_tsel[varkey_sel] / 3600 # 3600s/h #TODO: 1W = 1J/s, so does not make sense?
#solar influx increase for beta=6%
if 'ssr' in varkeys:
    print('ssr (solar influx) increase for beta=6%')
    data_xr_tsel['ssr'] = data_xr_tsel['ssr'] *.94


#GTSM specific addition for longitude overlap
if add_global_overlap: # assumes -180 to ~+179.75 (full global extent, but no overlap). Does not seem to mess up results for local models.
    if len(data_xr_tsel.longitude.values) != len(np.unique(data_xr_tsel.longitude.values%360)):
        raise Exception(f'add_global_overlap=True, but there are already overlapping longitude values: {data_xr_tsel.longitude}')
    overlap_ltor = data_xr_tsel.sel(longitude=data_xr_tsel.longitude<=-179)
    overlap_ltor['longitude'] = overlap_ltor['longitude'] + 360
    overlap_rtol = data_xr_tsel.sel(longitude=data_xr_tsel.longitude>=179)
    overlap_rtol['longitude'] = overlap_rtol['longitude'] - 360
    data_xr_tsel = xr.concat([data_xr_tsel,overlap_ltor,overlap_rtol],dim='longitude').sortby('longitude')

#GTSM specific addition for zerovalues during spinup
if zerostart:
    field_zerostart = data_xr_tsel.isel(time=[0,0])*0 #two times first field, set values to 0
    field_zerostart['time'] = [times_pd.index[0]-dt.timedelta(days=2),times_pd.index[0]-dt.timedelta(days=1)] #TODO: is one zero field not enough? (is replacing first field not also ok? (results in 1hr transition period)
    data_xr_tsel = xr.concat([field_zerostart,data_xr_tsel],dim='time')#.sortby('time')

encoding = {}
#encoding['time'] = {'units': 'hours since 1900-01-01 00:00:00'} #TODO: maybe add different reftime?
#for varkey in list(data_xr_tsel.data_vars.keys()):
#    encoding[varkey] = {'scale_factor':0.01,'add_offset':0} #TODO: maybe add, but not necessary since xarray uses encoding from first file and that is already quite efficient.

# ERA5 specific addition for combine validated and non-validated data
if 'expver' in data_xr_tsel.dims:
    data_xr_tsel = data_xr_tsel.mean(dim='expver')

#write to netcdf file
print('writing file (can take a while)')
file_out = os.path.join(dir_output, f'{file_out_prefix}_{tstart_str}to{tstop_str}_{mode}.nc')
data_xr_tsel.to_netcdf(file_out, encoding=encoding)

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
        elif 'expver' in data_xr_tsel[varkey].coords:
            data_xr_check[varkey].isel(time=0).isel(expver=0).plot(ax=ax1)
        else:
            data_xr_check[varkey].isel(time=0).plot(ax=ax1)
        fig.savefig(file_out.replace('.nc',f'_{varkey}'))

script_telapsed = (dt.datetime.now()-script_tstart)
print(f'elapsed time: {script_telapsed}')





