# -*- coding: utf-8 -*-
"""
Created on Fri Oct 14 19:58:36 2022

@author: veenstra
"""

import os
import re
from netCDF4 import Dataset
import xarray as xr
import datetime as dt
import glob
import pandas as pd
import warnings
import numpy as np
from dfm_tools.errors import OutOfRangeError


def file_to_list(file_nc):
    if isinstance(file_nc,list):
        file_nc_list = file_nc
    else:
        basename = os.path.basename(file_nc)
        dirname = os.path.dirname(file_nc)
        ext = os.path.splitext(file_nc)[-1]
        if '.*' in basename:
            def glob_re(pattern, strings):
                return list(filter(re.compile(pattern).search, strings))
            file_nc_list = glob_re(basename, glob.glob(os.path.join(dirname,f'*{ext}')))
        else:
            file_nc_list = glob.glob(file_nc)
        file_nc_list.sort()
    if len(file_nc_list)==0:
        raise FileNotFoundError('file(s) not found, empty file_nc_list')
    return file_nc_list


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
        ds[coord] = ds[coord].load().str.decode('utf-8',errors='ignore').str.strip() #.load() is essential to convert not only first letter of string.
        ds = ds.set_index({dim:coord})
        
        #drop duplicate indices (stations/crs/gs), this avoids "InvalidIndexError: Reindexing only valid with uniquely valued Index objects"
        duplicated_keepfirst = ds[dim].to_series().duplicated(keep='first')
        if duplicated_keepfirst.sum()>0:
            print(f'dropping {duplicated_keepfirst.sum()} duplicate "{coord}" labels to avoid InvalidIndexError')
            ds = ds[{dim:~duplicated_keepfirst}]

    #check dflowfm version/date and potentially raise warning about incorrect layers
    try:
        source_attr = ds.attrs['source'] # fails if no source attr present in dataset
        source_attr_version = source_attr.split(', ')[1]
        source_attr_date = source_attr.split(', ')[2]
        if pd.Timestamp(source_attr_date) < dt.datetime(2020,11,28):
            warnings.warn(UserWarning(f'Your model was run with a D-FlowFM version from before 28-10-2020 ({source_attr_version} from {source_attr_date}), the layers in the hisfile are incorrect. Check UNST-2920 and UNST-3024 for more information, it was fixed from OSS 67858.'))
    except KeyError: #no source attr present in hisfile, cannot check version
        pass
    except IndexError: #contains no ', '
        pass

    return ds


def preprocess_hirlam(ds):
    """
    add xy variables as longitude/latitude to avoid duplicate var/dim names (we rename it anyway otherwise)
    add xy as variables again with help of NetCDF4
    this function is hopefully temporary, necessary since >1D-variables cannot have the same name as one of its dimensions in xarray. Background and future solution: https://github.com/pydata/xarray/issues/6293
    "MissingDimensionsError: 'y' has more than 1-dimension and the same name as one of its dimensions ('y', 'x'). xarray disallows such variables because they conflict with the coordinates used to label dimensions."
    might be solved in https://github.com/pydata/xarray/releases/tag/v2023.02.0 (xr.concat), but not certain
    """
    
    print('hirlam workaround: adding dropped x/y variables again as longitude/latitude')
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


def preprocess_ERA5(ds):
    """
    Reduces the expver dimension in some of the ERA5 data (mtpr and other variables), which occurs in files with very recent data. The dimension contains the unvalidated data from the latest month in the second index in the expver dimension. The reduction is done with mean, but this is arbitrary, since there is only one valid value per timestep and the other one is nan.
    """
    if 'expver' in ds.dims:
        ds = ds.mean(dim='expver')
    return ds


def preprocess_woa(ds):
    """
    WOA time units is 'months since 0000-01-01 00:00:00' and calendar is not set (360_day is the only calendar that supports that unit in xarray)
    """
    ds.time.attrs['calendar'] = '360_day'
    ds = xr.decode_cf(ds) #decode_cf after adding 360_day calendar attribute
    return ds


def prevent_dtype_int(ds): #TODO: this is not used, maybe phase out?
    """
    Prevent writing to int, since it might mess up dataset (https://github.com/Deltares/dfm_tools/issues/239 and https://github.com/pydata/xarray/issues/7039)
    Since floats are used instead of ints, the disksize of the dataset will be larger
    """
    #TODO: alternatively remove scale_factor key from attrs, so it can be recomputed (seems to also work): ds[var].encoding.pop('scale_factor')
    #TODO: maybe add to preprocess_ERA5 (preferrably popping scale_factor attribute to keep file size small)
    for var in ds.data_vars:
        var_encoding = ds[var].encoding
        if 'dtype' in var_encoding.keys():
            if 'int' in str(var_encoding['dtype']):
                ds[var].encoding.pop('dtype') #remove dtype key from attrs
    return ds


def recompute_scaling_and_offset(ds:xr.Dataset) -> xr.Dataset:
    """
    Recompute add_offset and scale_factor for variables of dtype int. As suggested by https://github.com/ArcticSnow/TopoPyScale/issues/60#issuecomment-1459747654
    This is a proper fix for https://github.com/Deltares/dfm_tools/issues/239, https://github.com/pydata/xarray/issues/7039 and more
    However, rescaling causes minor semi-accuracy loss, but it does keep the disksize small (which does not happen when converting it to dtype(float32)).

    Parameters
    ----------
    ds : xr.Dataset
        DESCRIPTION.

    Returns
    -------
    ds : xr.Dataset
        DESCRIPTION.

    """

    for var in ds.data_vars:
        da = ds[var]
        
        if 'dtype' not in da.encoding: #check if dype is available, sometime encoding={}
            continue
        dtype = da.encoding['dtype']
        
        if 'int' not in str(dtype): #skip non-int vars TODO: make proper int-check?
            continue
        if 'scale_factor' not in da.encoding.keys() or 'add_offset' not in da.encoding.keys(): #prevent rescaling of non-scaled vars (like crs) #TODO: convert to set().issubset(set())
            continue
        
        n = dtype.itemsize * 8 #n=16 for int16
        vmin = float(da.min().values)
        vmax = float(da.max().values)
        
        # stretch/compress data to the available packed range
        scale_factor = (vmax - vmin) / (2 ** n - 1)
        
        # translate the range to be symmetric about zero
        add_offset =  vmin + 2 ** (n - 1) * scale_factor
        #add_offset = (vmax+vmin)/2 #difference with above is scale_factor/2
        
        da.encoding['scale_factor'] = scale_factor
        da.encoding['add_offset'] = add_offset
    return ds


def merge_meteofiles(file_nc:str, preprocess=None, 
                     chunks:dict = {'time':1},
                     time_slice:slice = slice(None,None),
                     add_global_overlap:bool = False, zerostart:bool = False) -> xr.Dataset:
    """
    for merging for instance meteo files
    x/y and lon/lat are renamed to longitude/latitude #TODO: is this desireable?

    Parameters
    ----------
    file_nc : str
        DESCRIPTION.
    preprocess : TYPE, optional
        DESCRIPTION. The default is None.
    chunks : dict, optional
        Prevents large chunks and memory issues. The default is {'time':1}.
    time_slice : slice, optional
        DESCRIPTION. The default is slice(None,None).
    add_global_overlap : bool, optional
        GTSM specific: extend data beyond -180 to 180 longitude. The default is False.
    zerostart : bool, optional
        GTSM specific: extend data with 0-value fields 1 and 2 days before time_slice.start. The default is False.

    Raises
    ------
    Exception
        DESCRIPTION.

    Returns
    -------
    TYPE
        DESCRIPTION.

    """
    #TODO: add ERA5 conversions and features from hydro_tools\ERA5\ERA52DFM.py (except for varRhoair_alt, request FM support for varying airpressure: https://issuetracker.deltares.nl/browse/UNST-6593)
    #TODO: request FM support for charnock (etc) separate meteo forcing (currently airpressure_windx_windy_charnock merged file is required): https://issuetracker.deltares.nl/browse/UNST-6453
    #TODO: provide extfile example with fmquantity/ncvarname combinations and cleanup FM code: https://issuetracker.deltares.nl/browse/UNST-6453
    #TODO: add coordinate conversion (only valid for models with multidimensional lat/lon variables like HARMONIE and HIRLAM). This should work: ds_reproj = ds.set_crs(4326).to_crs(28992)
    #TODO: add CMCC etc from gtsmip repos (mainly calendar conversion)
    #TODO: put conversions in separate function?
    #TODO: maybe add renaming like {'salinity':'so', 'water_temp':'thetao'} for hycom
    
    #hirlam workaround
    if preprocess == preprocess_hirlam:
        drop_variables = ['x','y'] #will be added again as longitude/latitude, this is a workaround (see dfmt.preprocess_hirlam for details)
    else:
        drop_variables = None
    #woa workaround
    if preprocess == preprocess_woa:
        decode_cf = False
    else:
        decode_cf = True        

    file_nc_list = file_to_list(file_nc)

    print(f'>> opening multifile dataset of {len(file_nc_list)} files (can take a while with lots of files): ',end='')
    dtstart = dt.datetime.now()
    data_xr = xr.open_mfdataset(file_nc_list,
                                #parallel=True, #TODO: speeds up the process, but often "OSError: [Errno -51] NetCDF: Unknown file format" on WCF
                                preprocess=preprocess,
                                chunks=chunks,
                                drop_variables=drop_variables, #necessary since dims/vars with equal names are not allowed by xarray, add again later and requested matroos to adjust netcdf format.
                                decode_cf=decode_cf)
    print(f'{(dt.datetime.now()-dtstart).total_seconds():.2f} sec')
    
    #rename variables
    if 'longitude' not in data_xr.variables: #TODO: make generic, comparable rename in rename_dims_dict in dfmt.open_dataset_extra()
        if 'lon' in data_xr.variables:
            data_xr = data_xr.rename({'lon':'longitude', 'lat':'latitude'})
        elif 'x' in data_xr.variables:
            data_xr = data_xr.rename({'x':'longitude', 'y':'latitude'})
        else:
            raise KeyError('no longitude/latitude, lon/lat or x/y variables found in dataset')

    #data_xr.attrs['comment'] = 'merged with dfm_tools from https://github.com/Deltares/dfm_tools' #TODO: add something like this or other attributes? (some might also be dropped now)
    
    #select time and do checks #TODO: check if calendar is standard/gregorian
    data_xr = data_xr.sel(time=time_slice)
    if data_xr.get_index('time').duplicated().any():
        print('dropping duplicate timesteps')
        data_xr = data_xr.sel(time=~data_xr.get_index('time').duplicated()) #drop duplicate timesteps
    
    #check if there are times selected
    if len(data_xr.time)==0:
        raise OutOfRangeError(f'ERROR: no times selected, ds_text={data_xr.time[[0,-1]].to_numpy()} and time_slice={time_slice}')
    
    #check if there are no gaps (more than one unique timestep)
    times_pd = data_xr['time'].to_series()
    timesteps_uniq = times_pd.diff().iloc[1:].unique()
    if len(timesteps_uniq)>1:
        raise Exception(f'ERROR: gaps found in selected dataset (are there sourcefiles missing?), unique timesteps (hour): {timesteps_uniq/1e9/3600}')
    
    #check if requested times are available in selected files (in times_pd)
    if time_slice.start not in times_pd.index:
        raise OutOfRangeError(f'ERROR: time_slice_start="{time_slice.start}" not in selected files, timerange: "{times_pd.index[0]}" to "{times_pd.index[-1]}"')
    if time_slice.stop not in times_pd.index:
        raise OutOfRangeError(f'ERROR: time_slice_stop="{time_slice.stop}" not in selected files, timerange: "{times_pd.index[0]}" to "{times_pd.index[-1]}"')
    
    data_xr = convert_meteo_units(data_xr)
    
    #convert 0to360 sourcedata to -180to+180
    convert_360to180 = (data_xr['longitude'].to_numpy()>180).any()
    if convert_360to180: #TODO: make more flexible for models that eg pass -180/+180 crossing (add overlap at lon edges).
        lon_newvar = (data_xr.coords['longitude'] + 180) % 360 - 180
        data_xr.coords['longitude'] = lon_newvar.assign_attrs(data_xr['longitude'].attrs) #this re-adds original attrs
        data_xr = data_xr.sortby(data_xr['longitude'])
    
    #GTSM specific addition for longitude overlap
    if add_global_overlap: # assumes -180 to ~+179.75 (full global extent, but no overlap). Does not seem to mess up results for local models.
        if len(data_xr.longitude.values) != len(np.unique(data_xr.longitude.values%360)):
            raise Exception(f'add_global_overlap=True, but there are already overlapping longitude values: {data_xr.longitude}')
        overlap_ltor = data_xr.sel(longitude=data_xr.longitude<=-179)
        overlap_ltor['longitude'] = overlap_ltor['longitude'] + 360
        overlap_rtol = data_xr.sel(longitude=data_xr.longitude>=179)
        overlap_rtol['longitude'] = overlap_rtol['longitude'] - 360
        data_xr = xr.concat([data_xr,overlap_ltor,overlap_rtol],dim='longitude').sortby('longitude')
    
    #GTSM specific addition for zerovalues during spinup
    #TODO: doing this drops all encoding from variables, causing them to be converted into floats. Also makes sense since 0 pressure does not fit into int16 range as defined by scalefac and offset
    #'scale_factor': 0.17408786412952254, 'add_offset': 99637.53795606793
    #99637.53795606793 - 0.17408786412952254*32768
    #99637.53795606793 + 0.17408786412952254*32767
    if zerostart:
        field_zerostart = data_xr.isel(time=[0,0])*0 #two times first field, set values to 0
        field_zerostart['time'] = [times_pd.index[0]-dt.timedelta(days=2),times_pd.index[0]-dt.timedelta(days=1)] #TODO: is one zero field not enough? (is replacing first field not also ok? (results in 1hr transition period)
        data_xr = xr.concat([field_zerostart,data_xr],dim='time',combine_attrs='no_conflicts') #combine_attrs argument prevents attrs from being dropped
    
    # converting from int16 to float32 (by removing dtype from encoding) or recompute scale_factor/add_offset is necessary for ERA5 dataset
    #data_xr = prevent_dtype_int(data_xr)
    data_xr = recompute_scaling_and_offset(data_xr)
    
    return data_xr


def convert_meteo_units(data_xr):
    
    #TODO: check conversion implementation with hydro_tools\ERA5\ERA52DFM.py
    #TODO: keep/update attrs
    #TODO: reduce code complexity
    
    def get_unit(data_xr_var):
        if 'units' in data_xr_var.attrs.keys():
            unit = data_xr_var.attrs["units"]
        else:
            unit = '-'
        return unit
    
    varkeys = data_xr.variables.mapping.keys()
    
    #convert Kelvin to Celcius
    for varkey_sel in ['air_temperature','dew_point_temperature','d2m','t2m']: # 2 meter dewpoint temparature / 2 meter temperature
        if varkey_sel not in varkeys:
            continue
        current_unit = get_unit(data_xr[varkey_sel])
        new_unit = 'C'
        print(f'converting {varkey_sel} unit from Kelvin to Celcius: [{current_unit}] to [{new_unit}]')
        data_xr[varkey_sel] = data_xr[varkey_sel] - 273.15
        data_xr[varkey_sel].attrs['units'] = new_unit
    #convert fraction to percentage
    for varkey_sel in ['cloud_area_fraction','tcc']: #total cloud cover
        if varkey_sel not in varkeys:
            continue
        current_unit = get_unit(data_xr[varkey_sel])
        new_unit = '%' #unit is soms al %
        print(f'converting {varkey_sel} unit from fraction to percentage: [{current_unit}] to [{new_unit}]')
        data_xr[varkey_sel] = data_xr[varkey_sel] * 100
        data_xr[varkey_sel].attrs['units'] = new_unit
    #convert kg/m2/s to mm/day
    for varkey_sel in ['mer','mtpr']: #mean evaporation rate / mean total precipitation rate
        if varkey_sel not in varkeys:
            continue
        current_unit = get_unit(data_xr[varkey_sel])
        new_unit = 'mm/day'
        print(f'converting {varkey_sel} unit from kg/m2/s to mm/day: [{current_unit}] to [{new_unit}]')
        data_xr[varkey_sel] = data_xr[varkey_sel] * 86400 # kg/m2/s to mm/day (assuming rho_water=1000)
        data_xr[varkey_sel].attrs['units'] = new_unit
    #convert J/m2 to W/m2
    for varkey_sel in ['ssr','strd']: #solar influx (surface_net_solar_radiation) / surface_thermal_radiation_downwards
        if varkey_sel not in varkeys:
            continue
        current_unit = get_unit(data_xr[varkey_sel])
        new_unit = 'W m**-2'
        print(f'converting {varkey_sel} unit from J/m2 to W/m2: [{current_unit}] to [{new_unit}]')
        data_xr[varkey_sel] = data_xr[varkey_sel] / 3600 # 3600s/h #TODO: 1W = 1J/s, so does not make sense?
        data_xr[varkey_sel].attrs['units'] = new_unit
    #solar influx increase for beta=6%
    if 'ssr' in varkeys:
        print('ssr (solar influx) increase for beta=6%')
        data_xr['ssr'] = data_xr['ssr'] *.94
    
    return data_xr


def Dataset_varswithdim(ds,dimname): #TODO: dit zit ook in xugrid, wordt nu gebruikt in hisfile voorbeeldscript en kan handig zijn, maar misschien die uit xugrid gebruiken?
    """
    empty docstring
    """
    if dimname not in ds.dims:
        raise KeyError(f'dimension {dimname} not in dataset, available are: {list(ds.dims)}')
    
    varlist_keep = []
    for varname in ds.variables.keys():
        if dimname in ds[varname].dims:
            varlist_keep.append(varname)
    ds = ds[varlist_keep]
    
    return ds


