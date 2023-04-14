# -*- coding: utf-8 -*-
"""
Created on Fri Oct 14 19:58:36 2022

@author: veenstra
"""

import os
import re
from netCDF4 import Dataset
import xarray as xr
import xugrid as xu
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
    TODO: alternatively remove scale_factor key from attrs, so it can be recomputed (seems to also work): ds[var].encoding.pop('scale_factor')
    TODO: maybe add to preprocess_ERA5 (preferrably popping scale_factor attribute to keep file size small)
    """
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
                                parallel=True, #speeds up the process
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

    varkeys = data_xr.variables.mapping.keys()
    #data_xr.attrs['comment'] = 'merged with dfm_tools from https://github.com/Deltares/dfm_tools' #TODO: add something like this or other attributes? (some might also be dropped now)
    
    #select time and do checks #TODO: check if calendar is standard/gregorian
    data_xr_tsel = data_xr.sel(time=time_slice)
    if data_xr_tsel.get_index('time').duplicated().any():
        print('dropping duplicate timesteps')
        data_xr_tsel = data_xr_tsel.sel(time=~data_xr_tsel.get_index('time').duplicated()) #drop duplicate timesteps
    
    #check if there are times selected
    if len(data_xr_tsel.time)==0:
        raise OutOfRangeError(f'ERROR: no times selected, ds_text={data_xr.time[[0,-1]].to_numpy()} and time_slice={time_slice}')
    
    #check if there are no gaps (more than one unique timestep)
    times_pd = data_xr_tsel['time'].to_series()
    timesteps_uniq = times_pd.diff().iloc[1:].unique()
    if len(timesteps_uniq)>1:
        raise Exception(f'ERROR: gaps found in selected dataset (are there sourcefiles missing?), unique timesteps (hour): {timesteps_uniq/1e9/3600}')
    
    #check if requested times are available in selected files (in times_pd)
    if time_slice.start not in times_pd.index:
        raise OutOfRangeError(f'ERROR: time_slice_start="{time_slice.start}" not in selected files, timerange: "{times_pd.index[0]}" to "{times_pd.index[-1]}"')
    if time_slice.stop not in times_pd.index:
        raise OutOfRangeError(f'ERROR: time_slice_stop="{time_slice.stop}" not in selected files, timerange: "{times_pd.index[0]}" to "{times_pd.index[-1]}"')
    
    #TODO: check conversion implementation with hydro_tools\ERA5\ERA52DFM.py. Also move to separate function?
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
    for varkey_sel in ['ssr','strd']: #solar influx (surface_net_solar_radiation) / surface_thermal_radiation_downwards
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
    
    #convert 0to360 sourcedata to -180to+180
    convert_360to180 = (data_xr['longitude'].to_numpy()>180).any()
    if convert_360to180:
        data_xr.coords['longitude'] = (data_xr.coords['longitude'] + 180) % 360 - 180
        data_xr = data_xr.sortby(data_xr['longitude'])
    
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
    
    # converting from int16 to float32 (by removing dtype from encoding) or recompute scale_factor/add_offset is necessary for ERA5 dataset
    #data_xr_tsel = prevent_dtype_int(data_xr_tsel)
    data_xr_tsel = recompute_scaling_and_offset(data_xr_tsel)
    
    #data_xr_tsel.time.encoding['units'] = 'hours since 1900-01-01 00:00:00' #TODO: maybe add different reftime?
    
    return data_xr_tsel


def Dataset_varswithdim(ds,dimname): #TODO: dit zit ook in xugrid, wordt nu gebruikt in hisfile voorbeeldscript en kan handig zijn, maar misschien die uit xugrid gebruiken?
    if dimname not in ds.dims:
        raise KeyError(f'dimension {dimname} not in dataset, available are: {list(ds.dims)}')
    
    varlist_keep = []
    for varname in ds.variables.keys():
        if dimname in ds[varname].dims:
            varlist_keep.append(varname)
    ds = ds[varlist_keep]
    
    return ds


def get_vertical_dimensions(uds): #TODO: maybe add layer_dimension and interface_dimension properties to xugrid?
    """
    get vertical_dimensions from grid_info of ugrid mapfile (this will fail for hisfiles). The info is stored in the layer_dimension and interface_dimension attribute of the mesh2d variable of the dataset (stored in uds.grid after reading with xugrid)
    
    processing cb_3d_map.nc
        >> found layer/interface dimensions in file: mesh2d_nLayers mesh2d_nInterfaces
    processing Grevelingen-FM_0*_map.nc
        >> found layer/interface dimensions in file: nmesh2d_layer nmesh2d_interface (these are updated in open_partitioned_dataset)
    processing DCSM-FM_0_5nm_0*_map.nc
        >> found layer/interface dimensions in file: mesh2d_nLayers mesh2d_nInterfaces
    processing MB_02_0*_map.nc
        >> found layer/interface dimensions in file: mesh2d_nLayers mesh2d_nInterfaces
    """
    
    if not hasattr(uds,'grid'): #early return in case of e.g. hisfile
        return None, None
        
    gridname = uds.grid.name
    grid_info = uds.grid.to_dataset()[gridname]
    if hasattr(grid_info,'layer_dimension'):
        return grid_info.layer_dimension, grid_info.interface_dimension
    else:
        return None, None


def remove_ghostcells(uds): #TODO: create JIRA issue: add domainno attribute to partitioned mapfiles or remove ghostcells from output (or make values in ghostcells the same as not-ghostcells)
    """
    Dropping ghostcells if there is a domainno variable present and there is a domainno in the filename.
    Not using most-occurring domainno in var, since this is not a valid assumption for merged datasets and might be invalid for a very small partition.
    
    """
    gridname = uds.grid.name
    varn_domain = f'{gridname}_flowelem_domain'
    
    #check if dataset has domainno variable, return uds if not present
    if varn_domain not in uds.data_vars:
        print('[nodomainvar] ',end='')
        return uds
    
    #derive domainno from filename, return uds if not present
    fname = uds.encoding['source']
    if '_' not in fname: #safety escape in case there is no _ in the filename
        print('[nodomainfname] ',end='')
        return uds
    fname_splitted = fname.split('_')
    part_domainno_fromfname = fname_splitted[-2] #this is not valid for rstfiles (date follows after partnumber), but they cannot be read with xugrid anyway since they are mapformat=1
    if not part_domainno_fromfname.isnumeric() or len(part_domainno_fromfname)!=4:
        print('[nodomainfname] ',end='')
        return uds
    
    #drop ghostcells
    part_domainno_fromfname = int(part_domainno_fromfname)
    da_domainno = uds[varn_domain]
    idx = np.flatnonzero(da_domainno == part_domainno_fromfname)
    uds = uds.isel({uds.grid.face_dimension:idx})
    return uds


def remove_periodic_cells(uds): #TODO: implement proper fix: https://github.com/Deltares/xugrid/issues/63
    """
    For global models with grids that go "around the back". Temporary fix to drop all faces that are larger than grid_extent/2 (eg 360/2=180 degrees in case of GTSM)
    
    """
    #print('>> remove_periodic_cells() on dataset: ',end='')
    #dtstart = dt.datetime.now()
    face_node_x = uds.grid.face_node_coordinates[:,:,0]
    grid_extent = uds.grid.bounds[2] - uds.grid.bounds[0]
    face_node_maxdx = np.nanmax(face_node_x,axis=1) - np.nanmin(face_node_x,axis=1)
    bool_face = face_node_maxdx < grid_extent/2
    if bool_face.all(): #early return for when no cells have to be removed (might increase performance)
        return uds
    uds = uds.sel({uds.grid.face_dimension:bool_face})
    #print(f'{(dt.datetime.now()-dtstart).total_seconds():.2f} sec')
    return uds


def open_partitioned_dataset(file_nc, chunks={'time':1}, remove_ghost=True, remove_periodic=False, **kwargs): 
    """
    using xugrid to read and merge partitions, with some additional features (remaning old layerdim, timings, set zcc/zw as data_vars)

    Parameters
    ----------
    file_nc : TYPE
        DESCRIPTION.
    chunks : TYPE, optional
        chunks={'time':1} increases performance significantly upon reading, but causes memory overloads when performing sum/mean/etc actions over time dimension (in that case 100/200 is better). The default is {'time':1}.

    Raises
    ------
    Exception
        DESCRIPTION.

    Returns
    -------
    ds_merged_xu : TYPE
        DESCRIPTION.
    
    file_nc = 'p:\\1204257-dcsmzuno\\2006-2012\\3D-DCSM-FM\\A18b_ntsu1\\DFM_OUTPUT_DCSM-FM_0_5nm\\DCSM-FM_0_5nm_0*_map.nc' #3D DCSM
    file_nc = 'p:\\archivedprojects\\11206813-006-kpp2021_rmm-2d\\C_Work\\31_RMM_FMmodel\\computations\\model_setup\\run_207\\results\\RMM_dflowfm_0*_map.nc' #RMM 2D
    file_nc = 'p:\\1230882-emodnet_hrsm\\GTSMv5.0\\runs\\reference_GTSMv4.1_wiCA_2.20.06_mapformat4\\output\\gtsm_model_0*_map.nc' #GTSM 2D
    file_nc = 'p:\\11208053-005-kpp2022-rmm3d\\C_Work\\01_saltiMarlein\\RMM_2019_computations_02\\computations\\theo_03\\DFM_OUTPUT_RMM_dflowfm_2019\\RMM_dflowfm_2019_0*_map.nc' #RMM 3D
    file_nc = 'p:\\archivedprojects\\11203379-005-mwra-updated-bem\\03_model\\02_final\\A72_ntsu0_kzlb2\\DFM_OUTPUT_MB_02\\MB_02_0*_map.nc'
    Timings (xu.open_dataset/xu.merge_partitions):
        - DCSM 3D 20 partitions  367 timesteps: 231.5/ 4.5 sec (decode_times=False: 229.0 sec)
        - RMM  2D  8 partitions  421 timesteps:  55.4/ 4.4 sec (decode_times=False:  56.6 sec)
        - GTSM 2D  8 partitions  746 timesteps:  71.8/30.0 sec (decode_times=False: 204.8 sec)
        - RMM  3D 40 partitions  146 timesteps: 168.8/ 6.3 sec (decode_times=False: 158.4 sec)
        - MWRA 3D 20 partitions 2551 timesteps:  74.4/ 3.4 sec (decode_times=False:  79.0 sec)
    
    """
    #TODO: FM-mapfiles contain wgs84/projected_coordinate_system variables. xugrid has .crs property, projected_coordinate_system/wgs84 should be updated to be crs so it will be automatically handled? >> make dflowfm issue (and https://github.com/Deltares/xugrid/issues/42)
    #TODO: add support for multiple grids via keyword? GTSM+riv grid also only contains only one grid, so no testcase available
    #TODO: speed up open_dataset https://github.com/Deltares/dfm_tools/issues/225 (also remove_ghost and remove_periodic)
    
    dtstart_all = dt.datetime.now()
    file_nc_list = file_to_list(file_nc)
    
    print(f'>> xu.open_dataset() with {len(file_nc_list)} partition(s): ',end='')
    dtstart = dt.datetime.now()
    partitions = []
    for iF, file_nc_one in enumerate(file_nc_list):
        print(iF+1,end=' ')
        ds = xr.open_dataset(file_nc_one, chunks=chunks, **kwargs)
        if 'nFlowElem' in ds.dims and 'nNetElem' in ds.dims: #for mapformat1 mapfiles: merge different face dimensions (rename nFlowElem to nNetElem) to make sure the dataset topology is correct
            print('[mapformat1] ',end='')
            ds = ds.rename({'nFlowElem':'nNetElem'})
        uds = xu.core.wrap.UgridDataset(ds)
        if remove_ghost: #TODO: this makes it way slower (at least for GTSM), but is necessary since values on overlapping cells are not always identical (eg in case of Venice ucmag)
            uds = remove_ghostcells(uds)
        if remove_periodic: #TODO: makes it also slower, check if bool/idx makes difference in performance?
            uds = remove_periodic_cells(uds)
        partitions.append(uds)
    print(': ',end='')
    print(f'{(dt.datetime.now()-dtstart).total_seconds():.2f} sec')
    
    if len(partitions) == 1: #do not merge in case of 1 partition
        return partitions[0]
    
    print(f'>> xu.merge_partitions() with {len(file_nc_list)} partition(s): ',end='')
    dtstart = dt.datetime.now()
    ds_merged_xu = xu.merge_partitions(partitions)
    print(f'{(dt.datetime.now()-dtstart).total_seconds():.2f} sec')
    
    #print variables that are dropped in merging procedure. Often only ['mesh2d_face_x_bnd', 'mesh2d_face_y_bnd'], which can be derived by combining node_coordinates (mesh2d_node_x mesh2d_node_y) and face_node_connectivity (mesh2d_face_nodes). >> can be removed from FM-mapfiles (email of 16-1-2023)
    varlist_onepart = list(partitions[0].variables.keys())
    varlist_merged = list(ds_merged_xu.variables.keys())
    varlist_dropped_bool = ~pd.Series(varlist_onepart).isin(varlist_merged)
    varlist_dropped = pd.Series(varlist_onepart).loc[varlist_dropped_bool]
    if varlist_dropped_bool.any():
        print(f'>> some variables dropped with merging of partitions: {varlist_dropped.tolist()}')
    
    print(f'>> dfmt.open_partitioned_dataset() total: {(dt.datetime.now()-dtstart_all).total_seconds():.2f} sec')
    return ds_merged_xu

