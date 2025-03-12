import os
import re
import xarray as xr
import datetime as dt
import glob
import pandas as pd
import logging
import numpy as np
from dfm_tools.interpolate_grid2bnd import _ds_sel_time_outside
from scipy.ndimage import distance_transform_edt
    
__all__ = [
    "preprocess_hisnc",
    "preprocess_ERA5",
    "preprocess_woa",
    "merge_meteofiles",
    "Dataset_varswithdim",
]

logger = logging.getLogger(__name__)


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
            # convert to string, since glob does not support pathlib.Path
            file_nc = str(file_nc)
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
            logger.warning('Your model was run with a D-FlowFM version from before 28-10-2020 '
                           f'({source_attr_version} from {source_attr_date}), the layers in the hisfile are incorrect. '
                           'Check UNST-2920 and UNST-3024 for more information, it was fixed from OSS 67858.')
    except KeyError: #no source attr present in hisfile, cannot check version
        pass
    except IndexError: #contains no ', '
        pass

    return ds


def preprocess_ERA5(ds):
    """
    Aligning ERA5 datasets before merging them. These operations are currently
    (2025) only required when (also) using previously retrieved ERA5 data.
    
    In recent datasets retrieved from ERA5 the time dimension and variable are
    now called valid_time. This is inconvenient since it causes issues when
    merging with previously retrieved datasets. However, it is not necessary
    for succesfully running a Delft3D FM simulation.
    
    Reducing the expver dimension: In the past, the expver dimension was
    present if you downloaded ERA5 data that consisted of a mix of ERA5 and
    ERA5T data. This dimension was also present in the data variables, so it
    broke code. Therefore this dimension is reduced with a mean operation.
    Any reduction operation would do the trick since there is only one valid
    value per timestep and the other one is nan. In datasets downloaded
    currently (2025) the expver dimension is not present anymore,
    but anexpver variable is present defining whether the data comes
    from ERA5 (1) or ERA5T (5).
    
    Adding expver coordinate if missing: The old datafiles did not contain an
    expver variable (sometimes did contain an expver dim). The new datafiles
    do contain an expver coordinate variable. Merging old and new files is only
    possible if the coordinates are the same, so add a expver coordinate
    variable to the old files with empty values.
    
    Removing scale_factor and add_offset: In the past, the ERA5 data was
    supplied as integers with a scaling and offset that was different for
    each downloaded file. This caused serious issues with merging files,
    since the scaling/offset from the first file was assumed to be valid
    for the others also, leading to invalid values. Only relevant for old
    files. More info at https://github.com/Deltares/dfm_tools/issues/239.
    """
    
    # datasets retrieved with new CDS have valid_time instead of time dim/var
    # https://forum.ecmwf.int/t/new-time-format-in-era5-netcdf-files/3796/5
    if 'valid_time' in ds.coords:
        ds = ds.rename({'valid_time':'time'})
    
    # datasets retrieved from feb 2025 onwards have different mer/mtpr varnames
    # convert back for backwards compatibility and clarity
    # https://github.com/Deltares/dfm_tools/issues/1140
    if 'avg_tprate' in ds.data_vars:
        ds = ds.rename_vars({'avg_tprate':'mtpr'})
    if 'avg_ie' in ds.data_vars:
        ds = ds.rename_vars({'avg_ie':'mer'})
    
    # reduce the expver dimension (not present in newly retrieved files)
    if 'expver' in ds.dims:
        ds = ds.mean(dim='expver')
    
    # add empty expver coordinate to old files if not present to prevent
    # "ValueError: coordinate 'expver' not present in all datasets"
    # when merging old datasets (without expver coord) with new datasets
    # has to be <U4 to avoid "NotImplementedError: Can not use auto rechunking
    # with object dtype"
    if 'expver' not in ds.variables:
        data_expver = np.empty(shape=(len(ds.time)), dtype='<U4')
        ds['expver'] = xr.DataArray(data=data_expver, dims='time')
        ds = ds.set_coords('expver')
    
    # drop scaling/offset encoding if present and converting to float32. Not
    # present in newly retrieved files, variables are zipped float32 instead
    for var in ds.data_vars.keys():
        list_attrs = ['dtype','scale_factor','add_offset']
        if not set(list_attrs).issubset(ds.variables[var].encoding.keys()):
            continue
        # the _FillValue will still be -32767 (int default)
        # this is no issue for float32
        ds[var].encoding.pop('scale_factor')
        ds[var].encoding.pop('add_offset')
        ds[var].encoding["dtype"] = "float32"
    
    return ds


def preprocess_woa(ds):
    """
    WOA time units is 'months since 0000-01-01 00:00:00' and calendar is not set (360_day is the only calendar that supports that unit in xarray)
    """
    ds.time.attrs['calendar'] = '360_day'
    ds = xr.decode_cf(ds) #decode_cf after adding 360_day calendar attribute
    return ds


def merge_meteofiles(file_nc:str,
                     time_slice:slice,
                     preprocess = None,
                     **kwargs) -> xr.Dataset:
    """
    Merging of meteo files. Variables/coordinates x/y and lon/lat are renamed
    to longitude/latitude.

    Parameters
    ----------
    file_nc : str
        DESCRIPTION.
    preprocess : TYPE, optional
        DESCRIPTION. The default is None.
    time_slice : slice
        slice(tstart,tstop).
    kwargs : dict, optional
        arguments for xr.open_mfdataset() like `chunks` to prevent large chunks and resulting memory issues.

    Returns
    -------
    data_xr : xr.Dataset
        Merged meteo dataset.

    """
    #TODO: add ERA5 conversions and features from hydro_tools\ERA5\ERA52DFM.py (except for varRhoair_alt, request FM support for varying airpressure: https://issuetracker.deltares.nl/browse/UNST-6593)
    #TODO: add coordinate conversion (only valid for models with multidimensional lat/lon variables like HARMONIE and HIRLAM). This should work: ds_reproj = ds.set_crs(4326).to_crs(28992)
    #TODO: add CMCC etc from gtsmip repos (mainly calendar conversion)
    #TODO: maybe add renaming like {'salinity':'so', 'water_temp':'thetao'} for hycom
       
    # woa workaround
    if preprocess == preprocess_woa:
        decode_cf = False
    else:
        decode_cf = True        

    if 'chunks' not in kwargs:
        # enable dask chunking
        kwargs['chunks'] = 'auto'
    if 'data_vars' not in kwargs:
        # avoid time dimension on other variables
        # enforce error in case of conflicting variables
        kwargs['data_vars'] = 'minimal'
    if 'coords' not in kwargs:
        # support number coordinate variable in some of the ERA5 datasets
        kwargs['coords'] = 'minimal'
    if 'join' not in kwargs:
        # forbid slightly changed lat/lon values
        # enforce alignment error if expver is not present in all datasets 
        kwargs['join'] = 'exact'

    file_nc_list = file_to_list(file_nc)
    print(f'>> opening multifile dataset of {len(file_nc_list)} files (can take a while with lots of files): ',end='')
    dtstart = dt.datetime.now()
    data_xr = xr.open_mfdataset(file_nc_list,
                                preprocess=preprocess,
                                decode_cf=decode_cf,
                                **kwargs)
    print(f'{(dt.datetime.now()-dtstart).total_seconds():.2f} sec')
    
    # rename variables
    # TODO: make generic, comparable rename in rename_dims_dict in dfmt.interpolate_grid2bnd.open_prepare_dataset()
    if 'longitude' not in data_xr.variables:
        if 'lon' in data_xr.variables:
            data_xr = data_xr.rename({'lon':'longitude', 'lat':'latitude'})
        elif 'x' in data_xr.variables:
            data_xr = data_xr.rename({'x':'longitude', 'y':'latitude'})
        else:
            raise KeyError('no longitude/latitude, lon/lat or x/y variables found in dataset')
    
    # check for duplicated timesteps
    if data_xr.get_index('time').duplicated().any():
        print('dropping duplicate timesteps')
        # drop duplicate timesteps
        data_xr = data_xr.sel(time=~data_xr.get_index('time').duplicated())
    
    # TODO: check if calendar is standard/gregorian
    # check available times and select outside bounds
    data_xr = _ds_sel_time_outside(
        data_xr,
        tstart=time_slice.start,
        tstop=time_slice.stop,
        )
    
    #check if there are no gaps (more than one unique timestep)
    times_pd = data_xr['time'].to_series()
    timesteps_uniq = times_pd.diff().iloc[1:].unique()
    if len(timesteps_uniq)>1:
        raise ValueError(
            'time gaps found in selected dataset (missing files?), '
            f'unique timesteps (hour): {timesteps_uniq/1e9/3600}'
            )
    
    data_xr = convert_meteo_units(data_xr)
    
    #convert 0to360 sourcedata to -180to+180
    convert_360to180 = (data_xr['longitude'].to_numpy()>180).any()
    if convert_360to180: #TODO: make more flexible for models that eg pass -180/+180 crossing (add overlap at lon edges).
        lon_newvar = (data_xr.coords['longitude'] + 180) % 360 - 180
        data_xr.coords['longitude'] = lon_newvar.assign_attrs(data_xr['longitude'].attrs) #this re-adds original attrs
        data_xr = data_xr.sortby(data_xr['longitude'])

    return data_xr


def convert_meteo_units(data_xr):
    #TODO: check conversion implementation with hydro_tools\ERA5\ERA52DFM.py
    #TODO: assert old unit instead of always converting
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
    #solar influx increase for beta=6% subtraction in DFM
    if 'ssr' in varkeys:
        print('ssr (solar influx) increase for beta=6% subtraction in DflowFM')
        data_xr['ssr'] = data_xr['ssr'] / 0.94
    
    return data_xr


def Dataset_varswithdim(ds,dimname): #TODO: dit zit ook in xugrid, wordt nu gebruikt in hisfile voorbeeldscript en kan handig zijn, maar misschien die uit xugrid gebruiken?
    """
    empty docstring
    """
    if dimname not in ds.dims:
        raise KeyError(f'dimension {dimname} not in dataset, available are: {list(ds.dims)}')
    
    varlist_keep = []
    for varname in ds.variables.keys():
        if dimname in ds.variables[varname].dims:
            varlist_keep.append(varname)
    ds = ds[varlist_keep]
    
    return ds


def _nearest(a):
    nans = np.isnan(a)
    if not nans.any():
        return a.copy()
    indices = distance_transform_edt(
        input=np.isnan(a),
        return_distances=False,
        return_indices=True,
    )
    return a[tuple(indices)]


def interpolate_na_multidim(da, dim, keep_attrs=True):
    """
    Interpolate_na for multiple dimensions at once. Since it 
    """
    arr = xr.apply_ufunc(
        _nearest,
        da,
        input_core_dims=[dim],
        output_core_dims=[dim],
        output_dtypes=[da.dtype],
        dask="parallelized",
        vectorize=True,
        keep_attrs=keep_attrs,
    ).transpose(*da.dims)
    return arr
