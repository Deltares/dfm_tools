# -*- coding: utf-8 -*-
"""
Created on Thu Aug 18 17:39:03 2022

@author: veenstra

"""

import os
import glob
import datetime as dt
import numpy as np
import pandas as pd
import xarray as xr
import xugrid as xu
from pathlib import Path
from scipy.spatial import KDTree
import warnings
import hydrolib.core.dflowfm as hcdfm
import geopandas

from dfm_tools.hydrolib_helpers import Dataset_to_TimeSeries, Dataset_to_T3D, Dataset_to_Astronomic, pointlike_to_DataFrame, PolyFile_to_geodataframe_points
from dfm_tools.errors import OutOfRangeError


def get_conversion_dict(ncvarname_updates={}):
    
    """
    get the conversion_dict, optionally with updated ncvarnames
    conversion_dict.keys() are the dflowfm quantities and the ncvarname the corresponding netcdf variable name/key
    alternative ncvarnames can be supplied via ncvarname_updates, e.g. get_conversion_dict(ncvarname_updates={'temperaturebnd':'tos'})
    
    interpolate_nc_to_bc() renames netcdf variable like this:
    data_xr = data_xr.rename({ncvarname:quantity})
    
    for CMCC:
    conversion_dict = { # mg/l is the same as g/m3: conversion is phyc in mol/m3 to newvar in g/m3
                       'tracerbndOXY'        : {'ncvarname': 'o2',          'unit': 'g/m3', 'conversion': 32.0 },
                       'tracerbndNO3'        : {'ncvarname': 'no3',         'unit': 'g/m3', 'conversion': 14.0 },
                       'tracerbndPO4'        : {'ncvarname': 'po4',         'unit': 'g/m3', 'conversion': 30.97 },
                       'tracerbndSi'         : {'ncvarname': 'si',          'unit': 'g/m3', 'conversion': 28.08},
                       'tracerbndPON1'       : {'ncvarname': 'phyc',        'unit': 'g/m3', 'conversion': 14.0}, # Caution: this empirical relation might not be applicable to your use case
                       'tracerbndPOP1'       : {'ncvarname': 'phyc',        'unit': 'g/m3', 'conversion': 30.97}, # Caution: this empirical relation might not be applicable to your use case
                       'tracerbndPOC1'       : {'ncvarname': 'phyc',        'unit': 'g/m3', 'conversion': 14.0 * (106 / 16)}, # Caution: this empirical relation might not be applicable to your use case
                       'salinitybnd'         : {'ncvarname': 'sos'},         #'1e-3'
                       'temperaturebnd'      : {'ncvarname': 'tos'},         #'degC'
                       'ux'                  : {'ncvarname': 'uo'},          #'m/s'
                       'uy'                  : {'ncvarname': 'vo'},          #'m/s'
                       'waterlevelbnd'       : {'ncvarname': 'zos'},         #'m' #steric
                       'tide'                : {'ncvarname': ''},            #'m' #tide (dummy entry)
                       }
    """
    # conversion_dict, 
    conversion_dict = { # mg/l is the same as g/m3: conversion is phyc in mmol/l to newvar in g/m3
                        'tracerbndOXY'        : {'ncvarname': 'o2',          'unit': 'g/m3', 'conversion': 32.0 / 1000.0}, 
                        'tracerbndNO3'        : {'ncvarname': 'no3',         'unit': 'g/m3', 'conversion': 14.0 / 1000.0},
                        'tracerbndPO4'        : {'ncvarname': 'po4',         'unit': 'g/m3', 'conversion': 30.97 / 1000.0},
                        'tracerbndSi'         : {'ncvarname': 'si',          'unit': 'g/m3', 'conversion': 28.08 / 1000.0},
                        'tracerbndPON1'       : {'ncvarname': 'phyc',        'unit': 'g/m3', 'conversion': 2. * 16. * 14. / (106. * 1000.0)}, # Caution: this empirical relation might not be applicable to your use case
                        'tracerbndPOP1'       : {'ncvarname': 'phyc',        'unit': 'g/m3', 'conversion': 2. * 30.97 / (106. * 1000.0)}, # Caution: this empirical relation might not be applicable to your use case
                        'tracerbndPOC1'       : {'ncvarname': 'phyc',        'unit': 'g/m3', 'conversion': 2. * 12. / 1000.0}, # Caution: this empirical relation might not be applicable to your use case
                        'tracerbndDON'        : {'ncvarname': 'phyc',        'unit': 'g/m3', 'conversion': 3.24 * 2. * 16. * 14. / (106. * 1000.0)}, # Caution: this empirical relation might not be applicable to your use case
                        'tracerbndDOP'        : {'ncvarname': 'phyc',        'unit': 'g/m3', 'conversion': 1.0 * 2. * 30.97 / (106. * 1000.0)}, # Caution: this empirical relation might not be applicable to your use case
                        'tracerbndDOC'        : {'ncvarname': 'phyc',        'unit': 'g/m3', 'conversion': (199. / 20.) * 3.24 * 2. * 16. * 12. / (106. * 1000.0)}, # Caution: this empirical relation might not be applicable to your use case
                        'tracerbndOpal'       : {'ncvarname': 'phyc',        'unit': 'g/m3', 'conversion': 0.5 * 0.13 * 28.08 / (1000.0)}, # Caution: this empirical relation might not be applicable to your use case
                        'salinitybnd'         : {'ncvarname': 'so'},          #'1e-3'
                        'temperaturebnd'      : {'ncvarname': 'thetao'},      #'degC'
                        'ux'                  : {'ncvarname': 'uo'},          #'m/s'
                        'uy'                  : {'ncvarname': 'vo'},          #'m/s'
                        'waterlevelbnd'       : {'ncvarname': 'zos'},         #'m' #steric
                        'tide'                : {'ncvarname': ''},            #'m' #tide (dummy entry)
                        }
    #do updates
    for k,v in ncvarname_updates.items():
        conversion_dict[k]['ncvarname'] = v
    return conversion_dict


def interpolate_tide_to_bc(tidemodel, file_pli, component_list=None, nPoints=None):
    """
    empty docstring
    """
    # translate dict from .\hydro_tools\FES\PreProcessing_FES_TideModel_imaginary.m
    #component_list = ['2N2','LABDA2','MF','MFM','P1','SSA','EPSILON2','M2','MKS2','MU2','Q1','T2','J1','M3','MM','N2','R2','K1','M4','MN4','N4','S1','K2','M6','MS4','NU2','S2','L2','M8','MSF','O1','S4','MSQM','SA']
    translate_dict = {'LA2':'LABDA2', #TODO: use value instead of key in bc file? Support using value instead of key in const_list also (like line above)
                      'EPS2':'EPSILON2', 
                      'Z0':'A0',
                      'MTM':'MFM', #Needs to be verified
                      }
    
    dir_pattern_dict = {'FES2014': Path(r'P:\metocean-data\licensed\FES2014','*.nc'), #ocean_tide_extrapolated
                        'FES2012': Path(r'P:\metocean-data\open\FES2012\data','*_FES2012_SLEV.nc'), #is eigenlijk ook licensed
                        'EOT20': Path(r'P:\metocean-data\open\EOT20\ocean_tides','*_ocean_eot20.nc'),
                        'GTSM4.1preliminary': Path(r'p:\1230882-emodnet_hrsm\GTSMv3.0EMODnet\EMOD_MichaelTUM_yearcomponents\GTSMv4.1_yeartide_2014_2.20.06\compare_fouhis_fouxyz_v3','gtsmv4.1_2014_*_withfu_v3_rasterized.nc'),
                        #TODO: add tpxo8 (tpxo9 is also available), catalog: https://opendap.deltares.nl/thredds/catalog/opendap/deltares/delftdashboard/tidemodels/tpxo80/catalog.html
                        }
    if tidemodel not in dir_pattern_dict.keys():
        raise KeyError(f'invalid tidemodel "{tidemodel}", options are: {list(dir_pattern_dict.keys())}')
    if tidemodel == 'GTSM4.1preliminary':
        warnings.warn(UserWarning(f'you are using tidemodel "{tidemodel}", beware that the dataset is preliminary so it is still quite coarse and may contain errors. Check your results carefully'))
    
    #Check whether the polyfile contains multiple polyline, in that case show a warning
    pli = hcdfm.PolyFile(file_pli)
    if len(pli.objects) > 1:
        warnings.warn(UserWarning(f"The polyfile {file_pli} contains multiple polylines. Only the first one will be used by DFLOW-FM for the boundary conditions."))
        #TODO when issue UNST-7012 is solved, remove this warning or add it in more places)
    
    dir_pattern = dir_pattern_dict[tidemodel]
    
    if component_list is None:
        file_list_nc = glob.glob(str(dir_pattern))
        dir_pattern_basename = os.path.basename(dir_pattern)
        replace = dir_pattern_basename.split('*')
        component_list = [os.path.basename(x).replace(replace[0],'').replace(replace[1],'') for x in file_list_nc] #TODO: make this less hard-coded
    component_list_upper_pd = pd.Series([x.upper() for x in component_list]).replace(translate_dict, regex=True)
    
    def extract_component(ds):
        #https://github.com/pydata/xarray/issues/1380
        if 'FES2012' in ds.encoding["source"]: #TODO: make more generic with regex, or just add tidemodel argument since they are quite specific
            compname = os.path.basename(ds.encoding["source"]).replace('_FES2012_SLEV.nc','')
            ds = ds.sel(lon=ds.lon<360) #drop last instance, since 0 and 360 are both present
            ds = ds.rename({'Ha':'amplitude','Hg':'phase'})
        elif 'eot20' in ds.encoding["source"]:
            compname = os.path.basename(ds.encoding["source"]).replace('_ocean_eot20.nc','')
            ds = ds.sel(lon=ds.lon<360) #drop last instance, since 0 and 360 are both present
        elif 'gtsm' in ds.encoding["source"].lower(): #TODO: make less hard coded
            compname = os.path.basename(ds.encoding["source"]).replace('gtsmv4.1_2014_','').replace('_withfu_v3_rasterized.nc','')
            ds = ds.rename({f'wl_amp{compname}':'amplitude',f'wl_phs{compname}':'phase'}) #TODO: adjust in rasterized dataset?
            ds['amplitude'] = ds['amplitude'].assign_attrs({'units':'m'}) #TODO: handle this rasterize function
            ds['phase'] = ds['phase'].assign_attrs({'units':'degrees'}) #TODO: handle this rasterize function
            ds = ds.rename({'x':'lon','y':'lat'})
        else:
            compname = os.path.basename(ds.encoding["source"]).replace('.nc','')
        compnumber = [component_list.index(compname)]
        ds = ds.assign(compno=compnumber)
        
        convert_360to180 = (ds['lon'].to_numpy()>180).any()
        if convert_360to180: # results in large chunks if it is done after concatenation, so do for each file before concatenation
            ds.coords['lon'] = (ds.coords['lon'] + 180) % 360 - 180
            ds = ds.sortby('lon')
        return ds
    
    #use open_mfdataset() with preprocess argument to open all requested FES files into one Dataset
    file_list_nc = [str(dir_pattern).replace('*',comp) for comp in component_list]
    data_xrsel = xr.open_mfdataset(file_list_nc, combine='nested', concat_dim='compno', preprocess=extract_component)
    data_xrsel = data_xrsel.rename({'lon':'longitude','lat':'latitude'})
    
    #derive uv phase components (using amplitude=1)
    data_xrsel_phs_rad = np.deg2rad(data_xrsel['phase'])
    #we need to compute u/v components for the phase to avoid zero-crossing interpolation issues
    data_xrsel['phase_u'] = 1*np.cos(data_xrsel_phs_rad)
    data_xrsel['phase_v'] = 1*np.sin(data_xrsel_phs_rad)
    data_xrsel['compnames'] = xr.DataArray(component_list_upper_pd,dims=('compno')) #TODO: convert to proper string variable
    
    #convert cm to m
    if data_xrsel['amplitude'].attrs['units'] == 'cm':
        data_xrsel['amplitude'] /= 100
        data_xrsel['amplitude'].attrs['units'] = 'm'
    
    data_interp = interp_regularnc_to_plipoints(data_xr_reg=data_xrsel, file_pli=file_pli, nPoints=nPoints)
    data_interp['phase_new'] = np.rad2deg(np.arctan2(data_interp['phase_v'],data_interp['phase_u']))
    
    ForcingModel_object = plipointsDataset_to_ForcingModel(plipointsDataset=data_interp)
    
    return ForcingModel_object


def open_dataset_extra(dir_pattern, quantity, tstart, tstop, conversion_dict=None, refdate_str=None, reverse_depth=False, chunks=None):
    """
    empty docstring
    """
    
    if conversion_dict is None:
        conversion_dict = get_conversion_dict()
    
    if quantity=='uxuyadvectionvelocitybnd': #T3Dvector
        quantity_list = ['ux','uy']
    else:
        quantity_list = [quantity]
    ncvarname_list = [conversion_dict[quan]['ncvarname'] for quan in quantity_list]
    
    #convert tstart/tstop from str/dt.datetime/pd.Timestamp to pd.Timestamp. WARNING: when supplying '05-04-2016', monthfirst is assumed, so 2016-05-04 will be the result (may instead of april).
    tstart = pd.Timestamp(tstart)
    tstop = pd.Timestamp(tstop)
    
    dir_pattern = [Path(str(dir_pattern).format(ncvarname=ncvarname)) for ncvarname in ncvarname_list]
    file_list_nc = []
    for dir_pattern_one in dir_pattern:
        file_list_nc = file_list_nc + glob.glob(str(dir_pattern_one))
    list_pattern_names = [x.name for x in dir_pattern]
    print(f'loading mfdataset of {len(file_list_nc)} files with pattern(s) {list_pattern_names}')
    
    try:
        data_xr = xr.open_mfdataset(file_list_nc, chunks=chunks) #TODO: does chunks argument solve "PerformanceWarning: Slicing is producing a large chunk."? {'time':1} is not a convenient chunking to use for timeseries extraction
    except xr.MergeError as e: #TODO: this except is necessary for CMCC, ux and uy have different lat/lon values, so renaming those of uy to avoid merging conflict
        def preprocess_CMCC_uovo(ds):
            if 'vo_' in os.path.basename(ds.encoding['source']):
                ds.coords['longitude'] = (ds.coords['longitude'] + 180) % 360 - 180 #normally this is done at convert_360to180, but inconvenient after renaming longitude variable
                ds = ds.rename({'longitude':'longitude_uy','latitude':'latitude_uy'})
                ds = ds.drop_vars(['vertices_longitude','vertices_latitude'])
            return ds
        print(f'catching "MergeError: {e}" >> WARNING: ux/uy have different latitude/longitude values, making two coordinates sets in Dataset.')
        data_xr = xr.open_mfdataset(file_list_nc, chunks=chunks, preprocess=preprocess_CMCC_uovo)
    
    #TODO: remove this commented code
    #rename variables with rename_dict derived from conversion_dict. duplicate keys are not possible, so phyc is always renamed to tracerbndOpal (last in conversion_dict)
    # rename_dict = {v['ncvarname']:k for k,v in conversion_dict.items()}
    # for ncvarn in data_xr.variables.mapping.keys():
    #     if ncvarn in rename_dict.keys():
    #         data_xr = data_xr.rename({ncvarn:rename_dict[ncvarn]})

    #renames ncvarnames to quantity names: proposal lisa, does not support ux/uy
    # for ncvarn in data_xr.variables.mapping.keys():
    #     if ncvarn == conversion_dict[quantity]['ncvarname']:
    #         data_xr = data_xr.rename({ncvarn:quantity})
    #         print(f'variable {ncvarn} renamed to {quantity}')
    
    for k,v in conversion_dict.items():
        ncvarn = v['ncvarname']
        if ncvarn in data_xr.variables.mapping.keys() and k in quantity_list: #k in quantity_list so phyc is not always renamed to tracerbndPON1 (first in conversion_dict)
            data_xr = data_xr.rename({ncvarn:k})
            print(f'variable {ncvarn} renamed to {k}')
    
    #rename dims time/depth/lat/lon/x/y #TODO: this has to be phased out some time, or made as an argument or merged with conversion_dict?
    rename_dims_dict = {'time_counter':'time', #time_counter instead of time for some CMCC files
                        'lev':'depth', #depth for CMEMS and many others, but lev for GFDL
                        'deptht':'depth', #deptht for some CMCC vars
                        'lon':'longitude','lat':'latitude',
                        'nav_lon':'longitude','nav_lat':'latitude', #nav_lon/nav_lat for some CMCC vars
                        'x':'j','y':'i', #x/y instead of j/i for some CMCC vars (non-regulargrid)
                        }
    for k,v in rename_dims_dict.items():
        if k in data_xr.dims and v not in data_xr.dims: 
            data_xr = data_xr.rename({k:v}) #TODO: can also do this for data_xr_var only?
            print(f'dimension {k} renamed to {v}')
    
    #get calendar and maybe convert_calendar, makes sure that nc_tstart/nc_tstop are of type pd._libs.tslibs.timestamps.Timestamp
    data_xr_calendar = data_xr['time'].dt.calendar
    if data_xr_calendar != 'proleptic_gregorian': #this is for instance the case in case of noleap (or 365_days) calendars from GFDL and CMCC
        units_copy = data_xr['time'].encoding['units']
        print(f'WARNING: calendar different than proleptic_gregorian found ({data_xr_calendar}), convert_calendar is called so check output carefully. It should be no issue for datasets with a monthly interval.')
        data_xr = data_xr.convert_calendar('standard') #TODO: does this not result in 29feb nan values in e.g. GFDL model? Check missing argument at https://docs.xarray.dev/en/stable/generated/xarray.Dataset.convert_calendar.html
        data_xr['time'].encoding['units'] = units_copy #put back dropped units
    
    #get timevar and compare requested dates
    timevar = data_xr['time']
    xr_tstartstop = pd.to_datetime(timevar.isel(time=[0,-1]).to_series())
    nc_tstart = xr_tstartstop.index[0]
    nc_tstop = xr_tstartstop.index[-1]
    if tstart < nc_tstart:
        raise OutOfRangeError(f'requested tstart {tstart} outside of available range {nc_tstart} to {nc_tstop}')
    if tstop > nc_tstop:
        raise OutOfRangeError(f'requested tstop {tstop} outside of available range {nc_tstart} to {nc_tstop}')
    
    #360 to 180 conversion
    convert_360to180 = (data_xr['longitude'].to_numpy()>180).any() #TODO: replace to_numpy() with load()
    latlon_ndims = len(data_xr['longitude'].shape)
    if convert_360to180: #TODO: make more flexible for models that eg pass -180/+180 crossing (add overlap at lon edges).
        data_xr.coords['longitude'] = (data_xr.coords['longitude'] + 180) % 360 - 180
        if latlon_ndims==1: #lon/lat has 1 dimension, .sortby() not possible if there are 2 dimensions
            data_xr = data_xr.sortby(data_xr['longitude'])
        else: #lon/lat is 2D #TODO: this can be removed
            print('WARNING: 2D latitude/longitude has more than one dim, continue without .sortby(). This is expected for e.g. CMCC')
    
    #retrieve var(s) (after potential longitude conversion) (also selecting relevant times)
    data_vars = list(data_xr.data_vars)
    bool_quanavailable = pd.Series(quantity_list).isin(data_vars)
    if not bool_quanavailable.all():
        quantity_list_notavailable = pd.Series(quantity_list).loc[~bool_quanavailable].tolist()
        raise KeyError(f'quantity {quantity_list_notavailable} not found in netcdf, available are: {data_vars}. Try updating conversion_dict to rename these variables.')
    data_xr_vars = data_xr[quantity_list].sel(time=slice(tstart,tstop))
    
    #optional conversion of units. Multiplications or other simple operatiors do not affect performance (dask.array(getitem) becomes dask.array(mul). With more complex operation it is better do do it on the interpolated array.
    for quan in quantity_list: #TODO: maybe do unit conversion before interp or is that computationally heavy?
        if 'conversion' in conversion_dict[quan].keys(): #if conversion is present, unit key must also be in conversion_dict
            print(f'> converting units from [{data_xr_vars[quan].attrs["units"]}] to [{conversion_dict[quan]["unit"]}]')
            #print(f'attrs are discarded:\n{data_xr_vars[quan].attrs}')
            data_xr_vars[quan] = data_xr_vars[quan] * conversion_dict[quan]['conversion'] #conversion drops all attributes of which units (which are changed anyway)
            data_xr_vars[quan].attrs['units'] = conversion_dict[quan]['unit'] #add unit attribute with resulting unit
    
    #optional refdate changing
    if refdate_str is not None:
        if 'long_name' in data_xr_vars.time.attrs: #for CMEMS it is 'hours since 1950-01-01', which would be wrong now #TODO: consider also removing attrs for depth/varname, since we would like to have salinitybnd/waterlevel instead of Salinity/sea_surface_height in xr plots?
            del data_xr_vars.time.attrs['long_name']
        data_xr_vars.time.encoding['units'] = refdate_str
    
    if 'depth' in data_xr_vars.coords:
        #make negative down
        if 'positive' in data_xr_vars['depth'].attrs.keys():
            if data_xr_vars['depth'].attrs['positive'] == 'down': #attribute appears in CMEMS, GFDL and CMCC, save to assume presence?
                data_xr_vars['depth'] = -data_xr_vars['depth']
        #optional reversing depth dimension for comparison to coastserv
        if reverse_depth:
            print('> reversing depth dimension')
            data_xr_vars = data_xr_vars.reindex({'depth':list(reversed(data_xr_vars['depth']))})
    
    return data_xr_vars


def interp_regularnc_to_plipoints(data_xr_reg, file_pli, nPoints=None, load=True):
    """
    load: interpolation errors are only raised upon loading, so do this per default
    """
    #TODO: make format of this dataset more in line with existing bnd-nc format and hisfile: https://issuetracker.deltares.nl/browse/UNST-6549
    data_xr_var = data_xr_reg #TODO: rename in script
    
    #load boundary file
    polyfile_object = hcdfm.PolyFile(file_pli)
    
    #check if polyobj names in plifile are unique
    polynames_pd = pd.Series([polyobj.metadata.name for polyobj in polyfile_object.objects])
    if polynames_pd.duplicated().any():
        raise ValueError(f'Duplicate polyobject names in polyfile {file_pli.name}, this is not allowed:\n{polynames_pd}')
    
    #create df of x/y/name of all plipoints in plifile
    data_pol_list = []
    for polyobj in polyfile_object.objects:
        data_pol_pd_one = pointlike_to_DataFrame(polyobj)
        data_pol_pd_one = data_pol_pd_one.iloc[:nPoints] #only use testset of points per polyobj in polyfile
        data_pol_pd_one['name'] = pd.Series(data_pol_pd_one.index).apply(lambda x: f'{polyobj.metadata.name}_{x+1:04d}')
        data_pol_list.append(data_pol_pd_one)
    data_pol_pd = pd.concat(data_pol_list)
    
    da_plipoints = xr.Dataset()
    da_plipoints['plipoint_x'] = xr.DataArray(data_pol_pd['x'], dims='plipoints')
    da_plipoints['plipoint_y'] = xr.DataArray(data_pol_pd['y'], dims='plipoints')
    da_plipoints['plipoint_name'] = xr.DataArray(data_pol_pd['name'].astype('S64'), dims='plipoints').str.decode('utf-8',errors='ignore').str.strip() #TODO: must be possible to do this less complex
    da_plipoints = da_plipoints.set_coords(['plipoint_x','plipoint_y','plipoint_name'])
    da_plipoints = da_plipoints.set_index({'plipoints':'plipoint_name'})
    
    #interpolation to lat/lon combinations
    print('> interp mfdataset to all PolyFile points (lat/lon coordinates)')
    #dtstart = dt.datetime.now()
    try:
        data_interp = data_xr_var.interp(longitude=da_plipoints['plipoint_x'], latitude=da_plipoints['plipoint_y'],
                                         method='linear', 
                                         kwargs={'bounds_error':True}, #error is only raised upon load(), so when the actual value retrieval happens
                                         )
    
    except ValueError as e: #Dimensions {'latitude', 'longitude'} do not exist. Expected one or more of Frozen({'time': 17, 'depth': 50, 'i': 292, 'j': 362}).
        #this is for eg CMCC model with multidimensional lat/lon variable
        #TODO: make nicer, without try except? eg latlon_ndims==1, but not sure if that is always valid >> add nonregular alternative for interp_regularnc_to_plipoints() and set kdtree to 1 (closest value) (uy stuff has to be dropped anyway)
        #TODO: maybe also spherical coordinate distance calculation instead of cartesian/eucledian
        #TODO: maybe use .sel(method='nearest'), but "KeyError: "no index found for coordinate 'longitude'""
        #TODO: interp for 2D also requested: https://github.com/pydata/xarray/issues/2281
        print(f'ValueError: {e}. Reverting to KDTree instead (nearest neigbour)')
        data_interp = xr.Dataset()
        for varone in list(data_xr_var.data_vars):
            path_lonlat_pd = data_pol_pd[['x','y']]
            if (varone=='uy') & (len(data_xr_var.data_vars)>1):
                data_lon_flat = data_xr_var['longitude_uy'].to_numpy().ravel()
                data_lat_flat = data_xr_var['latitude_uy'].to_numpy().ravel()
            else:
                data_lon_flat = data_xr_var['longitude'].to_numpy().ravel()
                data_lat_flat = data_xr_var['latitude'].to_numpy().ravel()
            data_lonlat_pd = pd.DataFrame({'x':data_lon_flat,'y':data_lat_flat})
            #KDTree, finds minimal eucledian distance between points (maybe haversine would be better)
            tree = KDTree(data_lonlat_pd) #alternatively sklearn.neighbors.BallTree: tree = BallTree(data_lonlat_pd)
            kdtree_k = 3 #TODO: nearest is probably just as wrong/right as weighted average, but raises error when using 1
            distance, data_lonlat_idx = tree.query(path_lonlat_pd, k=kdtree_k) #TODO: maybe add outofbounds treshold for distance
            #data_lonlat_pd.iloc[data_lonlat_idx]
            idx_i,idx_j = np.divmod(data_lonlat_idx, data_xr_var['longitude'].shape[1]) #get idx i and j by sort of counting over 2D array
            # import matplotlib.pyplot as plt
            # fig,ax = plt.subplots()
            # data_xr_var[varone].isel(time=0,depth=0).plot(ax=ax)
            # ax.plot(idx_j,idx_i,'xr')
            da_plipoints['da_idxi'] = xr.DataArray(idx_i, dims=('plipoints','nearestkpoints'))
            da_plipoints['da_idxj'] = xr.DataArray(idx_j, dims=('plipoints','nearestkpoints'))
            da_dist = xr.DataArray(distance, dims=('plipoints','nearestkpoints'))
            da_invdistweight = (1/da_dist)/(1/da_dist).sum(dim='nearestkpoints')
            da_varone_3k = data_xr_var[varone].isel(i=da_plipoints['da_idxi'],j=da_plipoints['da_idxj'])
            data_interp[varone] = (da_varone_3k * da_invdistweight).sum(dim='nearestkpoints')
            data_interp[varone].attrs = data_xr_var[varone].attrs #copy units and other attributes
    
    #time_passed = (dt.datetime.now()-dtstart).total_seconds()
    # print(f'>>time passed: {time_passed:.2f} sec')
    
    if not load:
        return data_interp
    
    print(f'> actual extraction of data from netcdf with .load() (for {len(data_pol_pd)} plipoints at once, this might take a while)')
    dtstart = dt.datetime.now()
    try:
        data_interp_loaded = data_interp.load() #loading data for all points at once is more efficient compared to loading data per point in loop 
    except ValueError as e: #generate a proper error with outofbounds requested coordinates, default is "ValueError: One of the requested xi is out of bounds in dimension 0" #TODO: improve error in xarray
        lonvar_vals = data_xr_var['longitude'].to_numpy()
        latvar_vals = data_xr_var['latitude'].to_numpy()
        data_pol_pd = data_interp[['plipoint_x','plipoint_y']].to_dataframe()
        bool_reqlon_outbounds = (data_pol_pd['plipoint_x'] <= lonvar_vals.min()) | (data_pol_pd['plipoint_x'] >= lonvar_vals.max())
        bool_reqlat_outbounds = (data_pol_pd['plipoint_y'] <= latvar_vals.min()) | (data_pol_pd['plipoint_y'] >= latvar_vals.max())
        reqlatlon_pd = pd.DataFrame({'longitude':data_pol_pd['plipoint_x'],'latitude':data_pol_pd['plipoint_y'],'lon outbounds':bool_reqlon_outbounds,'lat outbounds':bool_reqlat_outbounds})
        reqlatlon_pd_outbounds = reqlatlon_pd.loc[bool_reqlon_outbounds | bool_reqlat_outbounds]
        raise ValueError(f'{len(reqlatlon_pd_outbounds)} of requested pli points are out of bounds (valid longitude range {lonvar_vals.min()} to {lonvar_vals.max()}, valid latitude range {latvar_vals.min()} to {latvar_vals.max()}):\n{reqlatlon_pd_outbounds}')
    time_passed = (dt.datetime.now()-dtstart).total_seconds()
    print(f'>>time passed: {time_passed:.2f} sec')

    return data_interp_loaded


def interp_uds_to_plipoints(uds:xu.UgridDataset, gdf:geopandas.GeoDataFrame, nPoints:int=None) -> xr.Dataset:
    """
    To interpolate an unstructured dataset (like a *_map.nc file) read with xugrid to plipoint locations
    
    Parameters
    ----------
    uds : xu.UgridDataset
        dfm model output read using dfm_tools. Dims: mesh2d_nLayers, mesh2d_nInterfaces, time, mesh2d_nNodes, mesh2d_nFaces, mesh2d_nMax_face_nodes, mesh2d_nEdges.
    gdf : geopandas.GeoDataFrame (str/path is also supported)
        gdf with location, geometry (Point) and crs corresponding to model crs.
    nPoints : int, optional
        amount of points (None gives all). The default is None.

    Raises
    ------
    Exception
        DESCRIPTION.

    Returns
    -------
    ds : TYPE
        xr.Dataset with dims: plipoints, time, depth.

    """
    facedim = uds.grid.face_dimension
    
    if isinstance(gdf,(str,Path)): #TODO: align plipoints/gdf with other functions, now three input types are supported, but the interp_regularnc_to_plipoints requires paths to plifiles (and others also)
        gdf = PolyFile_to_geodataframe_points(hcdfm.PolyFile(gdf))
    
    gdf = gdf.iloc[:nPoints]
    ds = uds.ugrid.sel_points(x=gdf.geometry.x, y=gdf.geometry.y)
    #TODO: drop mesh2d_face_x and mesh2d_face_y variables
    
    if len(gdf)!=ds.dims[facedim]: #TODO: check this until https://github.com/Deltares/xugrid/issues/100 is solved, after that, make a testcase that checks only this if-statement
        ds_points = geopandas.points_from_xy(ds.x,ds.y)
        gdfpoint_inds_bool = pd.Series(index=range(len(gdf)))
        gdfpoint_inds_bool[:] = True
        for iR, gdf_row in gdf.iterrows():
            gdf_point = gdf_row.geometry
            if gdf_point in ds_points:
                gdfpoint_inds_bool.iloc[iR] = False
        gdf_stats = gdf.copy()
        gdf_stats['missing'] = gdfpoint_inds_bool
        raise ValueError(f'requested {len(gdf)} points but resulted in ds with {ds.dims[facedim]} points, missing points are probably outside of the uds model domain:\n{gdf_stats}')

    ds = ds.rename({facedim:'plipoints'}) # rename mesh2d_nFaces to plipoints
    
    ds['plipoint_name'] = xr.DataArray(gdf['plipoint_name'].tolist(), dims='plipoints') # change name of plipoint (node to gdf name)
    ds = ds.set_index({'plipoints':'plipoint_name'})
    return ds


def interp_hisnc_to_plipoints(data_xr_his, file_pli, kdtree_k=3, load=True):
    """
    interpolate stations in hisfile to points of polyfile
    """
    #KDTree, finds minimal eucledian distance between points (haversine would be better). Alternatively sklearn.neighbors.BallTree: tree = BallTree(data_lonlat_pd)
    
    datavars_list = list(data_xr_his.data_vars)
    if len(datavars_list)>5:
        print('WARNING: more than 5 data_vars, you might want to subset your data_xr_his')
    #read hisfile and make KDTree of xy of hisstations
    hisstations_pd = data_xr_his.stations.to_dataframe() #TODO: add check if stations are strings? (whether preprocess_hisnc was used)
    tree_nest2 = KDTree(hisstations_pd[['station_x_coordinate','station_y_coordinate']])
    
    #read polyfile and query k nearest hisstations (names)
    polyfile_object = hcdfm.PolyFile(file_pli)
    data_pol_list = []
    for polyobj in polyfile_object.objects:
        data_pol_pd_one = pointlike_to_DataFrame(polyobj)
        data_pol_pd_one['name'] = pd.Series(data_pol_pd_one.index).apply(lambda x: f'{polyobj.metadata.name}_{x+1:04d}')
        data_pol_list.append(data_pol_pd_one)
    data_pol_pd = pd.concat(data_pol_list)

    plicoords_distance2, plicoords_nestpointidx = tree_nest2.query(data_pol_pd[['x','y']], k=kdtree_k)
    da_plicoords_nestpointidx = xr.DataArray(plicoords_nestpointidx, dims=('plipoints','nearestkpoints'))
    da_plicoords_nestpointnames = data_xr_his.stations.isel(stations=da_plicoords_nestpointidx)
    
    #interpolate hisfile variables to plipoints
    data_interp = xr.Dataset()
    data_interp['plipoint_x'] = xr.DataArray(data_pol_pd['x'],dims=('plipoints'))
    data_interp['plipoint_y'] = xr.DataArray(data_pol_pd['y'],dims=('plipoints'))
    data_interp['plipoint_name'] = xr.DataArray(data_pol_pd['name'].astype('S64'), dims='plipoints').str.decode('utf-8',errors='ignore').str.strip() #TODO: must be possible to do this less complex
    data_interp = data_interp.set_coords(['plipoint_x','plipoint_y','plipoint_name'])
    data_interp = data_interp.set_index({'plipoints':'plipoint_name'})
    
    for varone in datavars_list:
        da_dist = xr.DataArray(plicoords_distance2, dims=('plipoints','nearestkpoints'))
        da_invdistweight = (1/da_dist)/(1/da_dist).sum(dim='nearestkpoints') #TODO: set weigths for invalid points to 0 (and increase others weights so sum remains 1)
        data_xr_his_pli_knearest = data_xr_his[varone].sel(stations=da_plicoords_nestpointnames)
        data_interp[varone] = (data_xr_his_pli_knearest * da_invdistweight).sum(dim='nearestkpoints')
        data_interp[varone].attrs = data_xr_his[varone].attrs #copy units and other attributes
        
    if not load:
        return data_interp
    
    print('loading dataset (might take a while, but increases performance of point loop significantly)')
    data_interp = data_interp.load()
    print('done')

    return data_interp
    

def plipointsDataset_to_ForcingModel(plipointsDataset):
    """
    empty docstring
    """
    quantity_list = list(plipointsDataset.data_vars)
    npoints = len(plipointsDataset.plipoints)
    
    #start conversion to Forcingmodel object
    print(f'Converting {npoints} plipoints to hcdfm.ForcingModel():',end='')
    dtstart = dt.datetime.now()
    ForcingModel_object = hcdfm.ForcingModel()
    for iP in range(npoints):
        print(f' {iP+1}',end='')
        
        #select data for this point, ffill nans, concatenating time column, constructing T3D/TimeSeries and append to hcdfm.ForcingModel()
        datablock_xr_onepoint = plipointsDataset.isel(plipoints=iP)
        plipoint_name = str(datablock_xr_onepoint.plipoints.to_numpy())
        
        for quan in quantity_list:
            datablock_xr_onepoint[quan].attrs['locationname'] = plipoint_name #TODO: is there a nicer way of passing this data?
            if datablock_xr_onepoint[quan].isnull().all(): # check if all values of plipoint are nan (on land)
                warnings.warn(UserWarning(f'Plipoint "{plipoint_name}" might be on land since it only contain nan values. Consider altering your plifile or using plipointsDataset.ffill(dim="plipoints").bfill(dim="plipoints"). Nan values replaced with 0 to avoid bc-writing errors')) #TODO: maybe fill along plipoints dimension, but beware on nans in deep water that are filled with neighbours
                datablock_xr_onepoint[quan] = datablock_xr_onepoint[quan].fillna(0)
        
        if 'depth' in plipointsDataset.dims:
            ts_one = Dataset_to_T3D(datablock_xr_onepoint)
        elif 'amplitude' in quantity_list:
            ts_one = Dataset_to_Astronomic(datablock_xr_onepoint)
        else:
            ts_one = Dataset_to_TimeSeries(datablock_xr_onepoint)
        ForcingModel_object.forcing.append(ts_one)
    print(f'. >> done in {(dt.datetime.now()-dtstart).total_seconds():.2f} sec')
    
    return ForcingModel_object

