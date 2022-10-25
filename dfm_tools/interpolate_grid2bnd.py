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
from pathlib import Path

from hydrolib.core.io.bc.models import (
    ForcingModel,
    QuantityUnitPair,
    Astronomic,
)
from hydrolib.core.io.polyfile.models import PolyFile

from dfm_tools.hydrolib_helpers import Dataset_to_TimeSeries, Dataset_to_T3D


def get_conversion_dict(ncvarname_updates={}):
    
    """
    get the conversion_dict, optionally with updated ncvarnames
    conversion_dict.keys() are the dflowfm quantities and the ncvarname the corresponding netcdf variable name/key
    alternative ncvarnames can be supplied via ncvarname_updates, e.g. get_conversion_dict(ncvarname_updates={'temperaturebnd':'tos'})
    
    interpolate_nc_to_bc() renames netcdf variable like this:
    data_xr = data_xr.rename({ncvarname:quantity})
    """
    # conversion_dict, 
    conversion_dict = { # mg/l is the same as g/m3: conversion is phyc in mmol/l to newvar in g/m3
                        'tracerbndOXY'        : {'ncvarname': 'o2',          'unit': 'g/m3', 'conversion': 32.0 / 1000.0}, 
                        'tracerbndNO3'        : {'ncvarname': 'no3',         'unit': 'g/m3', 'conversion': 14.0 / 1000.0},
                        'tracerbndPO4'        : {'ncvarname': 'po4',         'unit': 'g/m3', 'conversion': 30.97 / 1000.0},
                        'tracerbndSi'         : {'ncvarname': 'si',          'unit': 'g/m3', 'conversion': 28.08 / 1000.0},
                        'tracerbndPON1'       : {'ncvarname': 'phyc',        'unit': 'g/m3', 'conversion': 2. * 16. * 14. / (106. * 1000.0)},
                        'tracerbndPOP1'       : {'ncvarname': 'phyc',        'unit': 'g/m3', 'conversion': 2. * 30.97 / (106. * 1000.0)},
                        'tracerbndPOC1'       : {'ncvarname': 'phyc',        'unit': 'g/m3', 'conversion': 2. * 12. / 1000.0},
                        'tracerbndDON'        : {'ncvarname': 'phyc',        'unit': 'g/m3', 'conversion': 3.24 * 2. * 16. * 14. / (106. * 1000.0)},
                        'tracerbndDOP'        : {'ncvarname': 'phyc',        'unit': 'g/m3', 'conversion': 1.0 * 2. * 30.97 / (106. * 1000.0)},
                        'tracerbndDOC'        : {'ncvarname': 'phyc',        'unit': 'g/m3', 'conversion': (199. / 20.) * 3.24 * 2. * 16. * 12. / (106. * 1000.0)},
                        'tracerbndOpal'       : {'ncvarname': 'phyc',        'unit': 'g/m3', 'conversion': 0.5 * 0.13 * 28.08 / (1000.0)},
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
    """
    # translate dict from .\hydro_tools\FES\PreProcessing_FES_TideModel_imaginary.m
    #component_list = ['2N2','LABDA2','MF','MFM','P1','SSA','EPSILON2','M2','MKS2','MU2','Q1','T2','J1','M3','MM','N2','R2','K1','M4','MN4','N4','S1','K2','M6','MS4','NU2','S2','L2','M8','MSF','O1','S4','MSQM','SA']
    translate_dict = {'LA2':'LABDA2', #TODO: use value instead of key in bc file? Support using value instead of key in const_list also (like line above)
                      'EPS2':'EPSILON2', 
                      'Z0':'A0',
                      'MTM':'MFM', #Needs to be verified
                      }
    
    print('initialize ForcingModel()')
    ForcingModel_object = ForcingModel()
    
    dir_pattern_dict = {'FES2014': Path(r'P:\metocean-data\licensed\FES2014','*.nc'), #ocean_tide_extrapolated
                        'FES2012': Path(r'P:\metocean-data\open\FES2012\data','*_FES2012_SLEV.nc'), #is eigenlijk ook licensed
                        'EOT20': Path(r'P:\metocean-data\open\EOT20\ocean_tides','*_ocean_eot20.nc'),
                        }
    if tidemodel not in dir_pattern_dict.keys():
        raise Exception(f'invalid tidemodel "{tidemodel}", options are: {list(dir_pattern_dict.keys())}')
    dir_pattern = dir_pattern_dict[tidemodel]
    
    if component_list is None:
        file_list_nc = glob.glob(str(dir_pattern))
        component_list = [os.path.basename(x).replace('.nc','').replace('_FES2012_SLEV','').replace('_ocean_eot20','') for x in file_list_nc] #TODO: make this less hard-coded
    else:
        file_list_nc = [str(dir_pattern).replace('*',comp) for comp in component_list]
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
    
    #derive uv phase components (using amplitude=1)
    data_xrsel_phs_rad = np.deg2rad(data_xrsel['phase'])
    #we need to compute u/v components for the phase to avoid zero-crossing interpolation issues
    data_xrsel['phase_u'] = 1*np.cos(data_xrsel_phs_rad)
    data_xrsel['phase_v'] = 1*np.sin(data_xrsel_phs_rad)
    
    #load boundary file
    polyfile_object = PolyFile(file_pli)
    pli_PolyObjects = polyfile_object.objects

    for iPO, pli_PolyObject_sel in enumerate(pli_PolyObjects):
        print(f'processing PolyObject {iPO+1} of {len(pli_PolyObjects)}: name={pli_PolyObject_sel.metadata.name}')
        
        #create requestedlat/requestedlon DataArrays for proper interpolation in xarray (with new dimension name)
        path_lons = np.array([point.x for point in pli_PolyObject_sel.points])[:nPoints]
        path_lats = np.array([point.y for point in pli_PolyObject_sel.points])[:nPoints]
        da_lons = xr.DataArray(path_lons, dims='latloncombi')
        da_lats = xr.DataArray(path_lats, dims='latloncombi')
        
        #interpolation to lat/lon combinations
        print('> interp mfdataset with all lat/lon coordinates and compute phase_new')
        dtstart = dt.datetime.now()
        data_interp_allcoords = data_xrsel.interp(lat=da_lats,lon=da_lons,
                                                  method='linear', 
                                                  kwargs={'bounds_error':True})
        data_interp_allcoords['phase_new'] = np.rad2deg(np.arctan2(data_interp_allcoords['phase_v'],data_interp_allcoords['phase_u']))
        time_passed = (dt.datetime.now()-dtstart).total_seconds()
        # print(f'>>time passed: {time_passed:.2f} sec')
        
        
        print('> actual extraction of data from netcdf with .load() (for all PolyObject points at once, so this will take a while)')
        dtstart = dt.datetime.now()
        try:
            data_interp_amp_allcoords = data_interp_allcoords['amplitude'].load() #convert from cm to m
            data_interp_phs_allcoords = data_interp_allcoords['phase_new'].load()
        except ValueError: #generate a proper error with outofbounds requested coordinates, default is "ValueError: One of the requested xi is out of bounds in dimension 0"
            lonvar_vals = data_xrsel['lon'].to_numpy()
            latvar_vals = data_xrsel['lat'].to_numpy()
            bool_reqlon_outbounds = (path_lons <= lonvar_vals.min()) | (path_lons >= lonvar_vals.max())
            bool_reqlat_outbounds = (path_lats <= latvar_vals.min()) | (path_lats >= latvar_vals.max())
            reqlatlon_pd = pd.DataFrame({'longitude':path_lons,'latitude':path_lats,'lon outbounds':bool_reqlon_outbounds,'lat outbounds':bool_reqlat_outbounds})
            reqlatlon_pd_outbounds = reqlatlon_pd.loc[bool_reqlon_outbounds | bool_reqlat_outbounds]
            raise ValueError(f'{len(reqlatlon_pd_outbounds)} of {len(reqlatlon_pd)} requested pli points are out of bounds (valid longitude range {lonvar_vals.min()} to {lonvar_vals.max()}, valid latitude range {latvar_vals.min()} to {latvar_vals.max()}):\n{reqlatlon_pd_outbounds}')
        time_passed = (dt.datetime.now()-dtstart).total_seconds()
        print(f'>>time passed: {time_passed:.2f} sec')
        
        #convert cm to m
        if data_xrsel['amplitude'].attrs['units'] == 'cm':
            data_interp_amp_allcoords = data_interp_amp_allcoords/100
            data_xrsel['amplitude'].attrs['units'] = 'm'
        
        
        for iP, pli_Point_sel in enumerate(pli_PolyObject_sel.points[:nPoints]): #looping over plipoints within component loop, append to datablock_pd_allcomp
            print(f'processing Point {iP+1} of {len(pli_PolyObject_sel.points)}: ',end='')
            lonx, laty = pli_Point_sel.x, pli_Point_sel.y
            print(f'(x={lonx}, y={laty})')
            pli_PolyObject_name_num = f'{pli_PolyObject_sel.metadata.name}_{iP+1:04d}'
            
            data_interp_amp = data_interp_amp_allcoords.sel(latloncombi=iP)
            data_interp_phs = data_interp_phs_allcoords.sel(latloncombi=iP)
            
            datablock_pd_allcomp = pd.DataFrame({'component':component_list_upper_pd,'amp':data_interp_amp.to_numpy(),'phs':data_interp_phs.to_numpy()})
            if datablock_pd_allcomp['amp'].isnull().any():
                print('WARNING: only nans for this coordinate, this point might be on land')
            datablock_list = datablock_pd_allcomp.values.tolist()
            
            # Each .bc file can contain 1 or more timeseries, one for each support point:
            print('> constructing TimeSeries and appending to ForcingModel()')
            ts_one = Astronomic(name=pli_PolyObject_name_num,
                                quantityunitpair=[QuantityUnitPair(quantity="astronomic component", unit='-'),
                                                  QuantityUnitPair(quantity='waterlevelbnd amplitude', unit=data_xrsel['amplitude'].attrs['units']),
                                                  QuantityUnitPair(quantity='waterlevelbnd phase', unit=data_xrsel['phase'].attrs['units'])],
                                datablock=datablock_list, 
                                )
            
            ForcingModel_object.forcing.append(ts_one)        
    return ForcingModel_object


def interpolate_nc_to_bc(dir_pattern, file_pli, quantity, 
                         tstart, tstop, refdate_str=None, 
                         conversion_dict=None,
                         nPoints=None, #argument for testing
                         reverse_depth=False, #temporary argument to compare easier with old coastserv files
                         KDTree_invdistweigth=False):
    
    if conversion_dict is None:
        conversion_dict = get_conversion_dict()
    
    if quantity=='uxuy': #T3Dvector
        quantity_list = ['ux','uy']
        ncvarname_list = [conversion_dict[quan]['ncvarname'] for quan in quantity_list]
    else:
        quantity_list = [quantity]
        ncvarname_list = [conversion_dict[quan]['ncvarname'] for quan in quantity_list]
    
    print('initialize ForcingModel()')
    ForcingModel_object = ForcingModel()
    
    dir_pattern = [Path(str(dir_pattern).format(ncvarname=ncvarname)) for ncvarname in ncvarname_list]
    file_list_nc = []
    for dir_pattern_one in dir_pattern:
        file_list_nc = file_list_nc + glob.glob(str(dir_pattern_one))
    list_pattern_names = [x.name for x in dir_pattern]
    print(f'loading mfdataset of {len(file_list_nc)} files with pattern(s) {list_pattern_names}')
    
    dtstart = dt.datetime.now()
    try:
        data_xr = xr.open_mfdataset(file_list_nc,chunks={'time':1}) #TODO: does chunks argument solve "PerformanceWarning: Slicing is producing a large chunk."?
    except xr.MergeError as e: #TODO: this except is necessary for CMCC, ux and uy have different lat/lon values, so renaming those of uy to avoid merging conflict
        def preprocess_CMCC_uovo(ds):
            if 'vo_' in os.path.basename(ds.encoding['source']):
                ds.coords['longitude'] = (ds.coords['longitude'] + 180) % 360 - 180 #normally this is done at convert_360to180, but inconvenient after renaming longitude variable
                ds = ds.rename({'longitude':'longitude_uy','latitude':'latitude_uy'})
                ds = ds.drop_vars(['vertices_longitude','vertices_latitude'])
            return ds
        print(f'catching "MergeError: {e}" >> WARNING: ux/uy have different latitude/longitude values, making two coordinates sets in Dataset.')
        data_xr = xr.open_mfdataset(file_list_nc,chunks={'time':1},preprocess=preprocess_CMCC_uovo)
        
    #rename variables with rename_dict derived from conversion_dict
    rename_dict = {v['ncvarname']:k for k,v in conversion_dict.items()}
    for ncvarn in data_xr.variables.mapping.keys():
        if ncvarn in rename_dict.keys():
            data_xr = data_xr.rename({ncvarn:rename_dict[ncvarn]})
    
    #change refdate
    if refdate_str is not None:
        data_xr.time.encoding['units'] = refdate_str
    
    #get calendar and maybe convert_calendar, makes sure that nc_tstart/nc_tstop are of type pd._libs.tslibs.timestamps.Timestamp
    data_xr_calendar = data_xr['time'].dt.calendar
    if data_xr_calendar != 'proleptic_gregorian': #this is for instance the case in case of noleap (or 365_days) calendars from GFDL and CMCC
        units_copy = data_xr['time'].encoding['units']
        print('WARNING: calendar different than proleptic_gregorian found ({data_xr_calendar}), convert_calendar is called so check output carefully. It should be no issue for datasets with a monthly interval.')
        data_xr = data_xr.convert_calendar('standard') #TODO: does this not result in 29feb nan values in e.g. GFDL model? Check missing argument at https://docs.xarray.dev/en/stable/generated/xarray.Dataset.convert_calendar.html
        data_xr['time'].encoding['units'] = units_copy #put back dropped units
    time_passed = (dt.datetime.now()-dtstart).total_seconds()
    # print(f'>>time passed: {time_passed:.2f} sec')
    
    #get timevar and compare requested dates
    timevar = data_xr['time']
    xr_tstartstop = pd.to_datetime(timevar.isel(time=[0,-1]).to_series())
    nc_tstart = xr_tstartstop.index[0]
    nc_tstop = xr_tstartstop.index[-1]
    if tstart < nc_tstart:
        raise Exception(f'requested tstart {tstart} outside of available range {nc_tstart} to {nc_tstop}')
    if tstop > nc_tstop:
        raise Exception(f'requested tstop {tstop} outside of available range {nc_tstart} to {nc_tstop}')
    
    #rename coordinates
    if 'depth' not in data_xr.coords:
        if 'lev' in data_xr.coords: #depth for CMEMS and many others, but lev for GFDL, convert to depth #TODO: provide rename_dict as argument to this function or leave as is?
            data_xr = data_xr.rename({'lev':'depth'}) #TODO: can also do this for data_xr_var only?
            print('variable/coordinate lev renamed to depth')
    coordname_lon,coordname_lat = ['longitude','latitude'] #hardcode this
    if coordname_lon not in data_xr.coords:
        if 'lon' in data_xr.coords:
            data_xr = data_xr.rename({'lon':coordname_lon,'lat':coordname_lat})
        else:
            raise Exception(f'no lat/lon coords available in file: {data_xr.coords}')
    
    #360 to 180 conversion
    convert_360to180 = (data_xr[coordname_lon].to_numpy()>180).any()
    latlon_ndims = len(data_xr[coordname_lon].shape)
    if convert_360to180: #TODO: make more flexible for models that eg pass -180/+180 crossing (add overlap at lon edges).
        data_xr.coords[coordname_lon] = (data_xr.coords[coordname_lon] + 180) % 360 - 180
        if latlon_ndims==1: #lon/lat has 1 dimension, .sortby() not possible if there are 2 dimensions
            data_xr = data_xr.sortby(data_xr[coordname_lon])
        else: #lon/lat is 2D #TODO: this can be removed
            print('WARNING: 2D latitude/longitude has more than one dim, continue without .sortby(). This is expected for e.g. CMCC')
    
    #retrieve var(s) (after potential longitude conversion) (also selecting relevant times)
    data_vars = list(data_xr.data_vars)
    bool_quanavailable = pd.Series(quantity_list).isin(data_vars)
    if not bool_quanavailable.all():
        quantity_list_notavailable = pd.Series(quantity_list).loc[~bool_quanavailable].tolist()
        raise Exception(f'quantity {quantity_list_notavailable} not found, available are: {data_vars}. Try updating conversion_dict to rename these variables.')
    data_xr_var = data_xr[quantity_list].sel(time=slice(tstart,tstop))
    
    if coordname_lon not in data_xr_var.coords:
        raise Exception(f'{coordname_lon} not in variable coords: {data_xr_var.coords}.')
    
    #load boundary file
    polyfile_object = PolyFile(file_pli)
    pli_PolyObjects = polyfile_object.objects
    
    for iPO, pli_PolyObject_sel in enumerate(pli_PolyObjects):
        print(f'processing PolyObject {iPO+1} of {len(pli_PolyObjects)}: name={pli_PolyObject_sel.metadata.name}')
        
        #create requestedlat/requestedlon DataArrays for proper interpolation in xarray (with new dimension name)
        path_lons = np.array([point.x for point in pli_PolyObject_sel.points])[:nPoints]
        path_lats = np.array([point.y for point in pli_PolyObject_sel.points])[:nPoints]
        da_lons = xr.DataArray(path_lons, dims='latloncombi')
        da_lats = xr.DataArray(path_lats, dims='latloncombi')
        
        #interpolation to lat/lon combinations
        print('> interp mfdataset to all PolyObject points (lat/lon coordinates)')
        dtstart = dt.datetime.now()
        try:
            data_interp = data_xr_var.interp({coordname_lat:da_lats, coordname_lon:da_lons}, #also possible without dict: (latitude=da_lats, longitude=da_lons), but this is more flexible
                                             method='linear', 
                                             kwargs={'bounds_error':True}, #error is only raised upon load(), so when the actual value retrieval happens
                                             )
        except ValueError as e: #Dimensions {'latitude', 'longitude'} do not exist. Expected one or more of Frozen({'time': 17, 'depth': 50, 'i': 292, 'j': 362}).
            #this is for eg CMCC model with multidimensional lat/lon variable
            #TODO: make nicer, without try except? eg latlon_ndims==1, but not sure if that is always valid
            #TODO: maybe also spherical coordinate distance calculation instead of cartesian/eucledian
            print(f'ValueError: {e}. Reverting to KDTree instead (nearest neigbour)')
            from scipy.spatial import KDTree #TODO: move up
            data_interp = xr.Dataset()
            for varone in list(data_xr_var.data_vars):
                path_lonlat_pd = pd.DataFrame({'lon':da_lons,'lat':da_lats})
                if (varone=='uy') & (len(data_xr_var.data_vars)>1):
                    data_lon_flat = data_xr_var['longitude_uy'].to_numpy().ravel()
                    data_lat_flat = data_xr_var['latitude_uy'].to_numpy().ravel()
                else:
                    data_lon_flat = data_xr_var[coordname_lon].to_numpy().ravel()
                    data_lat_flat = data_xr_var[coordname_lat].to_numpy().ravel()
                data_lonlat_pd = pd.DataFrame({'lon':data_lon_flat,'lat':data_lat_flat})
                #KDTree, finds minimal eucledian distance between points (maybe haversine would be better)
                tree = KDTree(data_lonlat_pd) #alternatively sklearn.neighbors.BallTree: tree = BallTree(data_lonlat_pd)
                distance, data_lonlat_idx = tree.query(path_lonlat_pd, k=3) #TODO: maybe add outofbounds treshold for distance
                #data_lonlat_pd.iloc[data_lonlat_idx]
                idx_i,idx_j = np.divmod(data_lonlat_idx, data_xr_var['longitude'].shape[1]) #get idx i and j by sort of counting over 2D array
                # import matplotlib.pyplot as plt
                # fig,ax = plt.subplots()
                # data_xr_var[varone].isel(time=0,depth=0).plot(ax=ax)
                # ax.plot(idx_j,idx_i,'xr')
                da_idxi = xr.DataArray(idx_i, dims=('latloncombi','point3k'))
                da_idxj = xr.DataArray(idx_j, dims=('latloncombi','point3k'))
                if KDTree_invdistweigth:
                    da_dist = xr.DataArray(distance, dims=('latloncombi','point3k'))
                    da_invdistweight = (1/da_dist)/(1/da_dist).sum(dim='point3k')
                    da_varone_3k = data_xr_var[varone].isel(i=da_idxi,j=da_idxj)
                    data_interp[varone] = (da_varone_3k * da_invdistweight).sum(dim='point3k')
                    data_interp[varone].attrs = data_xr_var[varone].attrs #copy units and other attributes
                else:
                    data_interp[varone] = data_xr_var[varone].isel(i=da_idxi,j=da_idxj).isel(point3k=0) #take first of 3 closest points

        time_passed = (dt.datetime.now()-dtstart).total_seconds()
        # print(f'>>time passed: {time_passed:.2f} sec')
        
        print('> actual extraction of data from netcdf with .load() (for all PolyObject points at once, so this will take a while)')
        dtstart = dt.datetime.now()
        try:
            datablock_xr_allpoints = data_interp.load() #loading data for all points at once is more efficient compared to loading data per point in loop 
        except ValueError: #generate a proper error with outofbounds requested coordinates, default is "ValueError: One of the requested xi is out of bounds in dimension 0"
            lonvar_vals = data_xr[coordname_lon].to_numpy()
            latvar_vals = data_xr[coordname_lat].to_numpy()
            bool_reqlon_outbounds = (path_lons <= lonvar_vals.min()) | (path_lons >= lonvar_vals.max())
            bool_reqlat_outbounds = (path_lats <= latvar_vals.min()) | (path_lats >= latvar_vals.max())
            reqlatlon_pd = pd.DataFrame({'longitude':path_lons,'latitude':path_lats,'lon outbounds':bool_reqlon_outbounds,'lat outbounds':bool_reqlat_outbounds})
            reqlatlon_pd_outbounds = reqlatlon_pd.loc[bool_reqlon_outbounds | bool_reqlat_outbounds]
            raise ValueError(f'{len(reqlatlon_pd_outbounds)} of requested pli points are out of bounds (valid longitude range {lonvar_vals.min()} to {lonvar_vals.max()}, valid latitude range {latvar_vals.min()} to {latvar_vals.max()}):\n{reqlatlon_pd_outbounds}')
        time_passed = (dt.datetime.now()-dtstart).total_seconds()
        print(f'>>time passed: {time_passed:.2f} sec')
        
        #optional conversion of units and reversing depth dimension
        if 'conversion' in conversion_dict.keys(): #if conversion is present, unit key must also be in conversion_dict
            for quan in quantity_list:
                print(f'> converting units from [{datablock_xr_allpoints[quan].attrs["units"]}] to [{conversion_dict["unit"]}]')
                datablock_xr_allpoints[quan] = datablock_xr_allpoints[quan] * conversion_dict['conversion'] #conversion drops all attributes of which units (which are changed anyway)
                datablock_xr_allpoints[quan].attrs['units'] = conversion_dict['unit'] #add unit attribute with resulting unit
        
        if ('depth' in data_xr_var.coords) & reverse_depth:
            print('> reversing depth dimension')
            datablock_xr_allpoints = datablock_xr_allpoints.reindex({'depth':list(reversed(datablock_xr_allpoints['depth']))})
        
        for iP, pli_Point_sel in enumerate(pli_PolyObject_sel.points[:nPoints]):
            print(f'processing Point {iP+1} of {len(pli_PolyObject_sel.points)}: ',end='')
            lonx_print, laty_print = pli_Point_sel.x, pli_Point_sel.y
            print(f'(x={lonx_print}, y={laty_print})')
            pli_PolyObject_name_num = f'{pli_PolyObject_sel.metadata.name}_{iP+1:04d}'
            
            print('> select data for this point, ffill nans, concatenating time column, constructing T3D/TimeSeries and appending to ForcingModel()')
            dtstart = dt.datetime.now()
            datablock_xr_onepoint = datablock_xr_allpoints.isel(latloncombi=iP)
            
            for quan in quantity_list:
                datablock_xr_onepoint[quan].attrs['locationname'] = pli_PolyObject_name_num #TODO: is there a nicer way of passing this data?
                
                if np.isnan(datablock_xr_onepoint[quan].to_numpy()).all(): # check if only nan (out of bounds or land) # we can do .to_numpy() without performance loss, since data is already loaded in datablock_xr_allpoints
                    print('WARNING: only nans for this coordinate, this point might be on land')
            
            if 'depth' in data_xr_var.coords:
                ts_one = Dataset_to_T3D(datablock_xr_onepoint)
            else:
                ts_one = Dataset_to_TimeSeries(datablock_xr_onepoint)
            ForcingModel_object.forcing.append(ts_one)
            time_passed = (dt.datetime.now()-dtstart).total_seconds()
            # print(f'>>time passed: {time_passed:.2f} sec')
    
    return ForcingModel_object


