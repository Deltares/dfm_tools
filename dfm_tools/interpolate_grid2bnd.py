# -*- coding: utf-8 -*-
"""
Created on Thu Aug 18 17:39:03 2022

@author: veenstra

Some blocking hydrolib issues before this tool can replace coastserv:
-	uxuy is not yet programmed, add support for this merged array
    o	new bug issue: https://github.com/Deltares/HYDROLIB-core/issues/316
-	metadata header is not completely correct yet, but this seems like small fixes in hydrolib:
    o	new bug issue: https://github.com/Deltares/HYDROLIB-core/issues/317
-	Formatting of datablock writing in bc file can be improved (probably also improves perfomance):
    o	Existing issue: https://github.com/Deltares/HYDROLIB-core/issues/308
    o	Existing issue: https://github.com/Deltares/HYDROLIB-core/issues/313
-	some minor issues that do not seem blocking (my issues in the range of #305 to #322)

Non-hydrolib things to do (still missing compared to coastserv):
-	downloading part is not included (Bjorn improved that part so add to dfm_tools?) >> https://github.com/c-scale-community/use-case-hisea/blob/main/scripts/download/download_cmems_biogeochemistry.py and download_cmems_physics.py
-	FES is included but there are differences with the reference bc files for DCSM (uv or complex method is necessary for when phs crosses 0)
-	Waq variables incl unit conversion work for e.g. CMEMS and GFDL (and probably most other models), but CMCC has no lat/lon coords so currently crashes
-   other waq variables than NO3 (incl conversion) probably work, but were not checked yet

"""

import os
import glob
import datetime as dt
import numpy as np
import pandas as pd
import xarray as xr

from hydrolib.core.io.bc.models import (
    ForcingModel,
    QuantityUnitPair,
    Astronomic,
)
from hydrolib.core.io.polyfile.models import PolyFile

from dfm_tools.hydrolib_helpers import DataArray_to_TimeSeries, DataArray_to_T3D


def get_conversion_dict(model='CMEMS'):
    """
    interpolate_nc_to_bc() renames netcdf variable like this:
    data_xr = data_xr.rename({ncvarname:bcvarname})
    """
    conversion_dicts = {}
    # conversion_dict, contains ncvarname as array because uxuy relies on 2 CMEMS variables
    conversion_dicts['CMEMS'] = { # mg/l is the same as g/m3: conversion is phyc in mmol/l to newvar in g/m3
                                'OXY'        : {'ncvarname': 'o2',         'bcvarname': 'tracerbndOXY',  'unit': 'g/m3', 'conversion' : 32.0 / 1000.0}, 
                                'NO3'        : {'ncvarname': 'no3',        'bcvarname': 'tracerbndNO3',  'unit': 'g/m3', 'conversion' : 14.0 / 1000.0},
                                'PO4'        : {'ncvarname': 'po4',        'bcvarname': 'tracerbndPO4',  'unit': 'g/m3', 'conversion' : 30.97 / 1000.0},
                                'Si'         : {'ncvarname': 'si',         'bcvarname': 'tracerbndSi',   'unit': 'g/m3', 'conversion' : 28.08 / 1000.0},
                                'PON1'       : {'ncvarname': 'phyc',       'bcvarname': 'tracerbndPON1', 'unit': 'g/m3', 'conversion' : 2. * 16. * 14. / (106. * 1000.0)},
                                'POP1'       : {'ncvarname': 'phyc',       'bcvarname': 'tracerbndPOP1', 'unit': 'g/m3', 'conversion' : 2. * 30.97 / (106. * 1000.0)},
                                'POC1'       : {'ncvarname': 'phyc',       'bcvarname': 'tracerbndPOC1', 'unit': 'g/m3', 'conversion' : 2. * 12. / 1000.0},
                                'DON'        : {'ncvarname': 'phyc',       'bcvarname': 'tracerbndDON',  'unit': 'g/m3', 'conversion' : 3.24 * 2. * 16. * 14. / (106. * 1000.0)},
                                'DOP'        : {'ncvarname': 'phyc',       'bcvarname': 'tracerbndDOP',  'unit': 'g/m3', 'conversion' : 1.0 * 2. * 30.97 / (106. * 1000.0)},
                                'DOC'        : {'ncvarname': 'phyc',       'bcvarname': 'tracerbndDOC',  'unit': 'g/m3', 'conversion' : (199. / 20.) * 3.24 * 2. * 16. * 12. / (106. * 1000.0)},
                                'Opal'       : {'ncvarname': 'phyc',       'bcvarname': 'tracerbndOpal', 'unit': 'g/m3', 'conversion' : 0.5 * 0.13 * 28.08 / (1000.0)},
                                'salinity'   : {'ncvarname': 'so',         'bcvarname': 'salinitybnd'},    #'1e-3'
                                'temperature': {'ncvarname': 'thetao',     'bcvarname': 'temperaturebnd'}, #'degC'
                                'ux'         : {'ncvarname': 'uo',         'bcvarname': 'ux' },            #'m/s'
                                'uy'         : {'ncvarname': 'vo',         'bcvarname': 'uy' },            #'m/s'
                                'ux,uy'      : {'ncvarname': 'uo,vo',      'bcvarname': 'ux,uy' },         #'m/s'
                                'steric'     : {'ncvarname': 'zos',        'bcvarname': 'waterlevelbnd'},  #'m'
                                'tide'       : {'ncvarname': '',           'bcvarname': 'waterlevelbnd'},  #'m'
                                }
    conversion_dicts['HYCOM'] = {'salinity'  : {'ncvarname': 'salinity',   'bcvarname': 'salinitybnd'},
                                'temperature': {'ncvarname': 'water_temp', 'bcvarname': 'temperaturebnd'},
                                }
    
    conversion_dict_model = conversion_dicts[model]
    
    return conversion_dict_model


def interpolate_FES(dir_pattern, file_pli, component_list=None, convert_360to180=False, nPoints=None):
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
                         convert_360to180=False,
                         conversion_dict=None, #TODO: alternatively use rename_vars dict and use conversion_dict only for unit conversion. dict containing keys: ncvarname, bcvarname and optionally conversion and unit
                         nPoints=None, #argument for testing
                         reverse_depth=False, #temporary argument to compare easier with old coastserv files
                         ):
    
    if conversion_dict is None:
        conversion_dict_model = get_conversion_dict()
        conversion_dict = conversion_dict_model[quantity]
    ncvarname = conversion_dict['ncvarname'] #rename with origvarname/newvarname
    bcvarname = conversion_dict['bcvarname']
    if ',' in ncvarname:
        raise Exception('ERROR: combined variables not yet supported by hydrolib-core bc writer: https://github.com/Deltares/HYDROLIB-core/issues/316')    
    
    print('initialize ForcingModel()')
    ForcingModel_object = ForcingModel()
        
    file_list_nc = glob.glob(str(dir_pattern))
    print(f'loading mfdataset ({len(file_list_nc)} files with pattern "{dir_pattern.name}")')
    dtstart = dt.datetime.now()
    data_xr = xr.open_mfdataset(file_list_nc,chunks={'time':1}) # can also supply str(dir_pattern) #TODO: does chunks argument solve "PerformanceWarning: Slicing is producing a large chunk."?
    
    #change attributes
    data_xr = data_xr.rename({ncvarname:bcvarname}) #TODO: is this not a bit tricky?
    if refdate_str is not None:
        data_xr.time.encoding['units'] = refdate_str
    
    #get calendar and maybe convert_calendar, makes sure that nc_tstart/nc_tstop are of type pd._libs.tslibs.timestamps.Timestamp
    data_xr_calendar = data_xr['time'].dt.calendar
    if data_xr_calendar != 'proleptic_gregorian': #this is for instance the case in case of noleap (or 365_days) calendars from GFDL and CMCC
        print('WARNING: calendar different than proleptic_gregorian found ({data_xr_calendar}), convert_calendar is called so check output carefully. It should be no issue for datasets with a monthly interval.')
        data_xr = data_xr.convert_calendar('standard') #TODO: does this not result in 29feb nan values in e.g. GFDL model? 
    time_passed = (dt.datetime.now()-dtstart).total_seconds()
    # print(f'>>time passed: {time_passed:.2f} sec')
    
    #get timevar and compare requested dates
    timevar = data_xr['time']
    xr_tstartstop = pd.to_datetime(timevar.isel(time=[0,-1]).to_series())
    nc_tstart = xr_tstartstop.index[0]
    nc_tstop = xr_tstartstop.index[-1]
    if tstart < nc_tstart:
        raise Exception(f'requested tstart {tstart} before nc_tstart {nc_tstart}')
    if tstop > nc_tstop:
        raise Exception(f'requested tstop {tstop} after nc_tstop {nc_tstop}')
    
    if 'latitude' in data_xr.coords:
        lonvarname,latvarname = ['longitude','latitude']
    elif 'lat' in data_xr.coords:
        lonvarname,latvarname = ['lon','lat']
    else:
        raise Exception(f'no lat/lon coords available in file: {data_xr.coords}')
    if convert_360to180: #for FES since it ranges from 0 to 360 instead of -180 to 180 #TODO: make more flexible for models that eg pass -180/+180 crossing (add overlap at lon edges).
        data_xr.coords[lonvarname] = (data_xr.coords[lonvarname] + 180) % 360 - 180
        data_xr = data_xr.sortby(data_xr[lonvarname])
    
    if 'lev' in data_xr[bcvarname].coords: #depth for CMEMS and many others, but lev for GFDL, convert to depth #TODO: provide rename_dict as argument to this function or leave as is?
        data_xr = data_xr.rename({'lev':'depth'}) #TODO: can also do this for data_xr_var only?
        print('variable/coordinate lev renamed to depth')
    
    #retrieve var (after potential longitude conversion) (also selecting relevant times)
    data_xr_var = data_xr[bcvarname].sel(time=slice(tstart,tstop))
    
    if 'latitude' in data_xr_var.coords: #for CMEMS etc #TODO: do rename instead?
        coordname_lat = 'latitude'
        coordname_lon = 'longitude'
    elif 'lat' in data_xr_var.coords:
        coordname_lat = 'lat'
        coordname_lon = 'lon'
    else:
        raise Exception(f'latitude/longitude are not in variable coords: {data_xr_var.coords}. Extend this part of the code for e.g. lat/lon coords')
    
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
        data_interp = data_xr_var.interp({coordname_lat:da_lats, coordname_lon:da_lons}, #also possible without dict: (latitude=da_lats, longitude=da_lons), but this is more flexible
                                         method='linear', 
                                         kwargs={'bounds_error':True}, #error is only raised upon load(), so when the actual value retrieval happens
                                         )
        time_passed = (dt.datetime.now()-dtstart).total_seconds()
        # print(f'>>time passed: {time_passed:.2f} sec')
        
        print('> actual extraction of data from netcdf with .load() (for all PolyObject points at once, so this will take a while)')
        dtstart = dt.datetime.now()
        try:
            datablock_xr_allpoints = data_interp.load() #loading data for all points at once is more efficient compared to loading data per point in loop 
        except ValueError: #generate a proper error with outofbounds requested coordinates, default is "ValueError: One of the requested xi is out of bounds in dimension 0"
            lonvar_vals = data_xr[lonvarname].to_numpy()
            latvar_vals = data_xr[latvarname].to_numpy()
            bool_reqlon_outbounds = (path_lons <= lonvar_vals.min()) | (path_lons >= lonvar_vals.max())
            bool_reqlat_outbounds = (path_lats <= latvar_vals.min()) | (path_lats >= latvar_vals.max())
            reqlatlon_pd = pd.DataFrame({'longitude':path_lons,'latitude':path_lats,'lon outbounds':bool_reqlon_outbounds,'lat outbounds':bool_reqlat_outbounds})
            reqlatlon_pd_outbounds = reqlatlon_pd.loc[bool_reqlon_outbounds | bool_reqlat_outbounds]
            raise ValueError(f'{len(reqlatlon_pd_outbounds)} of requested pli points are out of bounds (valid longitude range {lonvar_vals.min()} to {lonvar_vals.max()}), valid latitude range {latvar_vals.min()} to {latvar_vals.max()}):\n{reqlatlon_pd_outbounds}')
        time_passed = (dt.datetime.now()-dtstart).total_seconds()
        print(f'>>time passed: {time_passed:.2f} sec')
        
        #optional conversion of units and reversing depth dimension
        if 'conversion' in conversion_dict.keys(): #if conversion is present, unit key must also be in conversion_dict
            print(f'> converting units from [{datablock_xr_allpoints.attrs["units"]}] to [{conversion_dict["unit"]}]')
            datablock_xr_allpoints = datablock_xr_allpoints * conversion_dict['conversion'] #conversion drops all attributes of which units (which are changed anyway)
            datablock_xr_allpoints.attrs['units'] = conversion_dict['unit'] #add unit attribute with resulting unit
        
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
            datablock_xr_onepoint.attrs['locationname'] = pli_PolyObject_name_num #TODO: is there a nicer way of passing this data?
            
            if np.isnan(datablock_xr_onepoint.to_numpy()).all(): # check if only nan (out of bounds or land) # we can do .to_numpy() without performance loss, since data is already loaded in datablock_xr_allpoints
                print('WARNING: only nans for this coordinate, this point might be on land')
            if 'depth' in data_xr_var.coords:
                ts_one = DataArray_to_T3D(datablock_xr_onepoint)#,locationname=pli_PolyObject_name_num,refdate_str=refdate_str,bcvarname=bcvarname)
            else:
                ts_one = DataArray_to_TimeSeries(datablock_xr_onepoint)#,locationname=pli_PolyObject_name_num,refdate_str=refdate_str,bcvarname=bcvarname)
            ForcingModel_object.forcing.append(ts_one)
            time_passed = (dt.datetime.now()-dtstart).total_seconds()
            # print(f'>>time passed: {time_passed:.2f} sec')
    
    return ForcingModel_object


