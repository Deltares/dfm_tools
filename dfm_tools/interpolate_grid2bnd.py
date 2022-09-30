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
from netCDF4 import date2num
import matplotlib.pyplot as plt

from hydrolib.core.io.bc.models import (
    ForcingModel,
    QuantityUnitPair,
    VerticalPositionType,
    TimeInterpolation,
    T3D,
    TimeSeries,
    Astronomic,
)
from hydrolib.core.io.polyfile.models import (
    #Description,
    #Metadata,
    #Point,
    PolyFile,
    #PolyObject,
)
from hydrolib.core.io.polyfile.parser import read_polyfile #TODO: should be replaced with PolyFile above


def get_conversion_dict():
    # conversion_dict, contains ncvarname as array because uxuy relies on 2 CMEMS variables
    conversion_dict = { # mg/l is the same as g/m3: conversion is phyc in mmol/l to newvar in g/m3
                        'OXY'        : {'ncvarname': 'o2',      'bcvarname': 'tracerbndOXY',  'unit': 'g/m3', 'conversion' : 32.0 / 1000.0}, 
                        'NO3'        : {'ncvarname': 'no3',     'bcvarname': 'tracerbndNO3',  'unit': 'g/m3', 'conversion' : 14.0 / 1000.0},
                        'PO4'        : {'ncvarname': 'po4',     'bcvarname': 'tracerbndPO4',  'unit': 'g/m3', 'conversion' : 30.97 / 1000.0},
                        'Si'         : {'ncvarname': 'si',      'bcvarname': 'tracerbndSi',   'unit': 'g/m3', 'conversion' : 28.08 / 1000.0},
                        'PON1'       : {'ncvarname': 'phyc',    'bcvarname': 'tracerbndPON1', 'unit': 'g/m3', 'conversion' : 2. * 16. * 14. / (106. * 1000.0)},
                        'POP1'       : {'ncvarname': 'phyc',    'bcvarname': 'tracerbndPOP1', 'unit': 'g/m3', 'conversion' : 2. * 30.97 / (106. * 1000.0)},
                        'POC1'       : {'ncvarname': 'phyc',    'bcvarname': 'tracerbndPOC1', 'unit': 'g/m3', 'conversion' : 2. * 12. / 1000.0},
                        'DON'        : {'ncvarname': 'phyc',    'bcvarname': 'tracerbndDON',  'unit': 'g/m3', 'conversion' : 3.24 * 2. * 16. * 14. / (106. * 1000.0)},
                        'DOP'        : {'ncvarname': 'phyc',    'bcvarname': 'tracerbndDOP',  'unit': 'g/m3', 'conversion' : 1.0 * 2. * 30.97 / (106. * 1000.0)},
                        'DOC'        : {'ncvarname': 'phyc',    'bcvarname': 'tracerbndDOC',  'unit': 'g/m3', 'conversion' : (199. / 20.) * 3.24 * 2. * 16. * 12. / (106. * 1000.0)},
                        'Opal'       : {'ncvarname': 'phyc',    'bcvarname': 'tracerbndOpal', 'unit': 'g/m3', 'conversion' : 0.5 * 0.13 * 28.08 / (1000.0)},
                        'salinity'   : {'ncvarname': 'so',      'bcvarname': 'salinitybnd'},    #'1e-3'
                        'temperature': {'ncvarname': 'thetao',  'bcvarname': 'temperaturebnd'}, #'degC'
                        'ux'         : {'ncvarname': 'uo',      'bcvarname': 'ux' },            #'m/s'
                        'uy'         : {'ncvarname': 'vo',      'bcvarname': 'uy' },            #'m/s'
                        'ux,uy'      : {'ncvarname': 'uo,vo',   'bcvarname': 'ux,uy' },         #'m/s'
                        'steric'     : {'ncvarname': 'zos',     'bcvarname': 'waterlevelbnd'},  #'m'
                        'tide'       : {'ncvarname': '',        'bcvarname': 'waterlevelbnd'},  #'m'
                        }
    
    return conversion_dict


def interpolate_FES(dir_pattern, file_pli, component_list=None, convert_360to180=False, nPoints=None, debug=False):
    """
    """
    #TODO: resulting amplitudes are slightly different, but the original code might make an 1 indexing mistake? c:\DATA\hydro_tools\FES\PreProcessing_FES_TideModel_imaginary.m
    # translate dict from .\hydro_tools\FES\PreProcessing_FES_TideModel_imaginary.m
    translate_dict = {'LA2':'LABDA2',
                      'MTM':'MFM', #Needs to be verified
                      }
    
    print('initialize ForcingModel()')
    ForcingModel_object = ForcingModel()
    
    if component_list is None:
        file_list_nc = glob.glob(str(dir_pattern))
        component_list = [os.path.basename(x).replace('.nc','') for x in file_list_nc]
    else:
        file_list_nc = [str(dir_pattern).replace('*',comp) for comp in component_list]
    component_list_upper_pd = pd.Series([x.upper() for x in component_list]).replace(translate_dict, regex=True)
    
    def extract_component(ds):
        #https://github.com/pydata/xarray/issues/1380
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
    polyfile_object = read_polyfile(file_pli,has_z_values=False)
    
    pli_PolyObjects = polyfile_object['objects']
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
        if debug: print(f'>>time passed: {time_passed:.2f} sec')
        
        
        print('> actual extraction of data from netcdf with .load() (for all PolyObject points at once, so this will take a while)')
        dtstart = dt.datetime.now()
        try:
            data_interp_amp_allcoords = data_interp_allcoords['amplitude'].load()/100 #convert from cm to m
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
        #if debug:
        print(f'>>time passed: {time_passed:.2f} sec')

        
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
                                                  QuantityUnitPair(quantity='waterlevelbnd amplitude', unit='m'),#unit=data_xr_amp.attrs['units']),
                                                  QuantityUnitPair(quantity='waterlevelbnd phase', unit=data_xrsel['phase'].attrs['units'])],
                                datablock=datablock_list, 
                                )
            
            ForcingModel_object.forcing.append(ts_one)        
    return ForcingModel_object


def interpolate_nc_to_bc(dir_pattern, file_pli, quantity, 
                         tstart, tstop, refdate_str, 
                         convert_360to180=False,
                         nPoints=None, debug=False,
                         reverse_depth=False): #temporary argument to compare easier with old coastserv files
    
    conversion_dict = get_conversion_dict()
    ncvarname = conversion_dict[quantity]['ncvarname']
    bcvarname = conversion_dict[quantity]['bcvarname']
    if ',' in ncvarname:
        raise Exception('ERROR: combined variables not yet supported by hydrolib-core bc writer: https://github.com/Deltares/HYDROLIB-core/issues/316')
    
    print('initialize ForcingModel()')
    ForcingModel_object = ForcingModel()
        
    file_list_nc = glob.glob(str(dir_pattern))
    print(f'loading mfdataset ({len(file_list_nc)} files with pattern "{dir_pattern.name}")')
    dtstart = dt.datetime.now()
    data_xr = xr.open_mfdataset(file_list_nc) # can also supply str(dir_pattern)
    
    #get calendar and maybe convert_calendar, makes sure that nc_tstart/nc_tstop are of type pd._libs.tslibs.timestamps.Timestamp
    data_xr_calendar = data_xr['time'].dt.calendar
    if data_xr_calendar != 'proleptic_gregorian': #this is for instance the case in case of noleap (or 365_days) calendars from GFDL and CMCC
        print('WARNING: calendar different than proleptic_gregorian found ({data_xr_calendar}), convert_calendar is called so check output carefully. It should be no issue for datasets with a monthly interval.')
        data_xr = data_xr.convert_calendar('standard') #TODO: does this not result in 29feb nan values in e.g. GFDL model? 
    time_passed = (dt.datetime.now()-dtstart).total_seconds()
    if debug: print(f'>>time passed: {time_passed:.2f} sec')
    
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
    
    #retrieve var (after potential longitude conversion) (also selecting relevant times)
    data_xr_var = data_xr[ncvarname].sel(time=slice(tstart,tstop))
    
    #check if depth coordinate is present in variable (not only in file)
    if 'depth' in data_xr_var.coords: #depth for CMEMS and many others
        has_depth = True
        depthvarname = 'depth'
    elif 'lev' in data_xr_var.coords: #lev for GFDL
        has_depth = True
        depthvarname = 'lev'
    else:
        has_depth = False
    
    #load boundary file
    #polyfile_object = PolyFile(file_pli,has_z_values=False) #TODO ISFIXED: should work with hydrolib-core>0.3.0. also without has_z_values argument
    polyfile_object = read_polyfile(file_pli,has_z_values=False) #TODO: this warning can be suppressed (or how to fix): "UserWarning: White space at the start of the line is ignored." https://github.com/Deltares/HYDROLIB-core/issues/320
    
    pli_PolyObjects = polyfile_object['objects']
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
        if 'latitude' in data_xr_var.coords: #for CMEMS etc
            coordname_lat = 'latitude'
            coordname_lon = 'longitude'
        elif 'lat' in data_xr_var.coords:
            coordname_lat = 'lat'
            coordname_lon = 'lon'
        else:
            raise Exception(f'latitude/longitude are not in variable coords: {data_xr_var.coords}. Extend this part of the code for e.g. lat/lon coords')
        data_interp = data_xr_var.interp({coordname_lat:da_lats, coordname_lon:da_lons}, #also possible without dict: (latitude=da_lats, longitude=da_lons), but this is more flexible
                                         method='linear', 
                                         kwargs={'bounds_error':True}, #error is only raised upon load(), so when the actual value retrieval happens
                                         )
        time_passed = (dt.datetime.now()-dtstart).total_seconds()
        if debug: print(f'>>time passed: {time_passed:.2f} sec')
        
        print('> actual extraction of data from netcdf with .load() (for all PolyObject points at once, so this will take a while)')
        dtstart = dt.datetime.now()
        try:
            datablock_raw_allcoords = data_interp.load() #loading data for all points at once is more efficient compared to loading data per point in loop 
        except ValueError: #generate a proper error with outofbounds requested coordinates, default is "ValueError: One of the requested xi is out of bounds in dimension 0"
            lonvar_vals = data_xr[lonvarname].to_numpy()
            latvar_vals = data_xr[latvarname].to_numpy()
            bool_reqlon_outbounds = (path_lons <= lonvar_vals.min()) | (path_lons >= lonvar_vals.max())
            bool_reqlat_outbounds = (path_lats <= latvar_vals.min()) | (path_lats >= latvar_vals.max())
            reqlatlon_pd = pd.DataFrame({'longitude':path_lons,'latitude':path_lats,'lon outbounds':bool_reqlon_outbounds,'lat outbounds':bool_reqlat_outbounds})
            reqlatlon_pd_outbounds = reqlatlon_pd.loc[bool_reqlon_outbounds | bool_reqlat_outbounds]
            raise ValueError(f'{len(reqlatlon_pd_outbounds)} of requested pli points are out of bounds (valid longitude range {lonvar_vals.min()} to {lonvar_vals.max()}), valid latitude range {latvar_vals.min()} to {latvar_vals.max()}):\n{reqlatlon_pd_outbounds}')
        time_passed = (dt.datetime.now()-dtstart).total_seconds()
        #if debug:
        print(f'>>time passed: {time_passed:.2f} sec')
        
        print('> optional conversion of units and reversing of depth dimension')
        #optional conversion of units and reversing depth dimension
        if 'conversion' in conversion_dict[quantity].keys(): #if conversion is present, unit key must also be in conversion_dict
            print('converting units')
            datablock_raw_allcoords = datablock_raw_allcoords * conversion_dict[quantity]['conversion'] #conversion drops all attributes of which units (which are changed anyway)
            datablock_raw_allcoords.attrs['units'] = conversion_dict[quantity]['unit'] #add unit attribute with resulting unit
        
        if has_depth & reverse_depth:
            datablock_raw_allcoords = datablock_raw_allcoords.reindex({depthvarname:list(reversed(datablock_raw_allcoords[depthvarname]))})
            
        for iP, pli_Point_sel in enumerate(pli_PolyObject_sel.points[:nPoints]):
            print(f'processing Point {iP+1} of {len(pli_PolyObject_sel.points)}: ',end='')
            lonx_print, laty_print = pli_Point_sel.x, pli_Point_sel.y
            print(f'(x={lonx_print}, y={laty_print})')
            pli_PolyObject_name_num = f'{pli_PolyObject_sel.metadata.name}_{iP+1:04d}'
            
            
            #TODO: steps from here to T3D/Timeseries could be bound method of those objects? (eg T3D.from_xarray_dataarray(datablock_xr)). Required information:
            # has_depth = 1 # can be avoided if separate for Timeseries and T3D, or inferred
            # depthvarname = 1 # can be avoided if overwritten with 'depth'
            # pli_PolyObject_name_num = 1
            # refdate_str = 1
            # bcvarname = 1
            # datablock_xr = 1
            print('> select data for point, ffill nans and concatenating time column')
            dtstart = dt.datetime.now()
            datablock_xr = datablock_raw_allcoords.isel(latloncombi=iP)
            if has_depth:
                #get depth variable and values
                depth_array = datablock_xr[depthvarname].to_numpy()
                if datablock_xr[depthvarname].attrs['positive'] == 'down': #attribute appears in CMEMS, GFDL and CMCC, save to assume presence?
                    depth_array = -depth_array
                #ffill data
                datablock_xr = datablock_xr.bfill(dim=depthvarname).ffill(dim=depthvarname) #fill nans back and forward (corresponds to vertical extrapolation for CMEMS). Deep values are filled if order is shallow to deep
                datablock_np = datablock_xr.to_numpy()
            else:
                datablock_np = datablock_xr.to_numpy()[:,np.newaxis]
                
            # check if only nan (out of bounds or land):
            if np.isnan(datablock_np).all():
                print('WARNING: only nans for this coordinate, this point might be on land')
            
            timevar_sel = datablock_xr.time
            timevar_sel_rel = date2num(pd.DatetimeIndex(timevar_sel.to_numpy()).to_pydatetime(),units=refdate_str,calendar='standard')
            
            datablock_incltime = np.concatenate([timevar_sel_rel[:,np.newaxis],datablock_np],axis=1)
            time_passed = (dt.datetime.now()-dtstart).total_seconds()
            if debug: print(f'>>time passed: {time_passed:.2f} sec')
            
            
            # Each .bc file can contain 1 or more timeseries, in this case one for each support point
            print('> constructing TimeSeries and appending to ForcingModel()')
            dtstart = dt.datetime.now()
            if has_depth:
                list_QUP_perlayer = [QuantityUnitPair(quantity=bcvarname, unit=datablock_xr.attrs['units']) for iL in range(datablock_xr[depthvarname].size)] #TODO: verticalposition 1/2/3/n is not supported. https://github.com/Deltares/HYDROLIB-core/issues/317
                ts_one = T3D(name=pli_PolyObject_name_num,
                             verticalpositions=depth_array.tolist(), #TODO: should be "Vertical position specification = [..]" but is verticalPositions = [..]" (both possible?). https://github.com/Deltares/HYDROLIB-core/issues/317
                             verticalInterpolation='linear',
                             verticalPositionType=VerticalPositionType('ZBed'), #TODO: should be "Vertical position type = zdatum" but is "verticalPositionType = ZBed" (zdatum is niet beschikbaar). https://github.com/Deltares/HYDROLIB-core/issues/317
                             quantityunitpair=[QuantityUnitPair(quantity="time", unit=refdate_str)]+list_QUP_perlayer,
                             timeinterpolation=TimeInterpolation.linear, #TODO: not passed on to bc file. https://github.com/Deltares/HYDROLIB-core/issues/317
                             datablock=datablock_incltime.tolist(), #TODO: numpy array is not supported by TimeSeries. https://github.com/Deltares/HYDROLIB-core/issues/322
                             )
            else:
                ts_one = TimeSeries(name=pli_PolyObject_name_num,
                                    verticalposition=VerticalPositionType('ZBed'), #TODO: is not passed on to bc file and that makes sense, but it should raise error since it is not relevant for timeseries. https://github.com/Deltares/HYDROLIB-core/issues/321
                                    quantityunitpair=[QuantityUnitPair(quantity="time", unit=refdate_str),
                                                      QuantityUnitPair(quantity=bcvarname, unit=datablock_xr.attrs['units'])],
                                    timeinterpolation=TimeInterpolation.linear,
                                    datablock=datablock_incltime.tolist(), 
                                    )
            
            ForcingModel_object.forcing.append(ts_one)
            time_passed = (dt.datetime.now()-dtstart).total_seconds()
            if debug: print(f'>>time passed: {time_passed:.2f} sec')
            
            
            if debug:
                print('> plotting')
                dtstart = dt.datetime.now()
                fig,ax = plt.subplots(figsize=(10,7))
                datablock_xr.T.plot() #uses plot.line() for 1D arrays and plot.pcolormesh() for 2D arrays: https://docs.xarray.dev/en/stable/generated/xarray.DataArray.plot.html
                ax.set_title(f'{quantity} {pli_PolyObject_name_num}')
                if 'depth' in data_xr_var.coords:
                    ax.set_ylim(0,200)
                    ax.invert_yaxis()
                fig.tight_layout()
                time_passed = (dt.datetime.now()-dtstart).total_seconds()
                if debug: print(f'>>time passed: {time_passed:.2f} sec')
            
    
    return ForcingModel_object


