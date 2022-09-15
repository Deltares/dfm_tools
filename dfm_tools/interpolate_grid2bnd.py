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
    conversion_dict = { #TODO: mg/l is the same as g/m3: conversion is phyc in mmol/l to newvar in g/m3
                        'OXY'        : {'ncvarname' : ['o2']       , 'unit': 'g/m3', 'conversion' : 32.0 / 1000.0}, 
                        'NO3'        : {'ncvarname' : ['no3']      , 'unit': 'g/m3', 'conversion' : 14.0 / 1000.0},
                        'PO4'        : {'ncvarname' : ['po4']      , 'unit': 'g/m3', 'conversion' : 30.97 / 1000.0},
                        'Si'         : {'ncvarname' : ['si']       , 'unit': 'g/m3', 'conversion' : 28.08 / 1000.0},
                        'PON1'       : {'ncvarname' : ['phyc']     , 'unit': 'g/m3', 'conversion' : 2. * 16. * 14. / (106. * 1000.0)},
                        'POP1'       : {'ncvarname' : ['phyc']     , 'unit': 'g/m3', 'conversion' : 2. * 30.97 / (106. * 1000.0)},
                        'POC1'       : {'ncvarname' : ['phyc']     , 'unit': 'g/m3', 'conversion' : 2. * 12. / 1000.0},
                        'DON'        : {'ncvarname' : ['phyc']     , 'unit': 'g/m3', 'conversion' : 3.24 * 2. * 16. * 14. / (106. * 1000.0)},
                        'DOP'        : {'ncvarname' : ['phyc']     , 'unit': 'g/m3', 'conversion' : 1.0 * 2. * 30.97 / (106. * 1000.0)},
                        'DOC'        : {'ncvarname' : ['phyc']     , 'unit': 'g/m3', 'conversion' : (199. / 20.) * 3.24 * 2. * 16. * 12. / (106. * 1000.0)},
                        'Opal'       : {'ncvarname' : ['phyc']     , 'unit': 'g/m3', 'conversion' : 0.5 * 0.13 * 28.08 / (1000.0)},
                        'salinity'   : {'ncvarname' : ['so']       , 'bcvarname' : ['salinitybnd']},    #'1e-3'
                        'temperature': {'ncvarname' : ['thetao']   , 'bcvarname' : ['temperaturebnd']}, #'degC'
                        'uxuy'       : {'ncvarname' : ['uo', 'vo'] , 'bcvarname' : ['ux', 'uy'] },      #'m/s'
                        'steric'     : {'ncvarname' : ['zos']      , 'bcvarname' : ['waterlevelbnd']},  #'m'
                        }
    
    return conversion_dict


def interpolate_FES(dir_pattern, file_pli, component_list=None, convert_360to180=False, nPoints=None, debug=False):
    """
    """
    #TODO: resulting amplitudes are slightly different, but the original code might make an 1 indexing mistake? c:\DATA\hydro_tools\FES\PreProcessing_FES_TideModel_imaginary.m
    #TODO: add A0?
    # translate dict from .\hydro_tools\FES\PreProcessing_FES_TideModel_imaginary.m
    translate_dict = {'LA2':'LABDA2',
                      'MTM':'MFM', #Needs to be verified
                      }
    
    print('initialize ForcingModel()')
    ForcingModel_object = ForcingModel()
    
    if component_list is None:
        file_list_nc = glob.glob(str(dir_pattern))
        component_list = [os.path.basename(x).replace('.nc','') for x in file_list_nc] #TODO: add sorting, manually? Add A0? translate dict for component names?
    else:
        file_list_nc = [str(dir_pattern).replace('*',comp) for comp in component_list]
    component_list_upper_pd = pd.Series([x.upper() for x in component_list]).replace(translate_dict, regex=True)
    
    #load boundary file and derive extents for range selection for nc dataset (speeds up process significantly)
    polyfile_object = read_polyfile(file_pli,has_z_values=False)
    pli_PolyObjects = polyfile_object['objects']
    xmin,xmax,ymin,ymax = None,None,None,None
    for iPO, pli_PolyObject_sel in enumerate(pli_PolyObjects):
        path_lons = np.array([point.x for point in pli_PolyObject_sel.points])#[:nPoints]
        path_lats = np.array([point.y for point in pli_PolyObject_sel.points])#[:nPoints]
        xmin = np.min(path_lons.min(),xmax)
        xmax = np.max(path_lons.max(),xmax)
        ymin = np.min(path_lats.min(),ymax)
        ymax = np.max(path_lats.max(),ymax)

    file_list_nc = [str(dir_pattern).replace('*',comp) for comp in component_list]
    #use mfdataset (currently twice as slow as looping over separate datasets, apperantly interp is more difficult but this might be solvable. Also, it does not really matter since we need to compute u/v components of phase anyway, which makes it also slow)
    def extract_component(ds): #TODO: there might be some performance improvement possible. eg by immediately selecting data on lat/lon range coords and then computing u/v
        #https://github.com/pydata/xarray/issues/1380
        if convert_360to180: #for FES since it ranges from 0 to 360 instead of -180 to 180
            ds.coords['lon'] = (ds.coords['lon'] + 180) % 360 - 180
            ds = ds.sortby(ds['lon'])
        ds = ds.sel(lon=slice(xmin-1,xmax+1),lat=slice(ymin-1,ymax+1)) #selection based on pli file extent +-1degree (-180 to +180 degrees convention)
        compname = os.path.basename(ds.encoding["source"]).replace('.nc','')
        compnumber = [component_list.index(compname)]
        ds = ds.assign(compno=compnumber)
        return ds
    
    data_xrsel = xr.open_mfdataset(file_list_nc, combine='nested', concat_dim='compno', preprocess=extract_component)
    lonvar_vals = data_xrsel['lon'].to_numpy() #TODO for extent check (now only selected pli-extent of FES dataset), maybe move outside of point loop since all points are now available
    latvar_vals = data_xrsel['lat'].to_numpy()

    #derive uv phase components (using amplitude=1)
    data_xrsel_phs_rad = np.deg2rad(data_xrsel['phase'])
    #we need to compute u/v components for the phase to avoid zero-crossing interpolation issues
    data_xrsel['phase_u'] = 1*np.cos(data_xrsel_phs_rad)
    data_xrsel['phase_v'] = 1*np.sin(data_xrsel_phs_rad)
    
    data_xr_amp = data_xrsel['amplitude']
    data_xr_phs = data_xrsel['phase']
    data_xr_phs_u = data_xrsel['phase_u']
    data_xr_phs_v = data_xrsel['phase_v']
    
    #load boundary file
    polyfile_object = read_polyfile(file_pli,has_z_values=False)
    
    pli_PolyObjects = polyfile_object['objects']
    for iPO, pli_PolyObject_sel in enumerate(pli_PolyObjects):
        print(f'processing PolyObject {iPO+1} of {len(pli_PolyObjects)}: name={pli_PolyObject_sel.metadata.name}')

        #create DataArrays of lat/lon combinations
        path_lons = np.array([point.x for point in pli_PolyObject_sel.points])[:nPoints]
        path_lats = np.array([point.y for point in pli_PolyObject_sel.points])[:nPoints]
        #make requestedlat/requestedlon DataArrays for proper interpolation in xarray (with new dimension name)
        da_lons = xr.DataArray(path_lons, dims='latloncombi')
        da_lats = xr.DataArray(path_lats, dims='latloncombi')
        
        bool_reqlon_outbounds = (path_lons <= lonvar_vals.min()) | (path_lons >= lonvar_vals.max())
        bool_reqlat_outbounds = (path_lats <= latvar_vals.min()) | (path_lats >= latvar_vals.max())
        if bool_reqlon_outbounds.any():
            raise Exception(f'some of requested lonx outside or on lon bounds ({lonvar_vals.min(),lonvar_vals.max()}):\n{path_lons[bool_reqlon_outbounds]}')
        if bool_reqlat_outbounds.any():
            raise Exception(f'some of requested laty outside or on lat bounds ({latvar_vals.min(),latvar_vals.max()}):\n{path_lats[bool_reqlat_outbounds]}')
        
        #interpolation to lat/lon combinations
        print('> interp mfdataset with all lat/lon coordinates. And actual extraction of data from netcdf with da.as_numpy() (for all PolyObject points at once, so this will take a while)')
        dtstart = dt.datetime.now()
        data_interp_amp_allcoords = data_xr_amp.interp(lat=da_lats,lon=da_lons).as_numpy()/100 #convert from cm to m #TODO: also add assume_sorted and bounds_error arguments
        data_interp_phs_u_allcoords = data_xr_phs_u.interp(lat=da_lats,lon=da_lons).as_numpy()
        data_interp_phs_v_allcoords = data_xr_phs_v.interp(lat=da_lats,lon=da_lons).as_numpy()
        
        time_passed = (dt.datetime.now()-dtstart).total_seconds()
        if debug: print(f'>>time passed: {time_passed:.2f} sec')

        
        for iP, pli_Point_sel in enumerate(pli_PolyObject_sel.points[:nPoints]): #looping over plipoints within component loop, append to datablock_pd_allcomp
            print(f'processing Point {iP+1} of {len(pli_PolyObject_sel.points)}: ',end='')
            lonx, laty = pli_Point_sel.x, pli_Point_sel.y
            print(f'(x={lonx}, y={laty})')
            pli_PolyObject_name_num = f'{pli_PolyObject_sel.metadata.name}_{iP+1:04d}'
            
            data_interp_amp = data_interp_amp_allcoords.sel(latloncombi=iP)
            data_interp_phs_u = data_interp_phs_u_allcoords.sel(latloncombi=iP)
            data_interp_phs_v = data_interp_phs_v_allcoords.sel(latloncombi=iP)
            
            data_interp_phs = np.rad2deg(np.arctan2(data_interp_phs_v,data_interp_phs_u)) #np.arctan2(y,x)
            datablock_pd_allcomp = pd.DataFrame({'component':component_list_upper_pd,'amp':data_interp_amp.to_numpy(),'phs':data_interp_phs.to_numpy()})
            if datablock_pd_allcomp['amp'].isnull().any():
                print('WARNING: only nans for this coordinate') #TODO: this can happen on land, raise exception or warning?
            datablock_list = datablock_pd_allcomp.values.tolist()
            
            # Each .bc file can contain 1 or more timeseries, one for each support point:
            print('> constructing TimeSeries and appending to ForcingModel()')
            ts_one = Astronomic(name=pli_PolyObject_name_num,
                                quantityunitpair=[QuantityUnitPair(quantity="astronomic component", unit='-'),
                                                  QuantityUnitPair(quantity='waterlevelbnd amplitude', unit='m'),#unit=data_xr_amp.attrs['units']),
                                                  QuantityUnitPair(quantity='waterlevelbnd phase', unit=data_xr_phs.attrs['units'])],
                                datablock=datablock_list, 
                                )
            
            ForcingModel_object.forcing.append(ts_one)        
    return ForcingModel_object


def interpolate_nc_to_bc(dir_pattern, file_pli, quantity, 
                         tstart, tstop, refdate_str, 
                         convert_360to180=False,
                         nPoints=None, debug=False):
    
    """
    nPolyObjects = None #None gives all PolyObjects
    nPoints = 2 #None gives all Points in PolyObject
    """
    nPolyObjects = None
    
    conversion_dict = get_conversion_dict()
    varname_file = conversion_dict[quantity]['ncvarname'][0] #TODO: [1] is also necessary for uxuy
    
    print('initialize ForcingModel()')
    ForcingModel_object = ForcingModel()
        
    file_list_nc = glob.glob(str(dir_pattern))
    print(f'loading mfdataset ({len(file_list_nc)} files with pattern "{dir_pattern.name}")')
    dtstart = dt.datetime.now()
    data_xr = xr.open_mfdataset(file_list_nc)# TODO: can also supply str(dir_pattern)
    
    #get calendar and maybe convert_calendar, makes sure that nc_tstart/nc_tstop are of type pd._libs.tslibs.timestamps.Timestamp
    data_xr_calendar = data_xr['time'].dt.calendar
    if data_xr_calendar != 'proleptic_gregorian': #this is for instance the case in case of noleap (or 365_days) calendars from GFDL and CMCC
        print('WARNING: calendar different than proleptic_gregorian found ({data_xr_calendar}), convert_calendar is called so check output carefully. It should be no issue for datasets with a monthly interval.')
        data_xr = data_xr.convert_calendar('standard') #TODO: does this not result in 29feb nan values in e.g. ☻GFDL model? 
    time_passed = (dt.datetime.now()-dtstart).total_seconds()
    if debug: print(f'>>time passed: {time_passed:.2f} sec')
    
    #get timevar and compare requested dates
    timevar = data_xr['time'] #TODO: mabye only request times on time selection (or first index times and then convert to series)
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
        print(data_xr)
        raise Exception(f'no lat/lon coords available in file: {data_xr.coords}')
    if convert_360to180: #for FES since it ranges from 0 to 360 instead of -180 to 180 #TODO: make more flexible for models that eg pass -180/+180 crossing (add overlap at lon edges).
        data_xr.coords[lonvarname] = (data_xr.coords[lonvarname] + 180) % 360 - 180
        data_xr = data_xr.sortby(data_xr[lonvarname])
    lonvar_vals = data_xr[lonvarname].to_numpy()
    latvar_vals = data_xr[latvarname].to_numpy()
    
    #retrieve var (after potential longitude conversion) (also selecting relevant times)
    data_xr_var = data_xr[varname_file].sel(time=slice(tstart,tstop))
    
    #check if depth coordinate is present in variable (not only in file)
    if 'depth' in data_xr_var.coords: #depth for CMEMS and many others
        has_depth = True
        depthvarname = 'depth'
    elif 'lev' in data_xr_var.coords: #lev for GFDL
        has_depth = True
        depthvarname = 'lev'
    else:
        has_depth = False
    
    if has_depth: #get depth variable and values
        vardepth = data_xr_var[depthvarname]
        depth_array = vardepth.to_numpy()[::-1] #TODO: flip array is not necessary, check datablock comment
        if vardepth.attrs['positive'] == 'down': #attribute appears in CMEMS, GFDL and CMCC, save to assume presence?
            depth_array = -depth_array
    
    #load boundary file
    #polyfile_object = PolyFile(file_pli,has_z_values=False) #TODO ISFIXED: should work with hydrolib-core>0.3.0. also without has_z_values argument
    polyfile_object = read_polyfile(file_pli,has_z_values=False) #TODO: this warning can be suppressed (or how to fix): "UserWarning: White space at the start of the line is ignored." https://github.com/Deltares/HYDROLIB-core/issues/320
    
    pli_PolyObjects = polyfile_object['objects']
    for iPO, pli_PolyObject_sel in enumerate(pli_PolyObjects[:nPolyObjects]):
        print(f'processing PolyObject {iPO+1} of {len(pli_PolyObjects)}: name={pli_PolyObject_sel.metadata.name}')
        
        #create requestedlat/requestedlon DataArrays for proper interpolation in xarray (with new dimension name)
        path_lons = np.array([point.x for point in pli_PolyObject_sel.points])[:nPoints]
        path_lats = np.array([point.y for point in pli_PolyObject_sel.points])[:nPoints]
        da_lons = xr.DataArray(path_lons, dims='latloncombi')
        da_lats = xr.DataArray(path_lats, dims='latloncombi')
        
        #check whether coordinates are in bounds. This is also done with da.interp(kwargs={'bounds_error':True}), but this is to give proper feedback about which coordinates are out of bounds
        bool_reqlon_outbounds = (path_lons <= lonvar_vals.min()) | (path_lons >= lonvar_vals.max())
        bool_reqlat_outbounds = (path_lats <= latvar_vals.min()) | (path_lats >= latvar_vals.max())
        if bool_reqlon_outbounds.any():
            raise Exception(f'some of requested lonx outside or on lon bounds ({lonvar_vals.min(),lonvar_vals.max()}):\n{path_lons[bool_reqlon_outbounds]}')
        if bool_reqlat_outbounds.any():
            raise Exception(f'some of requested laty outside or on lat bounds ({latvar_vals.min(),latvar_vals.max()}):\n{path_lats[bool_reqlat_outbounds]}')
        
        #interpolation to lat/lon combinations
        print('> interp mfdataset to all PolyObject points (lat/lon coordinates)')
        dtstart = dt.datetime.now()
        if 'latitude' in data_xr_var.coords: #for CMEMS etc
            coordname_lat = 'latitude'
            coordname_lon = 'longitude'
        elif 'lat' in data_xr_var.coords:
            coordname_lat = 'lat'
            coordname_lon = 'lon'
        else: #TODO: lat/latitude (should be flexible instead of hardcoded) 
            raise Exception(f'latitude/longitude are not in variable coords: {data_xr_var.coords}. Extend this part of the code for e.g. lat/lon coords')
        data_interp = data_xr_var.interp({coordname_lat:da_lats, coordname_lon:da_lons}, #also possible without dict: (latitude=da_lats, longitude=da_lons), but this is more flexible
                                         method='linear', 
                                         kwargs={'bounds_error':True}, #error is only raised upon as_numpy() or to_numpy() (when the actual value retrieval happens)
                                         assume_sorted=True, #TODO: assume_sorted increases performance?
                                         )
        time_passed = (dt.datetime.now()-dtstart).total_seconds()
        if debug: print(f'>>time passed: {time_passed:.2f} sec')
        
        print('> actual extraction of data from netcdf with da.as_numpy() (for all PolyObject points at once, so this will take a while)')
        dtstart = dt.datetime.now()
        #try:
        datablock_raw_allcoords = data_interp.as_numpy() #TODO: as_numpy() is maybe not the best method, but at least all data should be read now and we still want an object where sel works on
        #except ValueError as e:
        #    raise ValueError(f'{e}\nlons (valid range {lonvar_vals.min()},{lonvar_vals.max()}):\n{path_lons[bool_reqlon_outbounds]}\nlats (valid range {latvar_vals.min()},{latvar_vals.max()}):\n{path_lats[bool_reqlat_outbounds]}')
        time_passed = (dt.datetime.now()-dtstart).total_seconds()
        #if debug:
        print(f'>>time passed: {time_passed:.2f} sec')
        
        for iP, pli_Point_sel in enumerate(pli_PolyObject_sel.points[:nPoints]):
            print(f'processing Point {iP+1} of {len(pli_PolyObject_sel.points)}: ',end='')
            lonx_print, laty_print = pli_Point_sel.x, pli_Point_sel.y
            print(f'(x={lonx_print}, y={laty_print})')
            pli_PolyObject_name_num = f'{pli_PolyObject_sel.metadata.name}_{iP+1:04d}'
                        
            print('> selecting data for current coord from numpy array')
            dtstart = dt.datetime.now()
            datablock_xr = datablock_raw_allcoords.isel(latloncombi=iP)
            datablock_raw = datablock_xr.to_numpy()
            time_passed = (dt.datetime.now()-dtstart).total_seconds()
            if debug: print(f'>>time passed: {time_passed:.2f} sec')
            
            # check if only nan (out of bounds or land):
            if np.isnan(datablock_raw).all():
                print('WARNING: only nan values for this coordinate') #TODO: this can happen on land, raise exception or warning?

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
            
            print('> ffill nans, converting units and concatenating time column')
            dtstart = dt.datetime.now()
            if has_depth:
                datablock = pd.DataFrame(datablock_raw).fillna(method='ffill',axis=1).values #fill nans forward (corresponds to vertical extrapolation for CMEMS) #TODO: make depth axis flexible
                #if debug: print(datablock_raw), print(datablock)
                datablock = datablock[:,::-1] #flipping axis #TODO: this assumes depth as second dimension and that might not be true for other models. Flipping depth_vals and datablock is not required by kernel, so remove after validation is complete
            else:
                datablock = datablock_raw[:,np.newaxis]
            
            #conversion of units etc
            if 'conversion' in conversion_dict[quantity].keys():
                datablock = datablock * conversion_dict[quantity]['conversion']
            if 'bcvarname' in conversion_dict[quantity].keys():
                bcvarname = conversion_dict[quantity]['bcvarname'][0] #TODO: [1] is necessary for uxuy
            else: #TODO: only works for waq tracers, so add check or add to conversion_dict (latter also makes sense)
                bcvarname = f'tracerbnd{quantity}'
                #bcvarname = quantity
            if 'unit' in conversion_dict[quantity].keys():
                varunit = conversion_dict[quantity]['unit']
            else:
                varunit = data_xr_var.attrs['units']
            
            timevar_sel = datablock_xr.time
            timevar_sel_rel = date2num(pd.DatetimeIndex(timevar_sel.to_numpy()).to_pydatetime(),units=refdate_str,calendar='standard')
            datablock_incltime = np.concatenate([timevar_sel_rel[:,np.newaxis],datablock],axis=1)
            time_passed = (dt.datetime.now()-dtstart).total_seconds()
            if debug: print(f'>>time passed: {time_passed:.2f} sec')
            
            # Each .bc file can contain 1 or more timeseries, one for each support point:
            print('> constructing TimeSeries and appending to ForcingModel()')
            dtstart = dt.datetime.now()
            if has_depth:
                list_QUP_perlayer = [QuantityUnitPair(quantity=bcvarname, unit=varunit) for iL in range(vardepth.size)] #TODO: verticalposition 1/2/3/n is not supported. https://github.com/Deltares/HYDROLIB-core/issues/317
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
                                                      QuantityUnitPair(quantity=bcvarname, unit=varunit)],
                                    timeinterpolation=TimeInterpolation.linear,
                                    datablock=datablock_incltime.tolist(), 
                                    )
            
            ForcingModel_object.forcing.append(ts_one)
            time_passed = (dt.datetime.now()-dtstart).total_seconds()
            if debug: print(f'>>time passed: {time_passed:.2f} sec')
    
    return ForcingModel_object


