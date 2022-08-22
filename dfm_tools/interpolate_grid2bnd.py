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
-	downloading part is not included (Bjorn improved that part so add to dfm_tools?) 
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


def interpolate_FES(dir_pattern, file_pli, convert_360to180=False, nPoints=None, debug=False):
    """
    """
    #TODO: resulting amplitudes are slightly different, but the original code might make an 1 indexing mistake? c:\DATA\hydro_tools\FES\PreProcessing_FES_TideModel_imaginary.m

    # translate dict from .\hydro_tools\FES\PreProcessing_FES_TideModel_imaginary.m
    translate_dict = {'LA2':'LABDA2',
                      'MTM':'MFM', #Needs to be verified
                      }
    
    print('initialize ForcingModel()')
    ForcingModel_object = ForcingModel()
    
    file_list_nc = glob.glob(str(dir_pattern))
    component_list = [os.path.basename(x).replace('.nc','') for x in file_list_nc] #TODO: add sorting, manually? Add A0? translate dict for component names?
    component_list_upper_pd = pd.Series([x.upper() for x in component_list]).replace(translate_dict, regex=True)
    
    #use mfdataset (currently twice as slow as looping over separate datasets, apperantly interp is more difficult but this might be solvable. Also, it does not really matter since we need to compute u/v components of phase anyway, which makes it also slow)
    def extract_component(ds): #TODO: there might be some performance improvement possible
        #https://github.com/pydata/xarray/issues/1380
        compnumber = [0]
        ds = ds.assign(compno=compnumber)
        data_xr_phs_rad = np.deg2rad(ds['phase'])
        #we need to compute u/v components for the phase to avoid zero-crossing interpolation issues
        ds['phase_u'] = 1*np.cos(data_xr_phs_rad)
        ds['phase_v'] = 1*np.sin(data_xr_phs_rad)
        return ds
    
    data_xrall = xr.open_mfdataset(file_list_nc, combine='nested', concat_dim='compno', preprocess=extract_component)
    if convert_360to180: #for FES since it ranges from 0 to 360 instead of -180 to 180
        data_xrall.coords['lon'] = (data_xrall.coords['lon'] + 180) % 360 - 180
        data_xrall = data_xrall.sortby(data_xrall['lon'])
    
    lonvar_vals = data_xrall['lon'].to_numpy()
    latvar_vals = data_xrall['lat'].to_numpy()
    
    data_xr_amp = data_xrall['amplitude']
    data_xr_phs = data_xrall['phase']
    data_xr_phs_u = data_xrall['phase_u']
    data_xr_phs_v = data_xrall['phase_v']
    
    #load boundary file
    polyfile_object = read_polyfile(file_pli,has_z_values=False)
    
    pli_PolyObjects = polyfile_object['objects']
    for iPO, pli_PolyObject_sel in enumerate(pli_PolyObjects):
        print(f'processing PolyObject {iPO+1} of {len(pli_PolyObjects)}: name={pli_PolyObject_sel.metadata.name}')
        
        for iP, pli_Point_sel in enumerate(pli_PolyObject_sel.points[:nPoints]): #looping over plipoints within component loop, append to datablock_pd_allcomp
            print(f'processing Point {iP+1} of {len(pli_PolyObject_sel.points)}: ',end='')
            lonx, laty = pli_Point_sel.x, pli_Point_sel.y
            print(f'(x={lonx}, y={laty})')
            pli_PolyObject_name_num = f'{pli_PolyObject_sel.metadata.name}_{iP+1:04d}'
            
            if (lonx <= lonvar_vals.min()) or (lonx >= lonvar_vals.max()):
                raise Exception(f'requested lonx {lonx} outside or on lon bounds ({lonvar_vals.min(),lonvar_vals.max()})')
            if (laty <= latvar_vals.min()) or (laty >= latvar_vals.max()):
                raise Exception(f'requested laty {laty} outside or on lat bounds ({latvar_vals.min(),latvar_vals.max()})')
                    
            data_interp_amp = data_xr_amp.interp(lat=laty,lon=lonx)/100 #convert from cm to m
            data_interp_phs_u = data_xr_phs_u.interp(lat=laty,lon=lonx)
            data_interp_phs_v = data_xr_phs_v.interp(lat=laty,lon=lonx)
            
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
    data_xr = xr.open_mfdataset(file_list_nc)#, combine='by_coords', decode_times=False)
    
    #get calendar and maybe convert_calendar, makes sure that nc_tstart/nc_tstop are of type pd._libs.tslibs.timestamps.Timestamp
    data_xr_calendar = data_xr['time'].dt.calendar
    if data_xr_calendar != 'proleptic_gregorian': #this is for instance the case in case of noleap (or 365_days) calendars from GFDL and CMCC
        print('WARNING: calendar different than proleptic_gregorian found ({data_xr_calendar}), convert_calendar is called so check output carefully. It should be no issue for datasets with a monthly interval.')
        data_xr = data_xr.convert_calendar('standard')
    time_passed = (dt.datetime.now()-dtstart).total_seconds()
    if debug: print(f'>>time passed: {time_passed:.2f} sec')
    
    #get timevar and compare requested dates
    timevar = data_xr['time']
    nc_tstart = pd.to_datetime(timevar.to_series().index[0])
    nc_tstop = pd.to_datetime(timevar.to_series().index[-1])
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
    
    #retrieve var (after potential longitude conversion)
    data_xr_var = data_xr[varname_file]
    
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
        
        for iP, pli_Point_sel in enumerate(pli_PolyObject_sel.points[:nPoints]):
            print(f'processing Point {iP+1} of {len(pli_PolyObject_sel.points)}: ',end='')
    
            lonx, laty = pli_Point_sel.x, pli_Point_sel.y
            print(f'(x={lonx}, y={laty})')
            if (lonx <= lonvar_vals.min()) or (lonx >= lonvar_vals.max()):
                raise Exception(f'requested lonx {lonx} outside or on lon bounds ({lonvar_vals.min(),lonvar_vals.max()})')
            if (laty <= latvar_vals.min()) or (laty >= latvar_vals.max()):
                raise Exception(f'requested laty {laty} outside or on lat bounds ({latvar_vals.min(),latvar_vals.max()})')
            
            pli_PolyObject_name_num = f'{pli_PolyObject_sel.metadata.name}_{iP+1:04d}'
            
            print('> interp mfdataset with lat/lon coordinates')
            dtstart = dt.datetime.now()
            if 'latitude' in data_xr_var.coords:
                data_interp_alltimes = data_xr_var.interp(latitude=laty,longitude=lonx,
                                                          #kwargs={'bounds_error':True}, #TODO: bounds_error is ignored (catched by manual lonx/laty bounds check), but also nans are returned on land.
                                                          #assume_sorted=True, #TODO: assume_sorted increases performance?
                                                          )
            elif 'lat' in data_xr_var.coords:
                data_interp_alltimes = data_xr_var.interp(lat=laty,lon=lonx)
            else: #TODO: lat/latitude (should be flexible instead of hardcoded) 
                raise Exception(f'latitude/longitude are not in variable coords: {data_xr_var.coords}. Extend this part of the code for e.g. lat/lon coords')
            data_interp = data_interp_alltimes.sel(time=slice(tstart,tstop))
            time_passed = (dt.datetime.now()-dtstart).total_seconds()
            if debug: print(f'>>time passed: {time_passed:.2f} sec')
            
            # check if only nan (out of bounds or land):
            data_time0 = data_interp_alltimes.isel(time=0).to_numpy()
            if np.isnan(data_time0).all():
                print('WARNING: only nans on first time of this coordinate') #TODO: this can happen on land, raise exception or warning?
            
            if debug:
                print('> plotting')
                dtstart = dt.datetime.now()
                fig,ax = plt.subplots(figsize=(10,7))
                data_interp.T.plot() #uses plot.line() for 1D arrays and plot.pcolormesh() for 2D arrays: https://docs.xarray.dev/en/stable/generated/xarray.DataArray.plot.html
                ax.set_title(f'{quantity} {pli_PolyObject_name_num}')
                if 'depth' in data_xr_var.coords:
                    ax.set_ylim(0,200)
                    ax.invert_yaxis()
                fig.tight_layout()
                time_passed = (dt.datetime.now()-dtstart).total_seconds()
                if debug: print(f'>>time passed: {time_passed:.2f} sec')
            
            print('> converting data to numpy array, ffill nans and concatenating time column')
            dtstart = dt.datetime.now()
            datablock_raw = data_interp.to_numpy()
            
            if has_depth:
                datablock = pd.DataFrame(datablock_raw).fillna(method='ffill',axis=1).values #fill nans forward, is this efficient?
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
            
            timevar_sel = data_interp.time
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
                                    datablock=datablock.tolist(), 
                                    )
            
            ForcingModel_object.forcing.append(ts_one)
            time_passed = (dt.datetime.now()-dtstart).total_seconds()
            if debug: print(f'>>time passed: {time_passed:.2f} sec')
    
    return ForcingModel_object



