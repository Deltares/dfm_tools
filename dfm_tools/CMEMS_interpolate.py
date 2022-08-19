# -*- coding: utf-8 -*-
"""
Created on Thu Aug 18 22:39:03 2022

@author: veenstra
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
from hydrolib.core.io.polyfile.parser import read_polyfile

def get_varnames_dict(dictname='cmems'):
    
    varnameset_dict = {'cmems':
                       {#'':'bottomT', #TODO: update dict
                        'salinity':'so',
                        'temperature':'thetao',
                        #'':'uo',
                        #'':'vo',
                        'steric':'zos',
                        }
                       }
    varnames_dict = varnameset_dict[dictname]
    
    return varnames_dict
    


def interpolate_FES(dir_pattern, file_pli, nPoints=None, debug=False):
    """
    """
    
    nPolyObjects = None
    
    print('initialize ForcingModel()')
    ForcingModel_object = ForcingModel()
    
    file_list_nc = glob.glob(str(dir_pattern))
    component_list = [os.path.basename(x).replace('.nc','') for x in file_list_nc] #TODO: add sorting, manually? Add A0? translate dict for component names or not necessary?
    
    #load boundary file
    polyfile_object = read_polyfile(file_pli,has_z_values=False) #TODO REPORT: this warning can be suppressed (or how to fix): "UserWarning: White space at the start of the line is ignored."
    
    pli_PolyObjects = polyfile_object['objects']
    for iPO, pli_PolyObject_sel in enumerate(pli_PolyObjects[:nPolyObjects]):
        print(f'processing PolyObject {iPO+1} of {len(pli_PolyObjects)}: name={pli_PolyObject_sel.metadata.name}')
        
        for iP, pli_Point_sel in enumerate(pli_PolyObject_sel.points[:nPoints]):
            print(f'processing Point {iP+1} of {len(pli_PolyObject_sel.points)}: ',end='')
            
            lonx, laty = pli_Point_sel.x, pli_Point_sel.y
            #lonx,laty = 13.000000, 54.416667  #extra_rand_dcsm_0001 # is not in lat/lon bounds of selected CMEMS folder? >> probably on land
            #lonx,laty = -9.25, 43.00 #DCSM-FM_OB_all_20181108_0001
            lonx = lonx%360 #for FES since it ranges from 0 to 360 instead of -180 to 180
            print(f'(x={lonx}, y={laty})')
            pli_PolyObject_name_num = f'{pli_PolyObject_sel.metadata.name}_{iP+1:04d}'
            
            datablock_list = []
            for iC,component in enumerate(component_list):
                file_nc = os.path.join(os.path.dirname(dir_pattern),f'{component}.nc')
                data_xr = xr.open_dataset(file_nc)
                lonvar_vals = data_xr['lon'].to_numpy()
                latvar_vals = data_xr['lat'].to_numpy()
                data_xr_amp = data_xr['amplitude']
                data_xr_phs = data_xr['phase']
            
                if iC==0:
                    if (lonx <= lonvar_vals.min()) or (lonx >= lonvar_vals.max()):
                        raise Exception(f'requested lonx {lonx} outside or on lon bounds ({lonvar_vals.min(),lonvar_vals.max()})')
                    if (laty <= latvar_vals.min()) or (laty >= latvar_vals.max()):
                        raise Exception(f'requested laty {laty} outside or on lat bounds ({latvar_vals.min(),latvar_vals.max()})')
            
                data_interp_amp = data_xr_amp.interp(lat=laty,lon=lonx)/100 #convert from cm to m
                data_interp_phs = data_xr_phs.interp(lat=laty,lon=lonx)
                
                # check if only nan (out of bounds or land):
                if np.isnan(data_interp_amp.to_numpy()).all():
                    if iC==0:
                        print('WARNING: only nans for this coordinate') #TODO: this can happen on land, raise exception or warning?
                
                if 0:#debug and component=='m2':
                    print('> plotting')
                    dtstart = dt.datetime.now()
                    fig,(ax1,ax2) = plt.subplots(2,1,figsize=(10,7))
                    data_interp_amp.plot(ax=ax1)
                    data_interp_phs.plot(ax=ax2)
                    ax1.set_title(f'amplitude {component} {pli_PolyObject_name_num}')
                    ax2.set_title(f'phase {component} {pli_PolyObject_name_num}')
                    fig.tight_layout()
                    time_passed = (dt.datetime.now()-dtstart).total_seconds()
                    if debug: print(f'>>time passed: {time_passed:.2f} sec')
                
                datablock_list.append([component.upper(),data_interp_amp.to_numpy(),data_interp_phs.to_numpy()])
            
            # Each .bc file can contain 1 or more timeseries, one for each support point:
            print('> constructing TimeSeries and appending to ForcingModel()')
            """
            name: str = Field(alias="name")
            function: str = Field(alias="function")
            quantityunitpair: List[QuantityUnitPair]

            function: Literal["astronomic"] = "astronomic"
            factor: float = Field(1.0, alias="factor")

            """
            ts_one = Astronomic(name=pli_PolyObject_name_num,
                                quantityunitpair=[QuantityUnitPair(quantity="astronomic component", unit='-'),
                                                  QuantityUnitPair(quantity='waterlevelbnd amplitude', unit='m'),#unit=data_xr_amp.attrs['units']),
                                                  QuantityUnitPair(quantity='waterlevelbnd phase', unit=data_xr_phs.attrs['units'])],
                                datablock=datablock_list, 
                                )
            
            ForcingModel_object.forcing.append(ts_one)        
    return ForcingModel_object

    
def interpolate_nc_to_bc(dir_pattern, file_pli, 
                         modelvarname, varnames_dict, 
                         tstart, tstop, refdate_str, 
                         nPoints=None, debug=False):
    
    """
    nPolyObjects = None #None gives all PolyObjects
    nPoints = 2 #None gives all Points in PolyObject
    """
    nPolyObjects = None
    
    #get bcvarname and varname
    dict_modelvarname2bcvarname = {'salinity':'salinitybnd',
                                   'steric':'waterlevelbnd'} #TODO: is waterlevel not also waterlevelbnd?
    if modelvarname in dict_modelvarname2bcvarname.keys():
        bcvarname = dict_modelvarname2bcvarname[modelvarname]
    else:
        bcvarname = modelvarname
    varname = varnames_dict[modelvarname]
    
    print('initialize ForcingModel()')
    ForcingModel_object = ForcingModel()
    
    file_list_nc = glob.glob(str(dir_pattern))
    print(f'loading mfdataset ({len(file_list_nc)} files with pattern "{dir_pattern.name}")')
    dtstart = dt.datetime.now()
    data_xr = xr.open_mfdataset(file_list_nc)#, combine='by_coords', decode_times=False)
    time_passed = (dt.datetime.now()-dtstart).total_seconds()
    if debug: print(f'>>time passed: {time_passed:.2f} sec')
    timevar = data_xr['time']
    lonvar_vals = data_xr['longitude'].to_numpy()
    latvar_vals = data_xr['latitude'].to_numpy()
    nc_tstart = pd.to_datetime(timevar[0].to_numpy())
    nc_tstop = pd.to_datetime(timevar[-1].to_numpy())
    if tstart < nc_tstart:
        raise Exception(f'requested tstart {tstart} before nc_tstart {nc_tstart}')
    if tstop > nc_tstop:
        raise Exception(f'requested tstop {tstop} after nc_tstop {nc_tstop}')
    
    data_xr_var = data_xr[varname]
    
    if 'depth' in data_xr_var.coords:
        #get variable depths
        vardepth = data_xr_var['depth']
        depth_array = vardepth.to_numpy()[::-1] #fliplr because value array is also flipped #TODO: is this necessary for dflowfm?
        if '_CoordinateZisPositive' in vardepth.attrs.keys(): #correct for positive down to up
            if vardepth.attrs['_CoordinateZisPositive'] == 'down':
                depth_array = -depth_array
    
    #get variable units
    varunit = data_xr_var.attrs['units']
    
    #load boundary file
    #polyfile_object = PolyFile(file_pli,has_z_values=False) #TODO REPORT: (or ask) this does not work, since has_z_values is not passed on or is it not meant to work? (default seems False, so should also work without)
    polyfile_object = read_polyfile(file_pli,has_z_values=False) #TODO REPORT: this warning can be suppressed (or how to fix): "UserWarning: White space at the start of the line is ignored."
    """
    print(len(polyfile_object['objects'])) #1 #gives amount of polyobjects
    print(type(polyfile_object['objects'][0])) # hydrolib.core.io.polyfile.models.PolyObject
    print(polyfile_object['objects'][0]) #gives first polyobject
    print(polyfile_object['objects'][0].metadata) #Metadata(name='extra_rand_dcsm', n_rows=61, n_columns=2)
    print(polyfile_object['objects'][0].points) #gives all points (type==list)
    type(polyfile_object['objects'][0].points[0]) # hydrolib.core.io.polyfile.models.Point 
    lonx_array = [point.x for point in polyfile_object['objects'][0].points]
    lony_array = [point.y for point in polyfile_object['objects'][0].points]
    """
    pli_PolyObjects = polyfile_object['objects']
    for iPO, pli_PolyObject_sel in enumerate(pli_PolyObjects[:nPolyObjects]):
        print(f'processing PolyObject {iPO+1} of {len(pli_PolyObjects)}: name={pli_PolyObject_sel.metadata.name}')
        
        for iP, pli_Point_sel in enumerate(pli_PolyObject_sel.points[:nPoints]):
            print(f'processing Point {iP+1} of {len(pli_PolyObject_sel.points)}: ',end='')
    
            lonx, laty = pli_Point_sel.x, pli_Point_sel.y
            print(f'(x={lonx}, y={laty})')
            #lonx,laty = 13.000000, 54.416667  #extra_rand_dcsm_0001 # is not in lat/lon bounds of selected CMEMS folder? >> probably on land
            #lonx,laty = -9.25, 43.00 #DCSM-FM_OB_all_20181108_0001
            if (lonx <= lonvar_vals.min()) or (lonx >= lonvar_vals.max()):
                raise Exception(f'requested lonx {lonx} outside or on lon bounds ({lonvar_vals.min(),lonvar_vals.max()})')
            if (laty <= latvar_vals.min()) or (laty >= latvar_vals.max()):
                raise Exception(f'requested laty {laty} outside or on lat bounds ({latvar_vals.min(),latvar_vals.max()})')
            
            pli_PolyObject_name_num = f'{pli_PolyObject_sel.metadata.name}_{iP+1:04d}'
            
            print('> interp mfdataset with lat/lon coordinates')
            dtstart = dt.datetime.now()
            data_interp_alltimes = data_xr_var.interp(latitude=laty,longitude=lonx) #TODO: lat/latitude (should be flexible instead of hardcoded) #, kwargs={'bounds_error':True})#, assume_sorted=True) #TODO: bounds_error is ignored, but also nans are returned on land
            data_interp = data_interp_alltimes.sel(time=slice(tstart,tstop))
            time_passed = (dt.datetime.now()-dtstart).total_seconds()
            if debug: print(f'>>time passed: {time_passed:.2f} sec')
            
            # check if only nan (out of bounds or land):
            data_time0 = data_interp_alltimes.sel(time=data_xr['time'][1]).to_numpy()
            if np.isnan(data_time0).all():
                raise Exception('only nans on first time of this coordinate') #TODO: this can happen on land, raise exception or warning?
            
            if debug:
                print('> plotting')
                dtstart = dt.datetime.now()
                fig,ax = plt.subplots(figsize=(10,7))
                data_interp.T.plot()
                ax.set_title(f'{modelvarname} {pli_PolyObject_name_num}')
                if 'depth' in data_xr_var.coords:
                    ax.set_ylim(0,200)
                    ax.invert_yaxis()
                fig.tight_layout()
                time_passed = (dt.datetime.now()-dtstart).total_seconds()
                if debug: print(f'>>time passed: {time_passed:.2f} sec')
            
            print('> converting data to numpy array, ffill nans and concatenating time column')
            dtstart = dt.datetime.now()
            datablock_raw = data_interp.to_numpy()
            if 'depth' in data_xr_var.coords:
                datablock = pd.DataFrame(datablock_raw).fillna(method='ffill',axis=1).values #fill nans forward, is this efficient?
                #if debug: print(datablock_raw), print(datablock)
                datablock = datablock[:,::-1]#.flip(axis=1) #flipping axis #TODO: this assumes depth as second dimension, might not be true (maybe avoidable anyway)
            else:
                datablock = datablock_raw[:,np.newaxis]
                
            timevar_sel = data_interp.time
            timevar_sel_rel = date2num(pd.DatetimeIndex(timevar_sel.to_numpy()).to_pydatetime(),units=refdate_str,calendar='standard')
            datablock_incltime = np.concatenate([timevar_sel_rel[:,np.newaxis],datablock],axis=1)
            datablock_list = datablock_incltime.tolist()
            time_passed = (dt.datetime.now()-dtstart).total_seconds()
            if debug: print(f'>>time passed: {time_passed:.2f} sec')
            
            # Each .bc file can contain 1 or more timeseries, one for each support point:
            print('> constructing TimeSeries and appending to ForcingModel()')
            dtstart = dt.datetime.now()
            if 'depth' in data_xr_var.coords:
                """
                name: str = Field(alias="name")
                function: str = Field(alias="function")
                quantityunitpair: List[QuantityUnitPair]

                offset: float = Field(0.0, alias="offset")
                factor: float = Field(1.0, alias="factor")
                verticalpositions: List[float] = Field(alias="verticalPositions")
                verticalinterpolation: VerticalInterpolation = Field(alias="verticalInterpolation")
                verticalpositiontype: VerticalPositionType = Field(alias="verticalPositionType")
                """
                list_QUP_perlayer = [QuantityUnitPair(quantity=bcvarname, unit=varunit) for iL in range(vardepth.size)] #TODO REPORT: verticalposition 1/2/3/n is not supported
                ts_one = T3D(name=pli_PolyObject_name_num,
                             verticalpositions=depth_array.tolist(), #TODO REPORT: should be "Vertical position specification = [..]" but is verticalPositions = [..]"
                             verticalInterpolation='linear', #TODO REPORT: how about extrapolation? Also: this is not necessary in bc file since there is a default value specified in dflowfm so should not be required
                             verticalPositionType=VerticalPositionType('ZBed'), #TODO REPORT: should be "Vertical position type = zdatum" but is "verticalPositionType = ZBed"
                             quantityunitpair=[QuantityUnitPair(quantity="time", unit=refdate_str)]+list_QUP_perlayer,
                             timeinterpolation=TimeInterpolation.linear, #TODO REPORT: not passed on to bc file
                             datablock=datablock_list, 
                             )
            else:
                """
                name: str = Field(alias="name")
                function: str = Field(alias="function")
                quantityunitpair: List[QuantityUnitPair]

                timeinterpolation: TimeInterpolation = Field(alias="timeInterpolation")
                offset: float = Field(0.0, alias="offset")
                factor: float = Field(1.0, alias="factor")
                """
                ts_one = TimeSeries(name=pli_PolyObject_name_num,
                                    verticalposition=VerticalPositionType('ZBed'), #TODO REPORT: is not passed on to bc file, so should raise error
                                    quantityunitpair=[QuantityUnitPair(quantity="time", unit=refdate_str),
                                                      QuantityUnitPair(quantity=bcvarname, unit=varunit)],
                                    timeinterpolation=TimeInterpolation.linear,
                                    datablock=datablock_list, 
                                    )
            
            ForcingModel_object.forcing.append(ts_one)
            time_passed = (dt.datetime.now()-dtstart).total_seconds()
            if debug: print(f'>>time passed: {time_passed:.2f} sec')
    
    return ForcingModel_object