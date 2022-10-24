# -*- coding: utf-8 -*-
"""
Created on Thu Aug 18 16:37:38 2022

@author: veenstra

This script can be used to interpolate pre-downloaded CMEMS data (or other netcdf files) to a boundary locations in a pli-file
"""

import datetime as dt
from pathlib import Path
import matplotlib.pyplot as plt
plt.close('all')
from dfm_tools.interpolate_grid2bnd import get_conversion_dict, interpolate_tide_to_bc, interpolate_nc_to_bc
from dfm_tools.hydrolib_helpers import forcinglike_to_Dataset
from hydrolib.core.io.ext.models import Boundary, ExtModel

#TODO: add coordinate conversion of pli-coordinates (for nesting RD models)

model = 'CMEMS' #CMEMS GFDL CMCC HYCOM

if model == 'HYCOM': #not available in dcsm area
    list_plifiles = [Path(r'c:\DATA\dfm_tools_testdata\GLBu0.08_expt_91.2\bcline.pli')]
else:
    #copied plifile from DCSM folder: r'p:\1204257-dcsmzuno\data\CMEMS\bnd\NorthSeaAndBaltic_1993-2019_20210510'
    #list_plifiles = [Path(r'n:\My Documents\werkmap\hydrolib_test\DCSM\DCSM-FM_OB_all_20181108.pli')] #TODO: reading this file results in empty Polyfile, should raise an error. https://github.com/Deltares/HYDROLIB-core/issues/320
    list_plifiles = [Path(r'n:\My Documents\werkmap\hydrolib_test\DCSM\DCSM-FM_OB_all_20181108_nocomments.pli')]

dir_out = r'n:\My Documents\werkmap\hydrolib_test\DCSM'
bc_type = 'bc' #currently only 'bc' supported #TODO: add netcdf bc support. https://github.com/Deltares/HYDROLIB-core/issues/318

refdate_str = 'minutes since 2011-12-22 00:00:00 +00:00' # if None, xarray uses ds.time.encoding['units'] as refdate_str

nPoints = 3 #amount of Points to process per PolyObject in the plifile (for testing, use None for all Points)

#quantities should be in conversion_dict.keys(). waterlevelbnd is steric/zos, tide is tidal components from FES/EOT
list_quantities = ['waterlevelbnd','salinitybnd','tide','uxuy','temperaturebnd','tracerbndNO3']
list_quantities = ['uxuy']#,'temperaturebnd','waterlevelbnd','tide']

dtstart = dt.datetime.now()
ext_bnd = ExtModel()



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


def interpolate_nc_to_bc(dir_pattern, file_pli, quantity, 
                         tstart, tstop, refdate_str=None, 
                         conversion_dict=None, #rename_vars={}, #TODO: alternatively use rename_vars dict and use conversion_dict only for unit conversion. dict containing keys: ncvarname, bcvarname and optionally conversion and unit
                         nPoints=None, #argument for testing
                         reverse_depth=False, #temporary argument to compare easier with old coastserv files
                         ):
    
    if conversion_dict is None:
        conversion_dict = get_conversion_dict()
    
    if quantity=='uxuy': #T3Dvector
        print(f'combined variables ({quantity})')
        vector = True
        quantity_list = ['ux','uy']
        ncvarname_list = [conversion_dict[quan]['ncvarname'] for quan in quantity_list]
        ncvarname_joined = ','.join(ncvarname_list)
    else:
        vector = False
        quantity_list = [quantity]
    
    print('initialize ForcingModel()')
    ForcingModel_object = ForcingModel()
    
    if vector:
        dir_pattern = [Path(str(dir_pattern).format(ncvarname=ncvarname)) for ncvarname in ncvarname_list]
        file_list_nc = glob.glob(str(dir_pattern[0])) + glob.glob(str(dir_pattern[1]))
        print(f'loading mfdataset ({len(file_list_nc)} files with pattern "{dir_pattern[0].name}" and "{dir_pattern[1].name}")')
    else:
        file_list_nc = glob.glob(str(dir_pattern).format(ncvarname))
        print(f'loading mfdataset ({len(file_list_nc)} files with pattern "{dir_pattern.name}")')
    
    dtstart = dt.datetime.now()
    data_xr = xr.open_mfdataset(file_list_nc,chunks={'time':1}) # #TODO: does chunks argument solve "PerformanceWarning: Slicing is producing a large chunk."?
    
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
            data_xr = data_xr.rename({'lon':'longitude','lat':'latitude'})
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
            #TODO: kdtree k=3 and invdist weighing? Then also spherical coordinate distance calculation instead of cartesian/eucledian
            print(f'ValueError: {e}. Reverting to KDTree instead (nearest neigbour)')
            from scipy.spatial import KDTree #TODO: move up
            path_lonlat_pd = pd.DataFrame({'lon':da_lons,'lat':da_lats})
            data_lon_flat = data_xr_var[coordname_lon].to_numpy().ravel()
            data_lat_flat = data_xr_var[coordname_lat].to_numpy().ravel()
            data_lonlat_pd = pd.DataFrame({'lon':data_lon_flat,'lat':data_lat_flat})
            #KDTree, finds minimal eucledian distance between points (maybe haversine would be better)
            tree = KDTree(data_lonlat_pd) #alternatively sklearn.neighbors.BallTree: tree = BallTree(data_lonlat_pd)
            distance, data_lonlat_idx = tree.query(path_lonlat_pd, k=1) #TODO: maybe add outofbounds treshold for distance
            #data_lonlat_pd.iloc[data_lonlat_idx]
            idx_i,idx_j = np.divmod(data_lonlat_idx, data_xr_var['longitude'].shape[1]) #get idx i and j by sort of counting over 2D array
            # fig,ax = plt.subplots()
            # data_xr_var_tsel_dsel = data_xr_var.isel(time=0,depth=0)
            # data_xr_var_tsel_dsel.plot(ax=ax)#,x='longitude',y='latitude')
            # ax.plot(idx_j,idx_i,'xr')
            da_idxi = xr.DataArray(idx_i, dims='latloncombi')
            da_idxj = xr.DataArray(idx_j, dims='latloncombi')
            data_interp = data_xr_var.isel(i=da_idxi,j=da_idxj)

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





for file_pli in list_plifiles:
    for quantity in list_quantities:
        
        #TODO: take out of loop (but quantity is required)
        if model=='CMEMS': #2012-01-06 12:00:00 to 2013-01-03 12:00:00
            conversion_dict = get_conversion_dict()
            tstart = dt.datetime(2012, 1, 16, 12, 0)
            tstop = dt.datetime(2012, 4, 1, 12, 0)
            dir_sourcefiles_hydro = r'p:\1204257-dcsmzuno\data\CMEMS\nc\DCSM_allAvailableTimes' #CMEMS hydro: bottomT, so, thetao, uo, vo, zos (2012-01-06 12:00:00 to 2013-01-03 12:00:00) (daily values at noon, not at midnight)
            #dir_pattern_hydro = Path(dir_sourcefiles_hydro,f'{conversion_dict[quantity]["ncvarname"]}_2012*.nc') # later remove 2012 from string, but this is faster for testing #TODO: it is quite slow, maybe speed up possible?
            dir_sourcefiles_waq = r'p:\11206304-futuremares\python_scripts\ocean_boundaryCMEMS\data_monthly' #CMEMS waq: no3, o2, phyc, so4, si (2011-12-16 12:00:00 to 2019-01-16 12:00:00)
            #dir_pattern_waq = Path(dir_sourcefiles_waq,f'cmems_mod_glo_bgc_my_0.25_P1M-m_{conversion_dict[quantity]["ncvarname"]}_*.nc') 
            dir_pattern_hydro = Path(dir_sourcefiles_hydro,'{ncvarname}_2012*.nc') # later remove 2012 from string, but this is faster for testing #TODO: it is quite slow, maybe speed up possible?
            dir_pattern_waq = Path(dir_sourcefiles_waq,'cmems_mod_glo_bgc_my_0.25_P1M-m_{ncvarname}_*.nc') 
        elif model=='GFDL':
            conversion_dict = get_conversion_dict()
            tstart = dt.datetime(2012, 1, 16, 12, 0)
            tstop = dt.datetime(2012, 4, 1, 12, 0)
            dir_sourcefiles_hydro = None
            dir_pattern_hydro = None
            dir_sourcefiles_waq = r'p:\11206304-futuremares\data\CMIP6_BC\GFDL-ESM4' #GFDL waq: no3 (1850-01-16 12:00:00 to 2014-12-16 12:00:00)
            dir_pattern_waq = Path(dir_sourcefiles_waq,f'{conversion_dict[quantity]["ncvarname"]}_esm-hist.nc')
        elif model=='CMCC': #TODO: check method, now finding nearest points (so always has values)
            conversion_dict = get_conversion_dict(ncvarname_updates={'salinitybnd':'sos', 'temperaturebnd':'tos'})
            tstart = dt.datetime(2015, 6, 16, 12, 0)
            tstop = dt.datetime(2016, 12, 1, 12, 0)
            dir_sourcefiles_hydro = r'p:\11206304-futuremares\data\CMIP6_BC\CMCC-ESM2'
            dir_pattern_hydro = Path(dir_sourcefiles_hydro,f'{conversion_dict[quantity]["ncvarname"]}_Omon_CMCC-ESM2_ssp126_r1i1p1f1_gn_*.nc')
            dir_sourcefiles_waq = dir_sourcefiles_hydro #CMCC waq: (2015-01-16 12:00:00 to 2100-12-16 12:00:00)
            dir_pattern_waq = dir_pattern_hydro
        elif model=='HYCOM':
            if quantity not in ['salinitybnd','temperaturebnd']: #only contains two quantities
                continue
            conversion_dict = get_conversion_dict(ncvarname_updates={'salinitybnd':'salinity', 'temperaturebnd':'water_temp'})
            tstart = dt.datetime(2016, 4, 20, 0, 0) #HYCOM
            tstop = dt.datetime(2016, 5, 3, 0, 0)
            dir_sourcefiles_hydro = 'c:\\DATA\\dfm_tools_testdata\\GLBu0.08_expt_91.2' #HYCOM hydro: salinity/so, water_temp/thetao (2016-04-19 00:00:00 to 2016-05-06 00:00:00)
            dir_pattern_hydro = Path(dir_sourcefiles_hydro,'HYCOM_ST_GoO_*.nc')
            dir_sourcefiles_waq = None
            dir_pattern_waq = None
        else:
            raise Exception(f'invalid model: {model}')
        
        print(f'processing quantity: {quantity}')#/{conversion_dict[quantity]["ncvarname"]}')
        if quantity in ['tide']: #tide #TODO: choose flexible/generic component notation
            tidemodel = 'FES2014' #FES2014, FES2012, EOT20
            component_list = ['2n2','mf','p1','m2','mks2','mu2','q1','t2','j1','m3','mm','n2','r2','k1','m4','mn4','s1','k2','m6','ms4','nu2','s2','l2','m8','msf','o1','s4'] #None results in all FES components
            ForcingModel_object = interpolate_tide_to_bc(tidemodel=tidemodel, file_pli=file_pli, component_list=component_list, nPoints=nPoints)
            for forcingobject in ForcingModel_object.forcing: #add A0 component
                forcingobject.datablock.append(['A0',0.0,0.0])
        elif quantity in ['waterlevelbnd','salinitybnd','temperaturebnd','uxuy']: #hydro
            if dir_sourcefiles_hydro is None:
                continue
            ForcingModel_object = interpolate_nc_to_bc(dir_pattern=dir_pattern_hydro, file_pli=file_pli,
                                                       quantity=quantity, conversion_dict=conversion_dict,
                                                       tstart=tstart, tstop=tstop, refdate_str=refdate_str,
                                                       #reverse_depth=True, #to compare with coastserv files, this argument will be phased out
                                                       nPoints=nPoints)
        else: #waq
            if dir_pattern_waq is None:
                continue
            ForcingModel_object = interpolate_nc_to_bc(dir_pattern=dir_pattern_waq, file_pli=file_pli,
                                                       quantity=quantity, conversion_dict=conversion_dict,
                                                       tstart=tstart, tstop=tstop, refdate_str=refdate_str,
                                                       #reverse_depth=True, #to compare with coastserv files, this argument will be phased out
                                                       nPoints=nPoints)
        
        if 1: #plotting example data point
            for iF in [2]:#range(nPoints):
                forcingobject_one = ForcingModel_object.forcing[iF]
                forcingobject_one_xr = forcinglike_to_Dataset(forcingobject_one,convertnan=True)
                data_vars = list(forcingobject_one_xr.data_vars)
                fig,ax1 = plt.subplots(figsize=(10, 6))
                if hasattr(forcingobject_one,'vertpositions'): #3D quantity (time/depth dimensions)
                    if hasattr(forcingobject_one.quantityunitpair[1],'elementname'): #uxuy vector
                        plt.close()
                        fig, axes = plt.subplots(2,1,figsize=(10, 6),sharex=True,sharey=True)
                        forcingobject_one_xr[data_vars[0]].T.plot(ax=axes[0])
                        forcingobject_one_xr[data_vars[1]].T.plot(ax=axes[1])
                    else:
                        forcingobject_one_xr[data_vars[0]].T.plot(ax=ax1)
                elif quantity=='tide':
                    ax2 = ax1.twinx()
                    forcingobject_one_xr[data_vars[0]].plot(ax=ax1)
                    forcingobject_one_xr[data_vars[1]].plot(ax=ax2)
                else:
                    forcingobject_one_xr[data_vars[0]].plot(ax=ax1)
        
        file_bc_basename = file_pli.name.replace('.pli','')
        if quantity=='tide':
            file_bc_out = Path(dir_out,f'{quantity}_{file_bc_basename}_{tidemodel}.bc')
        else:
            file_bc_out = Path(dir_out,f'{quantity}_{file_bc_basename}_{model}.bc')
        print(f'writing ForcingModel to bc file with hydrolib ({file_bc_out.name})')
        if bc_type=='bc':
            ForcingModel_object.save(filepath=file_bc_out)
            #TODO: improve formatting of bc file to make nicer, save diskspace and maybe write faster: https://github.com/Deltares/HYDROLIB-core/issues/308 (and https://github.com/Deltares/HYDROLIB-core/issues/313)
        else:
            raise Exception(f'invalid bc_type: {bc_type}')
        
        #make paths relative (sort of) (also necessary for locationfile) /../ should also be supported? 
        #ForcingModel_object.filepath = Path(str(ForcingModel_object.filepath).replace(dir_out,'')) #TODO: convert to relative paths in ext file possible? This path is the same as file_bc_out
        
        #generate boundary object for the ext file (quantity, pli-filename, bc-filename)
        boundary_object = Boundary(quantity=quantity,
                                   locationfile=Path(dir_out,file_pli.name),
                                   forcingfile=ForcingModel_object,
                                   )
        ext_bnd.boundary.append(boundary_object)

file_ext_out = Path(dir_out,'example_bnd.ext')
ext_bnd.save(filepath=file_ext_out)

time_passed = (dt.datetime.now()-dtstart).total_seconds()
print(f'>>time passed: {time_passed:.2f} sec')

