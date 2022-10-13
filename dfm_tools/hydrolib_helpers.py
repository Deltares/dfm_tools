# -*- coding: utf-8 -*-
"""
Created on Tue Aug 23 13:36:44 2022

@author: veenstra
"""

import datetime as dt
import pandas as pd
import cftime
import numpy as np
import xarray as xr
from hydrolib.core.io.polyfile.models import PolyObject
#from netCDF4 import date2num #TODO: take from cftime?
from cftime import date2num

from hydrolib.core.io.bc.models import (
    ForcingModel,
    QuantityUnitPair,
    T3D,
    TimeSeries,
    Astronomic,
)

#TODO: maybe add DataFrame_to_forcingobject() and others


def DataArray_to_T3D(datablock_xr, 
                     name, #TODO: add name to DataArray attrs? (also add refdate?)
                     refdate_str, bcvarname, fill_na=True, 
                     depthvarname='depth'): #TODO depthvarname argument can be avoided if rename_variables with 'depth'
    """
    convert an xarray.DataArray with depth dimension to a hydrolib T3D object
    """
    
    #get depth variable and values
    depth_array = datablock_xr[depthvarname].to_numpy()
    if datablock_xr[depthvarname].attrs['positive'] == 'down': #attribute appears in CMEMS, GFDL and CMCC, save to assume presence?
        depth_array = -depth_array
    if fill_na: #ffill data
        datablock_xr = datablock_xr.bfill(dim=depthvarname).ffill(dim=depthvarname) #fill nans back and forward (corresponds to vertical extrapolation for CMEMS). Deep values are filled if order is shallow to deep
    
    #get datablock and concatenate with relative time data
    datablock_np = datablock_xr.to_numpy()
    timevar_sel_rel = date2num(pd.DatetimeIndex(datablock_xr.time.to_numpy()).to_pydatetime(),units=refdate_str,calendar='standard')
    datablock_incltime = np.concatenate([timevar_sel_rel[:,np.newaxis],datablock_np],axis=1)
    
    # Each .bc file can contain 1 or more timeseries, in this case one for each support point
    verticalpositions_idx = np.arange(datablock_xr[depthvarname].size)+1
    #list_QUP_perlayer = [QuantityUnitPair(quantity=bcvarname, unit=datablock_xr.attrs['units']) for iVP in verticalpositions_idx] #TODO SOLVED: verwarrende foutmelding bij niet opgeven verticalpositionindex (should be missing error instead of not valid error)
    list_QUP_perlayer = [QuantityUnitPair(quantity=bcvarname, unit=datablock_xr.attrs['units'], vertpositionindex=iVP) for iVP in verticalpositions_idx] #TODO SOLVED: verticalposition 1/2/3/n is not supported. https://github.com/Deltares/HYDROLIB-core/issues/317
    #TODO: instead of supplying list of QuantityUnitPairs, it is also possible to supply the three quanties as list separately
    ts_one = T3D(name=name,
                 vertpositions=depth_array.tolist(), #TODO SOLVED: should be "Vertical position specification = [..]" but is verticalPositions = [..]" (both possible?). https://github.com/Deltares/HYDROLIB-core/issues/317
                 vertinterpolation='linear', #TODO SOLVED: not providing this results in VerticalInterpolation.linear
                 vertPositionType='ZDatum', #TODO SOLVED: should be "Vertical position type = zdatum" but is "verticalPositionType = ZBed" (zdatum is niet beschikbaar). https://github.com/Deltares/HYDROLIB-core/issues/317
                 quantityunitpair=[QuantityUnitPair(quantity="time", unit=refdate_str)]+list_QUP_perlayer,
                 timeinterpolation='linear', #TODO SOLVED: not passed on to bc file. https://github.com/Deltares/HYDROLIB-core/issues/317
                 datablock=datablock_incltime.tolist(), #TODO: numpy array is not supported by TimeSeries. https://github.com/Deltares/HYDROLIB-core/issues/322
                 )
    return ts_one


def DataArray_to_TimeSeries(datablock_xr, name, refdate_str, bcvarname):
    """
    convert an xarray.DataArray without depth dimension to a hydrolib TimeSeries object
    """
    
    #get datablock and concatenate with relative time data
    datablock_np = datablock_xr.to_numpy()[:,np.newaxis]
    timevar_sel_rel = date2num(pd.DatetimeIndex(datablock_xr.time.to_numpy()).to_pydatetime(),units=refdate_str,calendar='standard')
    datablock_incltime = np.concatenate([timevar_sel_rel[:,np.newaxis],datablock_np],axis=1)
    
    # Each .bc file can contain 1 or more timeseries, in this case one for each support point
    ts_one = TimeSeries(name=name,
                        quantityunitpair=[QuantityUnitPair(quantity="time", unit=refdate_str), #TODO: quantity is not validated: https://github.com/Deltares/HYDROLIB-core/issues/357
                                          QuantityUnitPair(quantity=bcvarname, unit=datablock_xr.attrs['units'])],
                        timeinterpolation='linear',
                        datablock=datablock_incltime.tolist(), 
                        )
    return ts_one


def DataFrame_to_PolyObject(poly_pd,name,content=None): #TODO: make this method bound? Better: make polyobject/extmodel accept dataframe?
    """
    convert a pandas dataframe with x/y columns (and optional others like z/data/comment) to a hydrolib PolyObject
    """
    if 'z' in poly_pd.columns:
        nondata_cols = ['x','y','z']
    else:
        nondata_cols = ['x','y']
    poly_pd_xy = poly_pd[nondata_cols] #TODO: actually z is also a thing, but that becomes part of data in this method
    poly_pd_data = pd.DataFrame({'data':poly_pd.drop(nondata_cols,axis=1).values.tolist()})
    poly_pd_polyobj = pd.concat([poly_pd_xy,poly_pd_data],axis=1)
    pointsobj_list = poly_pd_polyobj.T.apply(dict).tolist() #TODO: maybe faster with list iteration
    polyobject = PolyObject(metadata={'name':name,'n_rows':poly_pd.shape[0],'n_columns':poly_pd.shape[1]}, points=pointsobj_list)
    if content is not None:
        polyobject.description = {'content':content}
    return polyobject


def forcinglike_to_DataFrame(forcingobj, convert_time=True): #TODO: remove this method?
    """
    convert a hydrolib forcing like object (like Timeseries, T3D, Harmonic, etc) to a pandas DataFrame. #TODO: astronomic is also supported, what more?
    
    Parameters
    ----------
    forcingobj : hydrolib ForcingModel.forcing object (Timeseries/T3D etc)
        DESCRIPTION.
    convert_time : boolean, optional
        Convert time column from unit (e.g. minutes since date) to datetime index and drop the time column. Has no effect if there is no time column in the forcingobject. The default is True.

    Returns
    -------
    df_data : pd.DataFrame
        DESCRIPTION
        
    Example
    -------
         file_bc = Path('p:\\11208053-004-kpp2022-rmm1d2d\\C_Work\\09_Validatie2018_2020\\dflowfm2d-rmm_vzm-j19_6-v2d\\boundary_conditions\\2018\\flow\\rmm_discharge_laterals_20171220_20190101_MET.bc')
         m = ForcingModel(file_bc)
         df_data_list = [forcingobject_to_dataframe(forcingobj, convert_time=True) for forcingobj in m.forcing]

    """
    raise Exception('do not use this method')
    
    #check if forcingmodel instead of T3D/TimeSeries is provided
    if isinstance(forcingobj, ForcingModel):
        raise Exception('ERROR: instead of supplying a ForcingModel, provide a ForcingObject (Timeseries/T3D etc), by doing something like ForcingModel.forcing[0]')
    
    allowed_instances = (T3D, TimeSeries, Astronomic)
    if not isinstance(forcingobj, allowed_instances):
        raise Exception(f'ERROR: supplied input is not one of: {allowed_instances}')
    
    #if hasattr(forcingobj.quantityunitpair[0],'verticalpositionindex'): #TODO: might be there always
    QUP_list = [(QUP.quantity,QUP.unit,QUP.vertpositionindex) for QUP in forcingobj.quantityunitpair] #TODO: generating MultiIndex can probably be more elegant (e.g. getting names from QUP list), but I do not know how
    columns_MI = pd.MultiIndex.from_tuples(QUP_list,names=dict(forcingobj.quantityunitpair[0]).keys())
    df_data = pd.DataFrame(forcingobj.datablock,columns=columns_MI)
    df_data.index.name = forcingobj.name
    colnames_quantity = df_data.columns.get_level_values(level=0)
    if convert_time and ('time' in colnames_quantity): #this converts time to a datetime index #TODO: do automatically if TimeSeries/T3D? (save 'encoding'/refdate somewhere)
        time_colid = colnames_quantity.get_loc('time')
        time_unit = df_data.columns.get_level_values(level=1)[time_colid]
        df_data.index = cftime.num2pydate(df_data.iloc[:,time_colid],units=time_unit)
        df_data.index.name = forcingobj.name #again with new index
        #timezone was converted to GMT, re-adjust timezone if needed
        timeunit_sincedatetimetz = time_unit.split('since ')[1]
        tzone_minutes = cftime._parse_date(timeunit_sincedatetimetz)[-1]
        df_data.index = df_data.index.tz_localize('GMT')
        df_data.index = df_data.index.tz_convert(dt.timezone(dt.timedelta(minutes=tzone_minutes)))
        #drop original time column
        df_data = df_data.drop(labels='time',level=0,axis=1)
    return df_data


def forcinglike_to_DataArray(forcingobj): #TODO: would be convenient to have this as a method of ForcingModel objects (Timeseries/T3D/etc), or maybe as method of the ForcingModel (returning a list of DataFrames): https://github.com/Deltares/HYDROLIB-core/issues/307
    """
    convert a hydrolib forcing like object (like Timeseries, T3D, Harmonic, etc) to a xarray DataArray.
    #TODO: add doc
    """
    
    #check if forcingmodel instead of T3D/TimeSeries is provided
    if isinstance(forcingobj, ForcingModel):
        raise Exception('ERROR: instead of supplying a ForcingModel, provide a ForcingObject (Timeseries/T3D etc), by doing something like ForcingModel.forcing[0]')
    
    allowed_instances = (T3D, TimeSeries, Astronomic)
    if not isinstance(forcingobj, allowed_instances):
        raise Exception(f'ERROR: supplied input is not one of: {allowed_instances}')
    
    var_quantity = forcingobj.quantityunitpair[1].quantity
    var_unit = forcingobj.quantityunitpair[1].unit
    if isinstance(forcingobj, T3D):
        dims = ('time','depth')
    elif isinstance(forcingobj, TimeSeries):
        dims = ('time')
    elif isinstance(forcingobj, Astronomic):
        dims = ('astronomic_component','quantity') #TODO: what are dims in case of Astronomic? actually these are two variables
        print('WARNING: format of DataArray is not final in case of astronomic component')
    
    datablock_all = np.array(forcingobj.datablock)
    if forcingobj.quantityunitpair[0].quantity in ['time','astronomic component']: #first column in bcfile is time
        datablock_data = datablock_all[:,1:] #TODO: convert repeating values to nan? (reverse of ffill/bfill)
        datablock_data = datablock_data.squeeze() #drop dimensions of len 1 in case of eg "waterlevelbnd"
        datablock_data = datablock_data.astype(float) #convert str to float in case of "astronomic component"
    else:
        datablock_data = datablock_all
    
    #add dimension values
    data_xr_var = xr.DataArray(datablock_data,name=var_quantity,dims=dims)#,coords=coords)
    if 'depth' in dims:
        data_xr_var['depth'] = forcingobj.vertpositions
    if 'time' in dims:
        time_unit = forcingobj.quantityunitpair[0].unit.lower() #TODO: save this as encoding/variable attribute (align with DataArray_to_*)
        data_xr_var['time'] = cftime.num2pydate(datablock_all[:,0], units=time_unit)
    if 'astronomic_component' in dims:
        data_xr_var['astronomic_component'] = datablock_all[:,0]
    
    #add attributes
    attr_dict = {'source':'hydrolib-core object converted to xarray.DataArray with dfm_tools',
                 'unit':var_unit,
                 }
    for key in attr_dict.keys():
        data_xr_var.attrs[key] = attr_dict[key]
    forcingobj_keys = forcingobj.__dict__.keys()
    for key in forcingobj_keys: #['comments','name','function','offset','factor','vertinterpolation','vertpositiontype','timeinterpolation']: 
        if key in ['datablock','quantityunitpair','vertpositions']: #skipping these since they are in the DataArray already
            continue
        data_xr_var.attrs[key] = forcingobj.__dict__[key]
    
    return data_xr_var


def pointlike_to_DataFrame(pointlike,drop_emptycols=True):
    """
    convert a hydrolib object with points (like PolyObject, XYZModel and possibly others) to a pandas DataFrame.

    Parameters
    ----------
    PolyObject : TYPE
        PolyObject, or actually any object.

    Returns
    -------
    poly_pd : TYPE
        DESCRIPTION.
        
    Example:
        polyfile_object = read_polyfile(file_ldb,has_z_values=False)
        poly_pd_list = [polyobject_to_dataframe(PO) for PO in polyfile_object['objects']]

    """
    """
    #x,y,z = PolyObject.xyz_coordinates #TODO: this might be possible in new hydrolib-core version: https://github.com/Deltares/HYDROLIB-core/issues/329
    xvals_pd = pd.DataFrame({'x':[p.x for p in PolyObject.points]})
    yvals_pd = pd.DataFrame({'y':[p.y for p in PolyObject.points]})
    zvals_pd = pd.DataFrame({'z':[p.z for p in PolyObject.points]})
    datavals_pd = pd.DataFrame([p.data for p in PolyObject.points])
    if zvals_pd['z'].isnull().all(): #ignore column if all values are None
        polyobject_pd = pd.concat([xvals_pd,yvals_pd,datavals_pd],axis=1)
    else:
        polyobject_pd = pd.concat([xvals_pd,yvals_pd,zvals_pd,datavals_pd],axis=1)
    """
    pointlike_pd = pd.DataFrame([dict(p) for p in pointlike.points])
    if 'data' in pointlike_pd.columns:
        #datavals_pd = pointlike_pd['data'].apply(pd.Series) #this is quite slow, so use line below instead. maybe lambda or faster approach?
        datavals_pd = pd.DataFrame([p.data for p in pointlike.points])
        pointlike_pd = pd.concat([pointlike_pd.drop(['data'],axis=1), datavals_pd],axis=1)
        
    if drop_emptycols:
        allempty_bool = pointlike_pd.isnull().all(axis=0)
        pointlike_pd = pointlike_pd.loc[:,~allempty_bool]
        
    return pointlike_pd


def parse_xy_to_datetime(pointlike_pd):
    datatimevals_pdstr = (pointlike_pd['x'].astype(int).apply(lambda x:f'{x:08d}') +
                          pointlike_pd['y'].astype(int).apply(lambda x:f'{x:06d}'))
    pointlike_pd.index = pd.to_datetime(datatimevals_pdstr)
    pointlike_pd_timeidx = pointlike_pd.drop(['x','y'],axis=1)
    
    return pointlike_pd_timeidx



