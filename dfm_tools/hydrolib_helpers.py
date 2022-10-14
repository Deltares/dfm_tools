# -*- coding: utf-8 -*-
"""
Created on Tue Aug 23 13:36:44 2022

@author: veenstra
"""

import pandas as pd
import cftime
import numpy as np
import xarray as xr
from hydrolib.core.io.polyfile.models import PolyObject
from cftime import date2num

from hydrolib.core.io.bc.models import (
    ForcingModel,
    QuantityUnitPair,
    T3D,
    TimeSeries,
    Astronomic,
)


def DataArray_to_T3D(datablock_xr):   
    """
    convert an xarray.DataArray with time and depth dimension to a hydrolib T3D object
    """
    #TODO: clean up these first lines of code and add description to docstring?
    locationname = datablock_xr.attrs['locationname']
    bcvarname = datablock_xr.name
    refdate_str = datablock_xr.time.encoding['units']
    
    if datablock_xr.dims != ('time','depth'): #check if both time and depth dimensions are present #TODO: add support for flipped dimensions (datablock_xr.T or something is needed)
        raise Exception(f"ERROR: datablock_xr provided to DataArray_to_T3D has dimensions {datablock_xr.dims} while ('time','depth') is expected")
    
    #get depth variable and values
    depth_array = datablock_xr['depth'].to_numpy()
    if 'positive' in datablock_xr['depth'].attrs.keys():
        if datablock_xr['depth'].attrs['positive'] == 'down': #attribute appears in CMEMS, GFDL and CMCC, save to assume presence?
            depth_array = -depth_array
    
    #ffill/bfill nan data along over depth dimension (corresponds to vertical extrapolation)
    datablock_xr = datablock_xr.bfill(dim='depth').ffill(dim='depth')
    
    #get datablock and concatenate with relative time data
    datablock_np = datablock_xr.to_numpy()
    timevar_sel_rel = date2num(pd.DatetimeIndex(datablock_xr.time.to_numpy()).to_pydatetime(),units=refdate_str,calendar='standard')
    datablock_incltime = np.concatenate([timevar_sel_rel[:,np.newaxis],datablock_np],axis=1)
    
    # Each .bc file can contain 1 or more timeseries, in this case one for each support point
    verticalpositions_idx = np.arange(datablock_xr['depth'].size)+1
    list_QUP_perlayer = [QuantityUnitPair(quantity=bcvarname, unit=datablock_xr.attrs['units'], vertpositionindex=iVP) for iVP in verticalpositions_idx]
    ts_one = T3D(name=locationname,
                 vertpositions=np.round(depth_array.tolist(),decimals=4).tolist(), # make decimals userdefined? .tolist() is necessary for np.round to work for some reason
                 vertinterpolation='linear',
                 vertPositionType='ZDatum',
                 quantityunitpair=[QuantityUnitPair(quantity="time", unit=refdate_str)]+list_QUP_perlayer,
                 timeinterpolation='linear',
                 datablock=datablock_incltime.tolist(), #TODO: numpy array is not supported by TimeSeries. https://github.com/Deltares/HYDROLIB-core/issues/322
                 )
    return ts_one


def DataArray_to_TimeSeries(datablock_xr):
    """
    convert an xarray.DataArray with time dimension to a hydrolib TimeSeries object
    """
    #TODO: clean up these first lines of code and add description to docstring?
    locationname = datablock_xr.attrs['locationname']
    bcvarname = datablock_xr.name
    refdate_str = datablock_xr.time.encoding['units']
    
    if datablock_xr.dims != ('time',):
        raise Exception(f"ERROR: datablock_xr provided to DataArray_to_TimeSeries has dimensions {datablock_xr.dims} while ('time') is expected")
    
    #get datablock and concatenate with relative time data
    datablock_np = datablock_xr.to_numpy()[:,np.newaxis]
    timevar_sel_rel = date2num(pd.DatetimeIndex(datablock_xr.time.to_numpy()).to_pydatetime(),units=refdate_str,calendar='standard')
    datablock_incltime = np.concatenate([timevar_sel_rel[:,np.newaxis],datablock_np],axis=1)
    
    # Each .bc file can contain 1 or more timeseries, in this case one for each support point
    ts_one = TimeSeries(name=locationname,
                        quantityunitpair=[QuantityUnitPair(quantity="time", unit=refdate_str),
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
    poly_pd_xy = poly_pd[nondata_cols]
    poly_pd_data = pd.DataFrame({'data':poly_pd.drop(nondata_cols,axis=1).values.tolist()})
    poly_pd_polyobj = pd.concat([poly_pd_xy,poly_pd_data],axis=1)
    pointsobj_list = poly_pd_polyobj.T.apply(dict).tolist() #TODO: maybe faster with list iteration
    polyobject = PolyObject(metadata={'name':name,'n_rows':poly_pd.shape[0],'n_columns':poly_pd.shape[1]}, points=pointsobj_list)
    if content is not None:
        polyobject.description = {'content':content}
    return polyobject


def forcinglike_to_DataArray(forcingobj): #TODO: would be convenient to have this as a method of ForcingModel objects (Timeseries/T3D/etc): https://github.com/Deltares/HYDROLIB-core/issues/307
    """
    convert a hydrolib forcing like object (like Timeseries, T3D, Harmonic, etc) to a xarray DataArray.
    #TODO: clean up code (maybe split for T3D/Timeseries/Astronomic/etc objects separately) and add doc
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
    data_xr_var = xr.DataArray(datablock_data, name=var_quantity, dims=dims)
    data_xr_var.attrs['locationname'] = forcingobj.name
    data_xr_var.attrs['units'] = var_unit
    if 'depth' in dims:
        data_xr_var['depth'] = forcingobj.vertpositions
        #data_xr_var['depth'].attrs['positive'] == 'up' #TODO: maybe add this attribute
    if 'time' in dims:
        time_unit = forcingobj.quantityunitpair[0].unit.lower()
        data_xr_var['time'] = cftime.num2pydate(datablock_all[:,0], units=time_unit)
        data_xr_var['time'].encoding['units'] = time_unit #check tz conversion if eg '+01:00' is present in time_unit
        data_xr_var['time'].encoding['calendar'] = 'standard'
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
        data_xr_var.attrs[key] = str(forcingobj.__dict__[key])
    
    return data_xr_var


def forcinglike_to_DataFrame(forcingobj):
    """
    convert a hydrolib forcing like object (like Timeseries, T3D, Astronomic, etc) to a pandas DataFrame. Mostly via the forcinglike_to_DataArray() method. #TODO: Astronomic is also supported, what more?
    
    Parameters
    ----------
    forcingobj : hydrolib ForcingModel.forcing object (Timeseries/T3D etc)
        DESCRIPTION.

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
    
    #check if forcingmodel instead of T3D/TimeSeries is provided
    if isinstance(forcingobj, ForcingModel):
        raise Exception('ERROR: instead of supplying a ForcingModel, provide a ForcingObject (Timeseries/T3D etc), by doing something like ForcingModel.forcing[0]')
    
    allowed_instances = (T3D, TimeSeries, Astronomic)
    if not isinstance(forcingobj, allowed_instances):
        raise Exception(f'ERROR: supplied input is not one of: {allowed_instances}')
    
    """ #TODO: old complex code, remove this once timezone is properly supported in forcinglike_to_DataArray()
    #convert_time : boolean, optional
    #    Convert time column from unit (e.g. minutes since date) to datetime index and drop the time column. Has no effect if there is no time column in the forcingobject. The default is True.
    QUP_list = [(QUP.quantity,QUP.unit,QUP.vertpositionindex) for QUP in forcingobj.quantityunitpair]
    columns_MI = pd.MultiIndex.from_tuples(QUP_list,names=dict(forcingobj.quantityunitpair[0]).keys())
    df_data = pd.DataFrame(forcingobj.datablock,columns=columns_MI)
    df_data.index.name = forcingobj.name
    colnames_quantity = df_data.columns.get_level_values(level=0)
    if convert_time and ('time' in colnames_quantity): #this converts time to a datetime index
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
    """
    if isinstance(forcingobj, (T3D, TimeSeries)):
        data_xr_var = forcinglike_to_DataArray(forcingobj)
        df_data = data_xr_var.to_pandas()
    else: #for Astronomic, but might also work well for other objects without time dimension
        columns = [f'{QUP.quantity} [{QUP.unit}]' for QUP in forcingobj.quantityunitpair]
        df_data = pd.DataFrame(forcingobj.datablock,columns=columns)
        df_data = df_data.set_index(columns[0])
       
    return df_data


def pointlike_to_DataFrame(pointlike,drop_emptycols=True):
    """
    convert a hydrolib object with points (like PolyObject, XYZModel and possibly others) to a pandas DataFrame.

    Parameters
    ----------
    pointlike : TYPE
        Hydrolib-core object with point objects.

    Returns
    -------
    pointlike_pd : TYPE
        DESCRIPTION.
        
    Example:
        polyfile_object = PolyFile(file_pli)
        data_pol_pd_list = [pointlike_to_DataFrame(polyobj) for polyobj in polyfile_object.objects]

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



