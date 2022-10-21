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
    VectorQuantityUnitPairs,
    T3D,
    TimeSeries,
    Astronomic,
)


def DataArray_to_T3D(datablock_xr):   
    """
    convert an xarray.DataArray with time and depth dimension to a hydrolib T3D object
    """
    #if isinstance(datablock_xr,(tuple,list)): #if multiple arguments, call itself multiple time and return ruple of results
    #    print('WARNING: part of the code is not valid yet') #TODO: needed: a T3Dvector to (T3D,T3D) converter and a (T3D,T3D) to T3Dvector converter
    #    ts_one_tuple = tuple([DataArray_to_T3D(x) for x in datablock_xr])
    #    return ts_one_tuple
    
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
                 vertinterpolation='linear', #TODO: make these parameters user defined (via attrs)
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


def T3Dvector_to_T3Dtuple(forcingobj):
    """
    converts a T3D with vector element quantities (like a ux/uy combination) to a tuple of T3D objects (one for each vector element)
    """
    allowed_instances = (T3D)
    if not isinstance(forcingobj, allowed_instances):
        raise Exception(f'ERROR: supplied input is not one of: {allowed_instances}')

    #if not isinstance(forcingobj.quantityunitpair[1],VectorQuantityUnitPairs):
    if not hasattr(forcingobj.quantityunitpair[1],'elementname'): #vector with elementname
        raise Exception('non-vector quantity supplied to T3Dvector_to_T3Dtuple()')
    
    datablock_all = np.array(forcingobj.datablock)
    datablock_time = datablock_all[:,:1] #select first column
    datablock_data = datablock_all[:,1:] #select all columns except first one. 
    
    QUP_time = forcingobj.quantityunitpair[0]
    var_quantity_list = forcingobj.quantityunitpair[1].elementname
    #var_unit = forcingobj.quantityunitpair[1].quantityunitpair[0].unit #assuming unit of all quantities in T3Dvector is the same
    list_T3D = []
    for quan in var_quantity_list:
        QUP_quan = [colQUP for colidx,colQUP in enumerate(forcingobj.quantityunitpair[1].quantityunitpair) if colQUP.quantity==quan]
        colidx_quan = [colidx for colidx,colQUP in enumerate(forcingobj.quantityunitpair[1].quantityunitpair) if colQUP.quantity==quan]
        datablock_data_onequan = datablock_data[:,colidx_quan] #subset every nquan column, starting at iQ (gives all columns in case of nquan=1 or in case of one column)
        
        datablock_incltime = np.concatenate([datablock_time,datablock_data_onequan],axis=1)
        
        # Each .bc file can contain 1 or more timeseries, in this case one for each support point
        ts_one = T3D(name=forcingobj.name,
                     offset=forcingobj.offset,
                     factor=forcingobj.factor,
                     vertpositions=forcingobj.vertpositions,
                     vertinterpolation=forcingobj.vertinterpolation,
                     vertpositiontype=forcingobj.vertpositiontype,
                     quantityunitpair=[QUP_time]+QUP_quan,
                     timeinterpolation=forcingobj.timeinterpolation,
                     datablock=datablock_incltime.tolist(),
                     )
        list_T3D.append(ts_one)
    
    return tuple(list_T3D)


def T3Dtuple_to_T3Dvector(*forcingobjects,vectorname='uxuyadvectionvelocitybnd'):
    """
    converts a tuple of T3D objects (one for each vector element) to a T3D with vector element quantities (like a ux/uy combination)
    """
    
    allowed_instances = (T3D)
    forcingobj0 = forcingobjects[0]
    forcingobj0_datablock_all = np.array(forcingobj0.datablock)
    forcingobj0_time = forcingobj0_datablock_all[:,:1] #select first column
    for iF,forcingobj in enumerate(forcingobjects):
        if not isinstance(forcingobj, allowed_instances):
            raise Exception(f'ERROR: supplied input (argument {iF+1}) is not one of: {allowed_instances}')
        #if isinstance(forcingobj.quantityunitpair[1],VectorQuantityUnitPairs):
        if hasattr(forcingobj.quantityunitpair[1],'elementname'): #vector with elementname
            raise Exception(f'vector quantity supplied to T3Dtuple_to_T3Dvector() (argument {iF+1}), supply tuple of non-vector T3D instead')
    
    nquan = len(forcingobjects)
    T3Dvector_datashape = (forcingobj0_datablock_all.shape[0],nquan*(forcingobj0_datablock_all.shape[1]-1))
    datablock_data_allobjs = np.concatenate([forcingobj0_time,np.zeros(T3Dvector_datashape)],axis=1)
    QUP_quan = np.empty(T3Dvector_datashape[1],dtype=object)
    for iF,forcingobj in enumerate(forcingobjects):
        forcingobj_datablock_all = np.array(forcingobj.datablock)
        
        #check if shape of datablock and contents of time are equal
        if not iF==0:
            forcingobj_time = forcingobj_datablock_all[:,:1] #select first column
            if forcingobj_datablock_all.shape != forcingobj0_datablock_all.shape:
                raise Exception(f'shape of datablock in forcingobj {iF+1} ({forcingobj_datablock_all.shape}) does not correspond to shape in first forcingobj ({forcingobj0_datablock_all.shape})')
            if not (forcingobj_time==forcingobj0_time).all():
                raise Exception(f'time array of datablock in forcingobj {iF+1} does not correspond to time array in first forcingobj')
            if forcingobj.name != forcingobj0.name:
                raise Exception(f'name of forcingobj {iF+1} ({forcingobj.name}) does not correspond to name of first forcingobj ({forcingobj0.name})')
            if forcingobj.offset != forcingobj0.offset:
                raise Exception(f'offset of forcingobj {iF+1} ({forcingobj.offset}) does not correspond to offset of first forcingobj ({forcingobj0.offset})')
            if forcingobj.factor != forcingobj0.factor:
                raise Exception(f'factor of forcingobj {iF+1} ({forcingobj.factor}) does not correspond to factor of first forcingobj ({forcingobj0.factor})')
            if forcingobj.vertpositions != forcingobj0.vertpositions:
                raise Exception(f'vertpositions of forcingobj {iF+1} ({forcingobj.vertpositions}) does not correspond to vertpositions of first forcingobj ({forcingobj0.vertpositions})')
            if forcingobj.vertinterpolation != forcingobj0.vertinterpolation:
               raise Exception(f'vertinterpolation of forcingobj {iF+1} ({forcingobj.vertinterpolation}) does not correspond to vertinterpolation of first forcingobj ({forcingobj0.vertinterpolation})')
            if forcingobj.vertpositiontype != forcingobj0.vertpositiontype:
               raise Exception(f'vertpositiontype of forcingobj {iF+1} ({forcingobj.vertpositiontype}) does not correspond to vertpositiontype of first forcingobj ({forcingobj0.vertpositiontype})')
            if forcingobj.timeinterpolation != forcingobj0.timeinterpolation:
               raise Exception(f'timeinterpolation of forcingobj {iF+1} ({forcingobj.timeinterpolation}) does not correspond to timeinterpolation of first forcingobj ({forcingobj0.timeinterpolation})')
        
        #fill datablock_data_allobjs
        forcingobj_datablock = forcingobj_datablock_all[:,1:]
        datablock_data_allobjs[:,iF+1::nquan] = forcingobj_datablock
        QUP_quan[iF::nquan] = forcingobj.quantityunitpair[1:]
    QUP_time = forcingobj0.quantityunitpair[0]
    
    quan_list = [forcingobj.quantityunitpair[1].quantity for forcingobj in forcingobjects]
    QUP_quan_vector = VectorQuantityUnitPairs(vectorname=vectorname,
                                              elementname=quan_list,
                                              quantityunitpair=QUP_quan.tolist())
    
    # Each .bc file can contain 1 or more timeseries, in this case one for each support point
    T3Dvector = T3D(name=forcingobj0.name,
                    offset=forcingobj0.offset,
                    factor=forcingobj0.factor,
                    vertpositions=forcingobj0.vertpositions,
                    vertinterpolation=forcingobj0.vertinterpolation,
                    vertpositiontype=forcingobj0.vertpositiontype,
                    quantityunitpair=[QUP_time,QUP_quan_vector],
                    timeinterpolation=forcingobj0.timeinterpolation,
                    datablock=datablock_data_allobjs.tolist(),
                    )
    
    return T3Dvector


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
    
    if hasattr(forcingobj.quantityunitpair[1],'elementname'): #vector with elementname
        raise Exception('vector quantity supplied to forcinglike_to_DataArray(), first split in separate T3D objects with T3Dvector_to_T3Dtuple()')
    
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
    datablock_data = datablock_all[:,1:] #select all columns except first one. TODO: convert repeating values to nan? (reverse of ffill/bfill)
    if forcingobj.quantityunitpair[0].quantity == 'astronomic component': #first column in bcfile is astronomic component
        datablock_data = datablock_data.astype(float) #convert str to float in case of "astronomic component"
    datablock_data = datablock_data.squeeze() #drop dimensions of len 1 in case of 1 dimension, eg "waterlevelbnd" (first subsetting over depth dimension)
    
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



