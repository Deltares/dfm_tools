# -*- coding: utf-8 -*-
"""
Created on Tue Aug 23 13:36:44 2022

@author: veenstra
"""

import pandas as pd
import cftime
import numpy as np
import xarray as xr
from cftime import date2num

try: #0.3.1 release
    from hydrolib.core.io.polyfile.models import PolyObject
    from hydrolib.core.io.bc.models import ForcingModel, QuantityUnitPair, VectorQuantityUnitPairs, T3D, TimeSeries, Astronomic
except: #main branch and next release
    from hydrolib.core.io.dflowfm.polyfile.models import PolyObject
    from hydrolib.core.io.dflowfm.bc.models import ForcingModel, QuantityUnitPair, VectorQuantityUnitPairs, T3D, TimeSeries, Astronomic


def Dataset_to_T3D(datablock_xr):  
    """
    convert an xarray.DataArray (is one data_var) or an xarray.Dataset (with one or two data_vars) with time and depth dimension to a hydrolib T3D object
    """
    
    if isinstance(datablock_xr,xr.DataArray):
        vector = False
        data_xr_var0 = datablock_xr
    elif isinstance(datablock_xr,xr.Dataset):
        data_vars = list(datablock_xr.data_vars)
        if len(data_vars)==1: 
            vector = False
            data_xr_var0 = datablock_xr[data_vars[0]]
        elif len(data_vars)==2:
            if not pd.Series(data_vars).isin(['ux','uy']).all():
                raise Exception(f'Dataset with 2 data_vars should contain only ux/uy data_vars, but contains {data_vars}')
            vector = True
            data_xr_var0 = datablock_xr[data_vars[0]]
            data_xr_var1 = datablock_xr[data_vars[1]]
        else:
            raise Exception(f'Dataset should contain 1 or 2 data_vars, but contains {len(data_vars)} variables')
    else:
        raise Exception(f'expected xarray.DataArray or xarray.Dataset, not {type(datablock_xr)}')
    
    #ffill/bfill nan data along over depth dimension (corresponds to vertical extrapolation)
    data_xr_var0 = data_xr_var0.bfill(dim='depth').ffill(dim='depth')
    if vector:
        data_xr_var1 = data_xr_var1.bfill(dim='depth').ffill(dim='depth')
    
    #TODO: clean up these first lines of code and add description to docstring?
    locationname = data_xr_var0.attrs['locationname']
    refdate_str = data_xr_var0.time.encoding['units']
    
    if data_xr_var0.dims != ('time','depth'): #check if both time and depth dimensions are present #TODO: add support for flipped dimensions (datablock_xr.T or something is needed)
        raise Exception(f"ERROR: data_var in provided data_xr has dimensions {data_xr_var0.dims} while ('time','depth') is expected")
    
    #get depth variable and values
    depth_array = data_xr_var0['depth'].to_numpy()
    
    #get datablock and concatenate with relative time data
    if vector:
        data_xr_var0_np = data_xr_var0.to_numpy()
        data_xr_var1_np = data_xr_var1.to_numpy()
        datablock_np = np.stack((data_xr_var0_np,data_xr_var1_np),2).reshape(data_xr_var0_np.shape[0],-1) #merge data with alternating rows
    else:
        datablock_np = data_xr_var0.to_numpy()
    
    timevar_sel_rel = date2num(pd.DatetimeIndex(data_xr_var0.time.to_numpy()).to_pydatetime(),units=refdate_str,calendar='standard')
    datablock_incltime = np.concatenate([timevar_sel_rel[:,np.newaxis],datablock_np],axis=1)
    
    # Each .bc file can contain 1 or more timeseries, in this case one for each support point
    verticalpositions_idx = np.arange(data_xr_var0['depth'].size)+1
    if vector: #vector T3D object
        QUP_quan_list = [QuantityUnitPair(quantity=quan, unit=data_xr_var0.attrs['units'], vertpositionindex=iVP) for iVP in verticalpositions_idx for quan in data_vars]
        QUP_quan_vector = VectorQuantityUnitPairs(vectorname='uxuyadvectionvelocitybnd', #TODO: vectorname from global attr? (then also support other vectors which is not necessary)
                                                  elementname=data_vars,
                                                  quantityunitpair=QUP_quan_list)
        quantityunitpair = [QuantityUnitPair(quantity="time", unit=refdate_str)]+[QUP_quan_vector]
    else: #normal T3D object
        QUP_quan_list = [QuantityUnitPair(quantity=data_xr_var0.name, unit=data_xr_var0.attrs['units'], vertpositionindex=iVP) for iVP in verticalpositions_idx]
        quantityunitpair=[QuantityUnitPair(quantity="time", unit=refdate_str)]+QUP_quan_list
    
    T3D_object = T3D(name=locationname,
                     #offset=0,
                     #factor=1,
                     vertpositions=depth_array.tolist(),
                     vertinterpolation='linear', #TODO: make these parameters user defined (via attrs)
                     vertPositionType='ZDatum',
                     quantityunitpair=quantityunitpair,
                     timeinterpolation='linear',
                     datablock=datablock_incltime.tolist(), #TODO: numpy array is not supported by TimeSeries. https://github.com/Deltares/HYDROLIB-core/issues/322
                     )
    
    return T3D_object


def Dataset_to_TimeSeries(datablock_xr):
    """
    convert an xarray.DataArray or xarray.Dataset with time dimension to a hydrolib TimeSeries object
    """
    if isinstance(datablock_xr,xr.DataArray):
        pass
    if isinstance(datablock_xr,xr.Dataset):
        data_vars = list(datablock_xr.data_vars)
        if len(data_vars)!=1:
            raise Exception('more than one variable supplied in Dataset, not yet possible') #TODO: add support for multiple quantities and for vectors
        datablock_xr = datablock_xr[data_vars[0]]
    
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
    TimeSeries_object = TimeSeries(name=locationname,
                                   quantityunitpair=[QuantityUnitPair(quantity="time", unit=refdate_str),
                                                     QuantityUnitPair(quantity=bcvarname, unit=datablock_xr.attrs['units'])],
                                   timeinterpolation='linear', #TODO: make userdefined via attrs?
                                   datablock=datablock_incltime.tolist(), 
                                   )
    return TimeSeries_object


def DataFrame_to_PolyObject(poly_pd,name,content=None):
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


def forcinglike_to_Dataset(forcingobj, convertnan=False): #TODO: would be convenient to have this as a method of ForcingModel objects (Timeseries/T3D/etc): https://github.com/Deltares/HYDROLIB-core/issues/307
    """
    convert a hydrolib forcing like object (like Timeseries, T3D, Harmonic, etc) to an xarray Dataset with one or more variables.
    
    convertnan: convert depths with the same values over time as the deepest layer to nan (these were created with .bfill() or .ffill()).
    """
    
    #check if forcingmodel instead of T3D/TimeSeries is provided
    if isinstance(forcingobj, ForcingModel):
        raise Exception('ERROR: instead of supplying a ForcingModel, provide a ForcingObject (Timeseries/T3D etc), by doing something like ForcingModel.forcing[0]')
    
    allowed_instances = (T3D, TimeSeries, Astronomic)
    if not isinstance(forcingobj, allowed_instances):
        raise Exception(f'ERROR: supplied input is not one of: {allowed_instances}')
    
    if isinstance(forcingobj, Astronomic):
        var_quantity_list = [x.quantity for x in forcingobj.quantityunitpair[1:]]
        var_unit = [x.unit for x in forcingobj.quantityunitpair[1:]]
    elif hasattr(forcingobj.quantityunitpair[1],'elementname'): #T3D with vector quantity
        var_quantity_list = forcingobj.quantityunitpair[1].elementname
        var_unit_one = forcingobj.quantityunitpair[1].quantityunitpair[0].unit
        var_unit = [var_unit_one,var_unit_one]
    else: #non-vector TimeSeries or T3D
        var_quantity_list = [forcingobj.quantityunitpair[1].quantity]
        var_unit = [forcingobj.quantityunitpair[1].unit]
    nquan = len(var_quantity_list)
    
    if isinstance(forcingobj, T3D):
        dims = ('time','depth')
    elif isinstance(forcingobj, TimeSeries):
        dims = ('time')
    elif isinstance(forcingobj, Astronomic):
        dims = ('astronomic_component')
    
    datablock_all = np.array(forcingobj.datablock)
    datablock_data = datablock_all[:,1:] #select all columns except first one (which is the time column)
    if isinstance(forcingobj, Astronomic):
        datablock_data = datablock_data.astype(float) #convert str to float in case of "astronomic component"
    
    data_xr = xr.Dataset()
    for iQ, var_quantity in enumerate(var_quantity_list):
        datablock_data_onequan = datablock_data[:,iQ::nquan]
        datablock_data_onequan = datablock_data_onequan.squeeze() #drop dimensions of len 1 in case of 1 dimension, eg "waterlevelbnd" (first subsetting over depth dimension)
        
        data_xr_var = xr.DataArray(datablock_data_onequan, name=var_quantity, dims=dims)
        if 'depth' in dims:
            data_xr_var['depth'] = forcingobj.vertpositions
            #data_xr_var['depth'].attrs['positive'] == 'up' #TODO: maybe add this attribute
            if convertnan: #convert ffilled/bfilled values back to nan
                deepestlayeridx = data_xr_var.depth.to_numpy().argmin()
                if deepestlayeridx==0: #sorted from deep to shallow layers
                    bool_nandepths = (data_xr_var==data_xr_var.shift(depth=-1)).all(dim='time')
                else: #sorted from shallow to deep layers
                    bool_nandepths = (data_xr_var==data_xr_var.shift(depth=1)).all(dim='time')
                data_xr_var = data_xr_var.where(~bool_nandepths)
        if 'time' in dims:
            time_unit = forcingobj.quantityunitpair[0].unit.lower()
            data_xr_var['time'] = cftime.num2pydate(datablock_all[:,0], units=time_unit)
            data_xr_var['time'].encoding['units'] = time_unit #check tz conversion if eg '+01:00' is present in time_unit
            data_xr_var['time'].encoding['calendar'] = 'standard'
        if 'astronomic_component' in dims:
            data_xr_var['astronomic_component'] = datablock_all[:,0]
        
        #add attributes
        data_xr_var.attrs['source'] = 'hydrolib-core object converted to xarray.Dataset with dfm_tools.xarray_helpers.forcinglike_to_Dataset()'
        data_xr_var.attrs['locationname'] = forcingobj.name
        data_xr_var.attrs['units'] = var_unit[iQ]
        forcingobj_keys = forcingobj.__dict__.keys()
        for key in forcingobj_keys: #['comments','name','function','offset','factor','vertinterpolation','vertpositiontype','timeinterpolation']: 
            if key in ['datablock','quantityunitpair','vertpositions']: #skipping these since they are in the DataArray already
                continue
            data_xr_var.attrs[key] = str(forcingobj.__dict__[key])
        
        #add DataArray to Dataset
        data_xr[var_quantity] = data_xr_var
    
    return data_xr


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



