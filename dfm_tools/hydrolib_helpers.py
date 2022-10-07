# -*- coding: utf-8 -*-
"""
Created on Tue Aug 23 13:36:44 2022

@author: veenstra
"""

import datetime as dt
import pandas as pd
import cftime
import numpy as np
import hydrolib

def forcingobject_to_dataframe(forcingobj, convert_time=True):
    """
    #TODO: would be convenient to have this as a method of ForcingModel objects (Timeseries/T3D/etc), or maybe as method of the ForcingModel (returning a list of DataFrames): https://github.com/Deltares/HYDROLIB-core/issues/307

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
         file_bc = Path(r'p:\11208053-004-kpp2022-rmm1d2d\C_Work\09_Validatie2018_2020\dflowfm2d-rmm_vzm-j19_6-v2d\boundary_conditions\2018\flow\rmm_discharge_laterals_20171220_20190101_MET.bc')
         m = ForcingModel(file_bc)
         df_data_list = [forcingobject_to_dataframe(forcingobj, convert_time=True) for forcingobj in m.forcing]

    """
    if isinstance(forcingobj, hydrolib.core.io.bc.models.ForcingModel):
        raise Exception('ERROR: instead of supplying a ForcingModel, provide a ForcingObject (Timeseries/T3D etc), by doing something like ForcingModel.forcing[0]')
    QUP_list = [(QUP.quantity,QUP.unit) for QUP in forcingobj.__dict__['quantityunitpair']] #TODO: generating MultiIndex can probably be more elegant (e.g. getting names from QUP list), but I do not know how
    columns_MI = pd.MultiIndex.from_tuples(QUP_list,names=['quantity','unit'])
    df_data = pd.DataFrame(forcingobj.__dict__['datablock'],columns=columns_MI)
    df_data.index.name = forcingobj.__dict__['name']
    colnames_quantity = df_data.columns.get_level_values(level=0)
    if convert_time and ('time' in colnames_quantity): #this converts time to a datetime index
        time_colid = colnames_quantity.get_loc('time')
        time_unit = df_data.columns.get_level_values(level=1)[time_colid]
        df_data.index = cftime.num2pydate(df_data.iloc[:,time_colid],units=time_unit)
        df_data.index.name = forcingobj.__dict__['name'] #again with new index
        #timezone was converted to GMT, re-adjust timezone if needed
        timeunit_sincedatetimetz = time_unit.split('since ')[1]
        tzone_minutes = cftime._parse_date(timeunit_sincedatetimetz)[-1]
        df_data.index = df_data.index.tz_localize('GMT')
        df_data.index = df_data.index.tz_convert(dt.timezone(dt.timedelta(minutes=tzone_minutes)))
        #drop original time column
        df_data = df_data.drop(labels='time',level=0,axis=1)
    return df_data


def polyobject_to_dataframe(PolyObject, dummy=None):
    """
    

    Parameters
    ----------
    PolyObject : TYPE
        DESCRIPTION.

    Returns
    -------
    poly_pd : TYPE
        DESCRIPTION.
        
    Example:
        polyfile_object = read_polyfile(file_ldb,has_z_values=False)
        poly_pd_list = [polyobject_to_dataframe(PO) for PO in polyfile_object['objects']]

    """
    xvals_pd = pd.DataFrame({'x':[p.x for p in PolyObject.points]})
    yvals_pd = pd.DataFrame({'y':[p.y for p in PolyObject.points]})
    datavals_pd = pd.DataFrame([p.data for p in PolyObject.points])
    poly_pd = pd.concat([xvals_pd,yvals_pd,datavals_pd],axis=1)
    poly_pd[poly_pd==dummy] = np.nan

    return poly_pd


def xyzmodel_to_dataframe(XYZModel):
    """
    
    """
    #TODO: more generic would be:
    #for key in data_xyz.points[0].dict().keys():
    #    data_pd_list.append(pd.DataFrame({key:[p[key] for p in data_xyz.points]}))
    # but p[key] is not possible, only p.x etc
    
    xvals_pd = pd.DataFrame({'x':[p.x for p in XYZModel.points]})
    yvals_pd = pd.DataFrame({'y':[p.y for p in XYZModel.points]})
    zvals_pd = pd.DataFrame({'z':[p.z for p in XYZModel.points]})
    comments_pd = pd.DataFrame({'comment':[p.comment for p in XYZModel.points]})
    xyz_pd = pd.concat([xvals_pd,yvals_pd,zvals_pd,comments_pd],axis=1)

    return xyz_pd


