# -*- coding: utf-8 -*-
"""
Created on Tue Aug 23 13:36:44 2022

@author: veenstra
"""

import datetime as dt
import pandas as pd
import cftime

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
