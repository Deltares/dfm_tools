# -*- coding: utf-8 -*-
"""
Created on Tue Aug 23 13:36:44 2022

@author: veenstra
"""


def forcingobject_to_dataframe(forcingobj, convert_time=True):
    import datetime as dt
    import pandas as pd
    import cftime
    
    #TODO: would be convenient to have this as a method of ForcingModel objects (Timeseries/T3D/etc), or maybe as method of the ForcingModel (returning a list of DataFrames): https://github.com/Deltares/HYDROLIB-core/issues/307
    QUP_list = [(QUP.quantity,QUP.unit) for QUP in forcingobj.__dict__['quantityunitpair']] #TODO: generating MultiIndex can probably be more elegant (e.g. getting names from QUP list), but I do not know how
    columns_MI = pd.MultiIndex.from_tuples(QUP_list,names=['quantity','unit'])
    df_data = pd.DataFrame(forcingobj.__dict__['datablock'],columns=columns_MI)
    if convert_time: #this converts time to a datetime index
        time_colid = df_data.columns.get_level_values(level=0).get_loc('time')
        time_unit = df_data.columns.get_level_values(level=1)[time_colid]
        df_data.index = cftime.num2pydate(df_data.iloc[:,time_colid],units=time_unit)
        df_data.index.name = forcingobj.__dict__['name']
        #timezone was converted to GMT, re-adjust timezone if needed
        timeunit_sincedatetimetz = time_unit.split('since ')[1]
        tzone_minutes = cftime._parse_date(timeunit_sincedatetimetz)[-1]
        df_data.index = df_data.index.tz_localize('GMT')
        df_data.index = df_data.index.tz_convert(dt.timezone(dt.timedelta(minutes=tzone_minutes)))
        #drop original time column
        df_data = df_data.drop(labels='time',level=0,axis=1)
    return df_data
