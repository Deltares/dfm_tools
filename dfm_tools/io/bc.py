# -*- coding: utf-8 -*-
"""
Created on Mon Apr 20 21:01:08 2020

@author: veenstra
"""






def write_bcfile(filepath_bc, data_pd, pliname, refdate, tzone=0, mode='w', data_format='%6.2f'):
    #import numpy as np
    
    data_pd_out = data_pd.copy()
    pd_columns = data_pd_out.columns.tolist()
    if 'time' in pd_columns:
        times_wrtref_min = (data_pd_out['time']-refdate).dt.total_seconds()/60
        if tzone >= 0:
            tzone_sign = '+'
        else:
            tzone_sign = '-'
        times_unit = 'minutes since %s %s%02d:00'%(refdate.strftime('%Y-%m-%d %H:%M:%S'), tzone_sign, tzone)
        data_pd_out['time'] = times_wrtref_min
    
    with open(filepath_bc,mode) as file_bc:
        file_bc.write('[forcing]\n')
        file_bc.write('Name                            = %s\n'%(pliname))
        file_bc.write('Function                        = timeseries\n')
        file_bc.write('Time-interpolation              = linear\n')
        for iC, pd_colname in enumerate(pd_columns):
            if pd_colname in ['time']:
                unit = times_unit
            elif pd_colname in ['dischargebnd','discharge','lateral_discharge']:
                unit = 'm3/s'
            else:
                unit = '-'
            
            file_bc.write('Quantity                        = %s\n'%(pd_colname))
            file_bc.write('Unit                            = %s\n'%(unit))
        
    #    for timestamp, line_time, line_value in zip(times_pd, times_wrtref_min, values):
    #        file_bc.write('%15d %15E\n'%(line_time, line_value))
    data_pd_out.to_csv(filepath_bc, header=None, index=None, sep='\t', mode='a', float_format=data_format)
    with open(filepath_bc,'a') as file_bc:
        file_bc.write('\n')



