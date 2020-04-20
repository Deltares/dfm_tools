# -*- coding: utf-8 -*-
"""
Created on Fri Apr 17 13:02:48 2020

@author: veenstra
"""
#WARNING, THIS IS WORK IN PROGRESS, NOT TESTED YET

def convert_xycoords(loc_in_xy,epsg_in,epsg_out): #default epsg_out value is RD
    from pyproj import Proj, transform
    import pandas as pd
    import numpy as np
    
    if epsg_in == epsg_out:
        loc_out_xy = loc_in_xy
    else:
        inProj = Proj(init='epsg:%i'%(epsg_in))
        outProj = Proj(init='epsg:%i'%(epsg_out))
        
        if isinstance(loc_in_xy, pd.DataFrame):
            x_in = loc_in_xy['RDx'].values
            y_in = loc_in_xy['RDy'].values
            x_out,y_out = transform(inProj,outProj,x_in,y_in)
            data_forpd = np.empty([len(x_out),2])
            data_forpd[:,0] = x_out
            data_forpd[:,1] = y_out
            loc_out_xy = pd.DataFrame(data=data_forpd, columns=['RDx','RDy'])
        else:
            x_in = loc_in_xy[0]
            y_in = loc_in_xy[1]
            x,y = transform(inProj,outProj,x_in,y_in)
            loc_out_xy = [x,y]
        
    return loc_out_xy









def write_bcfile(data_pd, dir_bcfile, pliname, refdate, tzone=0, mode='w', data_format='%6.2f'):
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
    
    with open(dir_bcfile,mode) as file_bc:
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
    data_pd_out.to_csv(dir_bcfile, header=None, index=None, sep='\t', mode='a', float_format=data_format)
    with open(dir_bcfile,'a') as file_bc:
        file_bc.write('\n')





