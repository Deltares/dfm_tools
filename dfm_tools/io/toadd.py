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






def read_noosfile(Variables, dir_filename):
    import datetime as dt
    import numpy as np
    
    nan_value = 9999.0 #hard coded in definition file bij *.noos bestanden
    with open(dir_filename) as f:
        for linenum, line in enumerate(f, 0):
            #print linenum, line
            #if linenum >15:
            #    break
            if '#' in line:
                pass
                #print linenum, line
            else:
                startdata = linenum
                #print startdata
                break
    with open(dir_filename) as f:
        for linenum, line in enumerate(f, 0):
            if len(line)>3: #check for last nonempty line
                enddata = linenum
    with open(dir_filename) as f:
        content = f.readlines()[startdata:enddata+1]
    data = np.empty([len(content), 2])
    for iL in range(len(content)):
        if len(content[iL])>2:
            if len(content[iL].split())==2 and str(nan_value) not in content[iL].split()[1]:
                data[iL,0] = content[iL].split()[0]
                data[iL,1] = float(content[iL].split()[1])
            else:
                data[iL,0] = content[iL].split()[0]
                data[iL,1] = np.nan
    location = 'NaN'
    name = 'NaN'
    positionWGS84xy = 'NaN'
    positionRDxy = 'NaN'
    
    times_fromfile = [dt.datetime.strptime(str(int(x)), '%Y%m%d%H%M') for x in data[:,0]]
    timestep_mins_fromfile = (times_fromfile[1]-times_fromfile[0]).total_seconds()/60
    values_fromfile = data[:,1]
    unit = 'm'
    
    return times_fromfile, values_fromfile





def write_noosfile(filename, data_fromhis):
    #copied from arup, not finished yet
    with open('%s.txt'%filename,'w') as file_noos:
        pass
    
    pd_timeval = pd.DataFrame({'times': data_fromhis.var_times})
    
    fig, ax = plt.subplots(figsize=(18,7))
    for iS, stat in enumerate(station):
        ax.plot(data_fromhis.var_times,data_fromhis[:,iS],'-',linewidth=0.8,label='%s (lat %.3f, lon %.3f)'%(stat,stat_coordsy[iS],stat_coordsx[iS]),markersize=2)
        with open('%s.txt'%filename,'a') as file_noos:
            file_noos.write('# ----------------------------------------------------------------\n')
            file_noos.write('# \n')
            file_noos.write('# Timeseries created by Python script getnoos_fromnc.py\n')
            file_noos.write('# Values retrieved from hisdata, locations snapped to cell centers\n')
            file_noos.write('# \n')
            file_noos.write('# Source           : %s\n'%(source_str))
            file_noos.write('# Analysis time    : 202003260000\n')
            file_noos.write('# Time zone        : GMT\n')
            file_noos.write('# Coordinate system: WGS84\n')
            file_noos.write('# Station          : %s\n'%(stat))
            file_noos.write('# x-coordinate     : %.4f\n'%(stat_coordsx[iS]))
            file_noos.write('# y-coordinate     : %.4f\n'%(stat_coordsy[iS]))
            file_noos.write('# \n')
            file_noos.write('# ----------------------------------------------------------------\n')
            file_noos.write('# \n')
            file_noos.write('# Variable     : waterlevel\n')
            file_noos.write('# long_name    : waterlevel\n')
            file_noos.write('# units        : m\n')
            file_noos.write('# missing value: -9999\n')
            file_noos.write('#\n')
        pd_timevalstat = pd.DataFrame({'times': data_fromhis.var_times, 'vals':data_fromhis[:,iS]})
        pd_timeval[stat] = data_fromhis[:,iS]
        pd_timevalstat.to_csv('%s.txt'%filename, header=None, index=None, sep='\t', mode='a', date_format='%Y%m%d%H%M', float_format='%6.3f')
    pd_timeval.to_csv('%s.csv'%filename, index=None, date_format='%Y%m%d%H%M', float_format='%6.3f')




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





