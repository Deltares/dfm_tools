# -*- coding: utf-8 -*-
"""
dfm_tools are post-processing tools for Delft3D FM
Copyright (C) 2020 Deltares. All rights reserved.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

All names, logos, and references to "Deltares" are registered trademarks of
Stichting Deltares and remain full property of Stichting Deltares at all times.
All rights reserved.


INFORMATION
This script is part of dfm_tools: https://github.com/openearth/dfm_tools
Check the README.rst on github for other available functions
Check the tests folder on github for example scripts (this is the dfm_tools pytest testbank)
Check the pptx and example figures in (created by the testbank): N:/Deltabox/Bulletin/veenstra/info dfm_tools

Created on Thu Apr  1 12:03:47 2021

@author: buckman
"""






def write_timfile(filename, datablocks, metadatas, converttime=False, refdate=None, tzone=0, float_format='%6.2f'):
    """
    

    Parameters
    ----------
    filename : TYPE
        DESCRIPTION.
    datablocks : (list of) MultiIndex pandas dataframe
        The DataFrame should have MultiIndex column names with quantity names on level 0 and unit names on level 1.
    metadatas : (list of) dict
        the metadata dictionary contains all metadata (file header)
    converttime : covert column 'datetime' to minutes since refdate, optional
        DESCRIPTION, The default in False
    refdate : datetime object, optional
        DESCRIPTION. The default is None.
    tzone : TYPE, optional
        DESCRIPTION. The default is 0.
    data_format : TYPE, optional
        DESCRIPTION. The default is '%6.2f'.

    Raises
    ------
    Exception
        DESCRIPTION.

    Returns
    -------
    None.

    """
    
    if type(datablocks) is not list:
        datablocks = [datablocks]
        metadatas = [metadatas]
    
    #clear file
    with open(filename, 'w') as file_tim:
        pass
    
    for metadata_dict, data_pd in zip(metadatas, datablocks):
        data_pd_out = data_pd.copy()
        if converttime:
            if refdate is None:
                raise Exception('datetime as units in input data, but no refdate provided')
            print('datetime values are replaced by minutes since')
            times_wrtref_min = (data_pd[('datetime')]-refdate).dt.total_seconds()/60
            if tzone >= 0:
                tzone_sign = '+'
            else:
                tzone_sign = '-'
            times_unit = 'minutes since %s %s%02d:00'%(refdate.strftime('%Y-%m-%d %H:%M:%S'), tzone_sign, tzone)
            
            #replace datetime values by minutes since
            data_pd_out[(times_unit)] = times_wrtref_min
            data_pd_out.pop('datetime')
            
            #move coverted datetime column to first column and ammend header
            col = data_pd_out.pop(times_unit)
            data_pd_out.insert(0,times_unit,col)
            metadata_dict['* Column 1: '] = ('Time (minutes) w.r.t. refdate='+str(refdate))
            
        #check for number of columns
        if len(metadata_dict) != len(data_pd_out.columns):
            print("WARNING: The number of data columns does not match the number of header columns for file "+filename)

        with open(filename, 'a') as file_tim:
            for key in metadata_dict.keys():
                file_tim.write('%s%s\n'%(key,metadata_dict[key]))
        data_pd_out.to_csv(filename, header=None, index=None, sep='\t', mode='a', float_format=float_format)
        with open(filename,'a') as file_tim:
            file_tim.write('\n')



def read_timfile(filename, converttime=False):
    """
    

    Parameters
    ----------
    filename : *.tim file name.
    converttime : TYPE, optional
        DESCRIPTION. The default is False. 
        Time column header must be specified in following format: 'Time (minutes) w.r.t. refdate=yyyy-MM-dd hh:mm:ss'

    Raises
    ------
    Exception
        DESCRIPTION.

    Returns
    -------
    datablocks : (list of) MultiIndex pandas dataframe
        The DataFrame should have MultiIndex column names with quantity names on level 0 and unit names on level 1.
    metadatas : (list of) dict
        the metadata dictionary contains all metadata (file header)

    """
    import pandas as pd
    import numpy as np
    import datetime as dt
    from netCDF4 import num2date
    
    with open(filename, 'r') as bc:
        page = bc.readlines()
    lines_header = np.array([],dtype=int)
    lines_data = np.array([],dtype=int)
    metadata = {}
    quantity = np.array([],dtype=str)
    for iL, line in enumerate(page):
        if '*' in line.lower():
            lines_header = np.append(lines_header,iL)
            line_part = line.partition(': ')
            quantity = np.append(quantity,line_part[-1].strip('\n'))
            meta_name = ''.join(line_part[0:-1])
            metadata[meta_name] = quantity[iL]
        if line.lower()[0].isnumeric():
            lines_data= np.append(lines_data,iL)
            
    line_datastart = len(lines_header)
    #read data block
    data_block_pd = pd.read_csv(filepath_or_buffer=filename, names=quantity, skiprows=line_datastart, delim_whitespace=True, engine='python')
    if converttime:
        time_id = np.where(data_block_pd.columns.get_level_values(0).str.contains('Time'))[0]
        if len(time_id)==1:
            time_minutes = data_block_pd.iloc[:,time_id].values[:,0]
            units_str = units_str = data_block_pd.columns.get_level_values(0)[time_id][0]
            
            time_units_list = units_str.split(' ')             
            try:
                refdate_str = '%s %s'%(time_units_list[3].replace('refdate=',''), time_units_list[4])
                refdate = dt.datetime.strptime(refdate_str,'%Y-%m-%d %H:%M:%S')
            except:
                raise Exception('Failure reading reference time. Check that time column header is in format: "Time (minutes) w.r.t. refdate=yyyy-MM-dd hh:mm:ss"')
            try:
                data_nc_times_pdtd = pd.to_timedelta(time_minutes, unit='m')
                data_nc_datetimes = (refdate + data_nc_times_pdtd)#.to_pydatetime()
                print('retrieving original timezone succeeded, no conversion to UTC/GMT applied')
            except:
                print('retrieving original timezone failed, using num2date output instead')
                data_nc_datetimes = num2date(times=time_minutes, units='m', only_use_cftime_datetimes=False, only_use_python_datetimes=True)
        data_block_pd['datetime'] = data_nc_datetimes
            
    metadata_list = []
    datablock_list = []
    
    metadata_list.append(metadata)
    datablock_list.append(data_block_pd)
                    

    return datablock_list, metadata_list