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

import pandas as pd
import numpy as np
import datetime as dt
from netCDF4 import num2date
import warnings

def write_timfile(filename, datablock, header, converttime=False, refdate=None, tzone=0, float_format='%6.2f'):
    """
    

    Parameters
    ----------
    filename : TYPE
        DESCRIPTION.
    datablockss : Pandas dataframe
        The DataFrame should contain time variable in the first column. Column names will not be parsed.
    header : list of file header lines
        Header lines will be printed to file
    converttime : boolean True/False, optional
        DESCRIPTION, Convert first data column from datetime to minutes since refdate. The default in False
    refdate : datetime object, optional
        DESCRIPTION. The default is None.
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
    
    warnings.warn(DeprecationWarning('the function dfm_tools.io.tim.write_timfile() is deprecated, please use the new hydrolib alternative (in development).')) #TODO: remove code (there is no hydrolib-core alternative available yet)
    
    if type(datablock) is not list:
        datablock = [datablock]
        #metadatas = [metadatas]
    
    #clear file
    with open(filename, 'w') as file_tim:
        pass
    
    data_pd_out = datablock.copy()
    if converttime:
        if refdate is None:
            raise Exception('datetime as units in input data, but no refdate provided')
        print('datetime values are replaced by minutes since')
        try:
            times_wrtref_min = (data_pd_out.iloc[:,0]-refdate).dt.total_seconds()/60
        except TypeError:
            raise Exception('Failure to convert time units. Please check that refdate is a valid datetime object and first column of dataset contains valid datetime objects.')
        
        #replace datetime values by minutes since
        data_pd_out.iloc[:,0] = times_wrtref_min
        
    with open(filename, 'a') as file_tim:
        for line in header:
            file_tim.write('%s\n'%(line))
    data_pd_out.to_csv(filename, header=None, index=None, sep='\t', mode='a', float_format=float_format)
    with open(filename,'a') as file_tim:
        file_tim.write('\n')



def read_timfile(filename, converttime=False, refdate=None):
    """
    

    Parameters
    ----------
    filename : *.tim file name.
    converttime : boolean True/False, optional
        DESCRIPTION. The default is False. First column must be time in minues since refdate.
    refdate : datetime object, optional
        DESCRIPTION. The default is None
    
    Raises
    ------
    Exception
        DESCRIPTION.

    Returns
    -------
    datablock : Pandas dataframe
        The DataFrame contains time variable in the first column. Column names are not parsed.
    header : list
        list of header lines from file

    """
    
    warnings.warn(DeprecationWarning('the function dfm_tools.io.tim.read_timfile() is deprecated, please use the new hydrolib alternative (in development).')) #TODO: remove code (there is no hydrolib-core alternative available yet)
    
    with open(filename, 'r') as bc:
        page = bc.readlines()

    header = np.array([],dtype=str)
    for iL, line in enumerate(page):    
        if '*' not in line:#line.lower()[0].isnumeric():
            line_datastart = iL
            break
        header = np.append(header,line.strip('\n'))
        
            
    #read data block
    data_block_pd = pd.read_csv(filepath_or_buffer=filename, header=None, skiprows=line_datastart, delim_whitespace=True, engine='python')
    if converttime:
        time_minutes = data_block_pd.iloc[:,0]
        data_nc_times_pdtd = pd.to_timedelta(time_minutes, unit='m')
        try:
            data_nc_datetimes = (refdate + data_nc_times_pdtd)#.to_pydatetime()
        except TypeError:
            raise Exception('Failure to convert time units. Please check that refdate is a valid datetime object and first column of dataset contains minutes since refdate.')
        print('Converting times to datetime format...')
        data_block_pd.iloc[:,0] = data_nc_datetimes
                    

    return data_block_pd, header
