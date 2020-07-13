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

Created on Mon Apr 20 21:01:08 2020

@author: veenstra
"""






def write_bcfile(filename, datablocks, metadatas, refdate=None, tzone=0, float_format='%6.2f'):
    """
    

    Parameters
    ----------
    filename : TYPE
        DESCRIPTION.
    datablocks : (list of) MultiIndex pandas dataframe
        The DataFrame should have MultiIndex column names with quantity names on level 0 and unit names on level 1.
    metadatas : (list of) dict
        the metadata dictionary contains all metadata, except for quantity and unit (they are in the dataframe columns.
    refdate : TYPE, optional
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
    #import numpy as np
    
    if type(datablocks) is not list:
        datablocks = [datablocks]
        metadatas = [metadatas]
    
    #clear file
    with open(filename, 'w') as file_bc:
        pass
    
    for metadata_dict, data_pd in zip(metadatas, datablocks):
        data_pd_out = data_pd.copy()
        pd_columns = data_pd.columns.tolist()
        if ('time','datetime') in pd_columns:
            if refdate is None:
                raise Exception('datetime as units in input data, but no refdate provided')
            time_colid = pd_columns.index(('time','datetime'))
            print('datetime values are replaced by minutes since')
            times_wrtref_min = (data_pd[('time','datetime')]-refdate).dt.total_seconds()/60
            if tzone >= 0:
                tzone_sign = '+'
            else:
                tzone_sign = '-'
            times_unit = 'minutes since %s %s%02d:00'%(refdate.strftime('%Y-%m-%d %H:%M:%S'), tzone_sign, tzone)
            
            #rename datetime unit
            pd_columns_lev1 = [y for x,y in pd_columns]
            pd_columns_lev1[time_colid] = times_unit
            data_pd_out.columns.set_levels(pd_columns_lev1,level=1,inplace=True)
            
            #replace datetime values by minutes since
            data_pd_out[('time',times_unit)] = times_wrtref_min
        pd_out_columns = data_pd_out.columns.tolist()
       
        with open(filename, 'a') as file_bc:
            file_bc.write('[forcing]\n')
            for key in metadata_dict.keys():
                file_bc.write('%-32s= %s\n'%(key,metadata_dict[key]))
            for iC, pd_colname in enumerate(pd_out_columns):
                quan = pd_colname[0]
                unit = pd_colname[1]
                """
                if pd_colname in ['time']:
                    unit = times_unit
                elif pd_colname in ['dischargebnd','discharge','lateral_discharge']:
                    unit = 'm3/s'
                else:
                    unit = '-'
                """
                file_bc.write('Quantity                        = %s\n'%(quan))
                file_bc.write('Unit                            = %s\n'%(unit))
            
        #    for timestamp, line_time, line_value in zip(times_pd, times_wrtref_min, values):
        #        file_bc.write('%15d %15E\n'%(line_time, line_value))
        data_pd_out.to_csv(filename, header=None, index=None, sep='\t', mode='a', float_format=float_format)
        with open(filename,'a') as file_bc:
            file_bc.write('\n')



def read_bcfile(filename):

    import pandas as pd
    import numpy as np
    
    with open(filename, 'r') as bc:
        page = bc.readlines()
    lines_blockstarts = np.array([],dtype=int)
    lines_quantity = np.array([],dtype=int)
    lines_unit = np.array([],dtype=int)
    lines_n = len(page)
    for iL, line in enumerate(page):
        if '[forcing]' in line.lower():
            lines_blockstarts = np.append(lines_blockstarts,iL)
        if 'quantity' in line.lower():
            lines_quantity = np.append(lines_quantity,iL)
        if 'unit' in line.lower():
            lines_unit = np.append(lines_unit,iL)
    if len(lines_quantity) != len(lines_unit):
        raise Exception('something wrong with bc-file, amount of quantity lines is not the same as amount of unit lines')
    
    metadata_list = []
    datablock_list = []
    for iB in range(len(lines_blockstarts)):
        line_blockstart = lines_blockstarts[iB]
        if iB == len(lines_blockstarts)-1:
            line_blockend = lines_n
        else:
            line_blockend = lines_blockstarts[iB+1]
        
        line_blockquantities = lines_quantity[(lines_quantity>line_blockstart) & (lines_quantity<line_blockend)]
        line_blockunits = lines_unit[(lines_unit>line_blockstart) & (lines_unit<line_blockend)]
        line_datastart = line_blockunits[-1]+1
        
        #get metadata (without quantities and units)
        metadata_str = page[line_blockstart:line_blockquantities[0]]
        metadata = {}
        for metadata_line in metadata_str:
            if metadata_line.startswith('['):
                pass#metadata[metadata_line.strip()] = ''
            else:
                meta_name = metadata_line.split('=')[0].strip()
                meta_val = metadata_line.split('=')[1].strip()
                metadata[meta_name] = meta_val
        metadata_list.append(metadata)
        
        #define column names of quantities and units
        colnames = []
        for line_quan,line_un in zip(line_blockquantities,line_blockunits):
            quan_val = page[line_quan].split('=')[1].strip()
            un_val = page[line_un].split('=')[1].strip()
            colnames.append((quan_val,un_val))
        
        #read data block
        data_block_pd = pd.read_csv(filepath_or_buffer=filename, names=colnames, skiprows=line_datastart, delim_whitespace=True, skipfooter=lines_n-line_blockend, engine='python')
        datablock_list.append(data_block_pd)

                    

    return datablock_list, metadata_list