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


Created on Mon Apr 20 20:34:30 2020

@author: veenstra
"""



def read_noosfile(file_noos, datetime_format='%Y%m%d%H%M', get_header=False, na_values=None):
    #import datetime as dt
    #import numpy as np
    import pandas as pd
    
    noosheader = []
    noosheader_dict = {}
    with open(file_noos) as f:
        for linenum, line in enumerate(f, 0):
            if '#' in line:
                noosheader.append(line)
                comment_stripped = line.strip('#').strip().split(': ')
                if len(comment_stripped) == 1:
                    if comment_stripped[0] != '':
                        noosheader_dict[comment_stripped[0]] = ''
                else:
                    noosheader_dict[comment_stripped[0].strip()] = comment_stripped[1].strip()
            else:
                startdata = linenum
                break

    
    content_pd = pd.read_csv(file_noos, header=startdata-1, delim_whitespace=True, names=['times_str','values'], na_values=na_values)
    noos_datetime = pd.to_datetime(content_pd['times_str'],format=datetime_format)
    noosdata_pd = pd.DataFrame({'datetime':noos_datetime, 'values':content_pd['values']})
    
    if get_header:
        return noosdata_pd, noosheader_dict
    else:
        return noosdata_pd





def write_noosfile(filename, pd_data, metadata=None, na_values=None, float_format='%8.3f'):
    import pandas as pd
    #import numpy as np
    
    with open('%s'%filename,'w') as file_noos:
        pass
    
    with open('%s'%filename,'a') as file_noos:
        file_noos.write('# ----------------------------------------------------------------\n')
        if metadata is not None:
            if type(metadata) == pd.DataFrame:
                col_headers = metadata.columns.tolist()
            elif type(metadata) == dict:
                col_headers = list(metadata.keys())
            else:
                raise Exception('metadata should be of type dict or pandas.DataFrame')
            #col_headerswidth = np.max([len(x) for x in col_headers])
            for iC, pdcol in enumerate(col_headers):
                if metadata[pdcol] == '':
                    file_noos.write('# {}\n'.format(pdcol))
                else:
                    file_noos.write('# {:30} : {}\n'.format(pdcol,metadata[pdcol]))
        if na_values is None:
            file_noos.write('# {:30} : {}\n'.format('missing values',''))
        else:
            file_noos.write('# {:30} : {}\n'.format('missing values',na_values))
        file_noos.write('# ----------------------------------------------------------------\n')
        
    pd_data.to_csv('%s'%filename, header=None, index=None, sep='\t', mode='a', date_format='%Y%m%d%H%M', float_format=float_format, na_rep=na_values)

