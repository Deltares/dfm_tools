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

Created on Wed Feb 12 01:55:57 2020

@author: veenstra
"""
import warnings
import numpy as np
import pandas as pd



class Polygon:
    def __init__(self, data, name, comments):
        #import numpy as np
        self.data = data
        self.name = name
        self.comments = comments
        raise DeprecationWarning('the function dfm_tools.io.polygon.Polygon() is deprecated, please use the new hydrolib alternative.') #TODO: remove this code
        #self.line_array = np.c_[self.x, self.y]
        
    def fromfile(file_pol, pd_output=False, tekmap_output=False):#, obj_output=False):
        """
        reads a pli boundary into an array
        

        Parameters
        ----------
        file_pol : TYPE
            DESCRIPTION.
        pd_output : TYPE, optional
            pd_output==True: create pandas array of tekal file, including comments
            pd_output==False (default): create list of numpy arrays and list of polygon names
            DESCRIPTION. The default is False.
        

        Raises
        ------
        Exception
            DESCRIPTION.

        Returns
        -------
        TYPE
            DESCRIPTION.

        """
        raise DeprecationWarning('the function dfm_tools.polygon.Polygon.fromfile() is deprecated, please use the new hydrolib alternative. Example script: https://github.com/Deltares/dfm_tools/blob/main/tests/examples/preprocess_hydrolib_readwritepol.py') #TODO: remove this definition

        with open(file_pol) as plifile:
            lines = plifile.readlines()
            pol_data_list = []
            pol_name_list = []
            pol_comment_list = []
            pol_comment_temp = []
            pol_data_pd_list = []
            #pol_object_list = []
            iLine=0
            while iLine<=len(lines)-1:
                line_str = lines[iLine].split()
                try:
                    float(line_str[0]) #this line should start with a string, so this should fail
                    raise Exception('ERROR: line_str[0] is not a string: (%s, on line %d), polygon file format might be wrong'%(line_str, iLine))
                except:
                    pass
                iLine=iLine+1
                if line_str[0].startswith('*'): #comments
                    if ('column' in lines[iLine-1].lower()) or ('kolom' in lines[iLine-1].lower()):
                        line_str_joined = ' '.join(line_str)
                        pol_comment_temp.append(line_str_joined)
                else:
                    #the start of a new polygon
                    pol_name_list.append(line_str[0]) #get name from name line_str
                    pol_comment_list.append(pol_comment_temp)
                    
                    #convert shape numbers to int
                    line_shp = [int(x) for x in lines[iLine].split()] #eg [1150, 10, 10] if it is a tekal mapfile
                    iLine=iLine+1
                    len_pol=int(line_shp[0])
                    
                    #retrieve data numbers as floats
                    data_pol_str = lines[iLine:iLine+len_pol]
                    iLine = iLine+len_pol
                    data_pol = np.array([x.split() for x in data_pol_str],dtype=float)
                    data_pol_nanbool = data_pol==999.999
                    data_pol[data_pol_nanbool] = np.nan
                    pol_data_list.append(data_pol)
                    if len(pol_comment_temp)==data_pol.shape[1]: #expected format
                        pol_data_pd = pd.DataFrame(data_pol,columns=pol_comment_temp)
                        if len(line_shp) == 2:
                            if 'date' in pol_comment_temp[0].lower() and 'time' in pol_comment_temp[1].lower():
                                #first convert columns to integers, since they were before
                                #pol_data_pd.iloc[:,0] = pol_data_pd.iloc[:,0].astype(int)
                                #pol_data_pd.iloc[:,1] = pol_data_pd.iloc[:,1].astype(int)
                                try:
                                    #add datetime column with parsed dates
                                    pol_data_datetime = pd.to_datetime(pol_data_pd.iloc[:,0]*1000000+pol_data_pd.iloc[:,1],format='%Y%m%d%H%M%S') #this does not work anymore
                                    #pol_data_datetime = (pol_data_pd.iloc[:,0]*1000000+pol_data_pd.iloc[:,1]).astype('datetime64[ns]').dt.round(freq='S')
                                    pol_data_pd.insert(0, 'datetime', pol_data_datetime)
                                except:
                                    warnings.warn('conversion from date/time column to datetime failed, incorrect format of first two columns?')
                                    print('conversion from date/time column to datetime failed, incorrect format of first two columns?')
                        elif len(line_shp) == 3: #tekal 2D map files
                            data_pol_map = data_pol.reshape([line_shp[0]//line_shp[2],line_shp[2],line_shp[1]])
                        elif len(line_shp) == 4: #tekal 2D map files, with fourth dummy in shape array
                            if line_shp[0]//line_shp[2]-line_shp[3]>=1e-9:
                                raise Exception('tekal map file shape header is not consistent with itself')
                            data_pol_map = data_pol.reshape([line_shp[0]//line_shp[2],line_shp[2],line_shp[1]])
                        else:
                            raise Exception('tekal file should have 2 or 3 dimensions, but it has this shape: %s'%(line_shp))
                    else:
                        pol_data_pd = pd.DataFrame(data_pol)
                    pol_data_pd_list.append(pol_data_pd)
                    
                    #pol_object = Polygon(data=data_pol, name=line_str[0], comments=pol_comment_temp)
                    #pol_object_list.append(pol_object)
                    
                    # reset comment list for next block
                    pol_comment_temp = []
        if pd_output is True:
            return pol_data_pd_list
        elif tekmap_output:
            if len(pol_comment_list) > 1:
                raise Exception('multiple tekal maps not anticipated yet')
            return data_pol_map, pol_comment_list[0]
        #elif obj_output is True:
        #    return pol_object_list
        else:
            return pol_data_list, pol_name_list, pol_comment_list
            
