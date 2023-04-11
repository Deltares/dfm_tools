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

"""

def write_timfile(filename, datablock, header, converttime=False, refdate=None, float_format='%6.2f'):
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
    
    raise DeprecationWarning('the function dfm_tools.write_timfile() is deprecated, please use the new hydrolib alternative: data_tim = hcdfm.TimModel(file_tim); tim_pd = dfmt.TimModel_to_DataFrame(data_tim).') #TODO: remove code


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
    
    raise DeprecationWarning('the function dfm_tools.read_timfile() is deprecated, please use the new hydrolib alternative: data_tim = hcdfm.TimModel(file_tim); tim_pd = dfmt.TimModel_to_DataFrame(data_tim).') #TODO: remove code
    
