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

