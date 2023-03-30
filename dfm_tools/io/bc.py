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

def write_bcfile(filename, datablocks, metadatas, refdate=None, tzone=0, float_format='%6.2f'): #TODO: remove this code, only raises an error
    
    raise DeprecationWarning('the function dfm_tools.io.bc.write_bcfile() is deprecated, please use the new hydrolib alternative. Example script: dfm_tools/tests/examples/CMEMS_interpolate_example.py')


def read_bcfile(filename, converttime=False): #TODO: remove this code, only raises an error

    raise DeprecationWarning('the function dfm_tools.io.bc.read_bcfile() is deprecated, please use the new hydrolib alternative. Example script: dfm_tools/tests/examples/hydrolib_readbc.py')
    