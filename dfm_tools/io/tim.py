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
    raise DeprecationWarning('the function dfm_tools.write_timfile() is deprecated, please use the new hydrolib alternative: https://github.com/Deltares/dfm_tools/blob/301-convert-timmodel-to-pandasdataframe/tests/examples/preprocess_hydrolib_readtim.py.') #TODO: remove code


def read_timfile(filename, converttime=False, refdate=None):
    raise DeprecationWarning('the function dfm_tools.read_timfile() is deprecated, please use the new hydrolib alternative: https://github.com/Deltares/dfm_tools/blob/301-convert-timmodel-to-pandasdataframe/tests/examples/preprocess_hydrolib_readtim.py.') #TODO: remove code
    
