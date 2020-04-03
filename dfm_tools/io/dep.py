"""
GNU GENERAL PUBLIC LICENSE
	      Version 3, 29 June 2007

dfm_tools are post-processing tools for Delft3D FM
Copyright (C) 2020 Deltares

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.


INFORMATION
This script is part of dfm_tools: https://github.com/openearth/dfm_tools
Check the README.rst on github for other available functions
Check the tests folder on github for example scripts (this is the dfm_tools pytest testbank)
Check the pptx and example figures in (created by the testbank): N:/Deltabox/Bulletin/veenstra/info dfm_tools

Author and source
Willem Ottevanger
willem.ottevanger@deltares.nl
https://svn.oss.deltares.nl/repos/openearthtools/trunk/python/OpenEarthTools/openearthtools/io/delft3d/dep.py

read/write Delft3D-FLOW *.dep files
"""


import numpy as np
class Dep(object):
    """Create a Delft3D dep file
    # Create a dep grid
    dep = dep()
    # Load a dep from file
    grid = Grid.fromfile('filename.grd')
    dep  = Dep.read('filename.dep',grid.shape)
    # Write dep to file
    Dep.write(dep,'filename.dep')
    """
    def __init__(self, *args, **kwargs):
        self.properties = {}
        self.shape = None
        self.dep   = None
        
    def copy(self):
        copy = Dep()
        copy.shape = self.shape
        copy.val = self.val.copy()
        
        return copy
        
    @staticmethod
    def read(filename, gridshape, **kwargs):
        dep = Dep()
        with open(filename, 'r') as f:
           strings = f.read()
           f.close
        dep.val = np.array([float(s) for s in strings.split()])
        dep.val[dep.val == -999.0] = np.nan
        dep.shape = (gridshape[0]+1, gridshape[1]+1)
        dep.val = np.reshape(dep.val, dep.shape)
        return dep

    def write(dep, filename, **kwargs):
        dep.val[np.isnan(dep.val)] = -999.0
        with open(filename, 'w') as f:
           for n in range(0,dep.shape[0]):
              for m in range(0,dep.shape[1]/12):
                 f.write(''.join(("   % .7E" % x) for x in dep.val[n,12*m:12*m+12]))
                 f.write('\n')
              if dep.shape[1]>dep.shape[1]/12*12:
                 f.write(''.join(("   % .7E" % x) for x in dep.val[n,dep.shape[1]/12*12:dep.shape[1]]))
                 f.write('\n')
           f.close   
        dep.val[dep.val == -999.0] = np.nan
