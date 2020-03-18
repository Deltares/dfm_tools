"""read/write Delft3D-FLOW *.dep files"""

__version__ = "$Revision: 7870 $"

#  Copyright notice
#   --------------------------------------------------------------------
#   Copyright (C) 2013 Deltares
#       Willem Ottevanger
#
#       willem.ottevanger@deltares.nl
#
#   This library is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   This library is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this library.  If not, see <http://www.gnu.org/licenses/>.
#   --------------------------------------------------------------------
#
# This tool is part of <a href="http://www.OpenEarth.eu">OpenEarthTools</a>.
# OpenEarthTools is an online collaboration to share and manage data and
# programming tools in an open source, version controlled environment.
# Sign up to recieve regular updates of this function, and to contribute
# your own tools.

# $Id: dep.py 7870 2012-12-31 13:33:52Z ottevan $
# $Date: 2012-12-31 14:33:52 +0100 (Mon, 31 Dec 2012) $
# $Author: ottevan $
# $Revision: 7870 $
# $HeadURL: https://svn.oss.deltares.nl/repos/openearthtools/trunk/python/OpenEarthTools/openearthtools/io/delft3d/dep.py $
# $Keywords: $

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
