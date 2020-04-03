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

Source
https://svn.oss.deltares.nl/repos/openearthtools/trunk/python/OpenEarthTools/openearthtools/io/delft3d/grid.py
"""


import numpy as np
import datetime

class Grid(object):
    """Create a Delft3D grid file
    # Create an empty grid
    grid = Grid()
    # Load a grid from file
    grid = Grid.fromfile('filename.grd')
    # Write grid to file
    Grid.write(grid,'filename.grd')
    """
    def __init__(self, **kwargs):
        self.properties = kwargs.get('properties', {})
        self.shape = kwargs.get('shape', None)
        self.x     = kwargs.get('x', None)
        self.y     = kwargs.get('y', None)

    @staticmethod
    def fromfile(filename, **kwargs):
        grid = Grid()
        grid.filename = filename
        grid.shape    = None
        rows = []
        with open(filename) as f:
            for line in f:
                line = line.strip()
                # skip comments and empty lines
                if line.startswith('*') or not line:
                    continue
                elif '=' in line:
                    key, value = line.split('=')
                    if 'Coordinate' in key:
                        grid.properties[key.strip()] = value.strip()
                    if 'ETA' in key:
                        row = value.split()
                        n, row = row[0], row[1:]
                        while len(row) < grid.shape[1]:
                            line = f.readline()
                            row.extend(line.split())
                        rows.append(row)
                # Read grid size
                elif grid.shape == None:
                    # line should contain size
                    # convert to nrow x ncolumns
                    grid.shape = tuple(np.array(line.split()[::-1], dtype='int'))
                    assert len(grid.shape) == 2, "Expected shape (2,), got {}, (subgrids not supported)".format(grid.shape)
                    # also read next line
                    line = f.readline()
                    grid.properties['xori'], grid.properties['yori'], grid.properties['alfori'] = np.array(line.split(), dtype='float')
        # rows now contains
        #[X
        # Y]
        data = np.array(rows, dtype="double")
        assert (data.shape[0], data.shape[1]) == (grid.shape[0]*2, grid.shape[1]), \
        "Expected shape of data:{} , got {}".format((grid.shape[0]*2, grid.shape[1]), (data.shape[0], data.shape[1]))
        X, Y = data[:grid.shape[0],:], data[grid.shape[0]:,:]
        grid.x = np.ma.MaskedArray(X,mask=X==0.) # apply standard nodatavalue of 0
        grid.y = np.ma.MaskedArray(Y,mask=Y==0.) # apply standard nodatavalue of 0
        return grid
                    
    def write(self, filename, **kwargs):

        f = open(filename,'w')
        
        mn = np.shape(self.x)

        f.write('* Created at ' + str(datetime.datetime.now()) + '\n')
        f.write('* by OpenEarth Tools \n')
        f.write('* $HeadURL: https://svn.oss.deltares.nl/repos/openearthtools/trunk/python/OpenEarthTools/openearthtools/io/delft3d/grid.py $\n')
        f.write('* $Revision: 16202 $\n')
        f.write('Coordinate System=' + self.properties['Coordinate System'] + '\n')
        f.write(str(mn[1]) + ' ' + str(mn[0]) + '\n')
        f.write(str(self.properties['xori']) + ' ' + str(self.properties['yori']) + ' ' + str(self.properties['alfori']) + '\n')
        
        nrow = int(np.ceil(mn[1]/5)) # max 5 m-values per line
        
        for n in range(mn[0]):
           ib   = 0
           ie   = min(5,mn[1])
           nstr = "%5g" % (n + 1)

           f.write(' ETA=' + nstr+' '+ ' '.join(("%f" % x) for x in self.x[n][ib:ie]) + '\n')
           for i in range(1,nrow):
              ib = i*5
              ie = min((i+1)*5,mn[1])
              f.write('            ' + ' '.join(("%f" % x) for x in self.x[n][ib:ie]) + '\n')

        for n in range(mn[0]):
           ib   = 0
           ie   = min(5,mn[1])
           nstr = "%5g" % (n + 1)

           f.write(' ETA=' + nstr+' '+ ' '.join(("%f" % x) for x in self.y[n][ib:ie]) + '\n')
           for i in range(1,nrow):
              ib = i*5
              ie = min((i+1)*5,mn[1])
              f.write('            ' + ' '.join(("%f" % x) for x in self.y[n][ib:ie]) + '\n')

        f.close()
