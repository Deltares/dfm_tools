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
