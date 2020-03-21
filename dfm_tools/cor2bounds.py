"""CF conventions subsidiaries: convert between corners, centers and bounds"""

__version__ = "$Revision: 9965 $"

#  Copyright notice
#   --------------------------------------------------------------------
#   Copyright (C) 2013 Deltares
#       Gerben J. de Boer
#
#       gerben.deboer@deltares.nl
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

# $Id: CF.py 9965 2014-01-08 13:59:08Z boer_g $
# $Date: 2014-01-08 14:59:08 +0100 (Wed, 08 Jan 2014) $
# $Author: boer_g $
# $Revision: 9965 $
# $HeadURL: https://svn.oss.deltares.nl/repos/openearthtools/trunk/python/OpenEarthTools/openearthtools/io/netcdf/CF.py $
# $Keywords: $

import numpy

def cor2bounds(cor):
   """bounds = cor2bounds(cor) calculates the bounds array
   where bounds has one extra trailing dimension which is 4 for
   2D arrays and 2 for vectors or [1,n] or [n,1] 2D arrays.
   Inverse of bounds2cor().
   
   >>> cor2bounds([1,2,3])
   array([[ 1.,  2.],
          [ 2.,  3.]])
   
   >>> cor2bounds([[1,2,3]])
   array([[[ 1.,  2.],
           [ 2.,  3.]]])
   
   >>> cor2bounds([[1],[2],[3]])
   array([[[ 1.,  2.]],
   <BLANKLINE>
          [[ 2.,  3.]]])
   
   >>> cor2bounds([[1,2,3],[4,5,6]])
   array([[[ 1.,  2.,  5.,  4.],
           [ 2.,  3.,  6.,  5.]]])
   """
   cor = numpy.asarray(cor)
   shp = numpy.shape(cor)
   if   len(shp)==1 and len(cor)<2:
       raise ValueError('at least 2 elements required')
   elif len(shp)==1:
       bnds = numpy.zeros(tuple([shp[0]-1,2])) # 1D vector
       bnds[:,0] = cor[ :-1]
       bnds[:,1] = cor[1:  ]
   elif len(shp)==2 and shp[0]==1:
       bnds = numpy.zeros(tuple([1,shp[1]-1,2])) # 2D vector
       bnds[:,:,0] = cor[:, :-1]
       bnds[:,:,1] = cor[:,1:  ]
   elif len(shp)==2 and shp[1]==1:
       bnds = numpy.zeros(tuple([shp[0]-1,1,2])) # 2D vector transposed
       bnds[:,:,0] = cor[ :-1,:]
       bnds[:,:,1] = cor[1:  ,:]            
   elif len(shp)==2:              
       bnds = numpy.zeros(tuple(numpy.array(shp)-1)+(4,)) # 2D matrix
       bnds[:,:,0] = cor[ :-1, :-1]
       bnds[:,:,1] = cor[ :-1,1:  ]
       bnds[:,:,2] = cor[1:  ,1:  ]
       bnds[:,:,3] = cor[1:  , :-1]
   elif len(shp)>2:       
       raise NotImplementedError('only 1D and 2D arrays implemented')
   return bnds
   
if __name__ == '__main__':
    import doctest
    doctest.testmod()