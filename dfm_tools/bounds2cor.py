"""CF conventions subsidiaries: convert between corners, centers and bounds"""

__version__ = "$Revision: 9965 $"

#  Copyright notice
#   --------------------------------------------------------------------
#   Copyright (C) 2014 Deltares
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

def bounds2cor(bnds):
   """cor = cor2bounds(bounds) calculates the corners array
   from bounds array where bounds has one extra trailing 
   dimension which is 4 for 2D arrays and 2 for vectors 
   or [1,n] or [n,1] 2D arrays. Inverse of cor2bounds().
   
   >>> bounds2cor([[ 1.,  2.],[ 2.,  3.]])
   array([ 1.,  2.,  3.])
   
   >>> bounds2cor([[[ 1.,  2.],[ 2.,  3.]]])
   array([[ 1.,  2.,  3.]])
   
   >>> bounds2cor([[[ 1.,  2.]],[[ 2.,  3.]]])
   array([[ 1.],
          [ 2.],
          [ 3.]])
   
   >>> bounds2cor([[[ 1.,  2.,  5.,  4.],[ 2.,  3.,  6.,  5.]]])
   array([[ 1.,  2.,  3.],
          [ 4.,  5.,  6.]])
   """
   bnds = numpy.asarray(bnds)
   shp = numpy.shape(bnds)

   if   len(shp)==2 and shp[1]==1:
       raise ValueError('at least 2 elements required')
   elif len(shp)==2:              
       cor = numpy.zeros((shp[0]+1)) # 1D vector
       cor[ :-1] = bnds[:,0]
       cor[1:  ] = bnds[:,1]
   elif len(shp)==3 and shp[0]==1 and shp[2]==2: # 2D vector
       cor = numpy.zeros((shp[0]+0,shp[1]+1))
       cor[:, :-1] = bnds[:,:,0]
       cor[:,1:  ] = bnds[:,:,1]
   elif len(shp)==3 and shp[1]==1 and shp[2]==2: # 2D vector transposed        
       cor = numpy.zeros((shp[0]+1,shp[1]+0))
       cor[ :-1,:] = bnds[:,:,0]
       cor[1:  ,:] = bnds[:,:,1]            
   elif len(shp)==3:              
       cor = numpy.zeros(tuple(numpy.array(shp)[0:2]+1)) # 2D matrix
       cor[ :-1, :-1] = bnds[:,:,0]
       cor[ :-1,1:  ] = bnds[:,:,1] # redundant
       cor[1:  ,1:  ] = bnds[:,:,2] # redundant
       cor[1:  , :-1] = bnds[:,:,3] # redundant
   elif len(shp)>3:       
       raise NotImplementedError('only 1D and 2D arrays implemented')
   return cor
   
if __name__ == '__main__':
    import doctest
    doctest.testmod()