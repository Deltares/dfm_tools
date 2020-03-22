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

def corner2center(cor):
   """cen = corner2center(cor) calculates the value of the center
   of the pixels by avareging the surrounding 4 corner values
   for arrays, or the surrounding 2 corner for vectors or [1,n]
   or [n,1 2D arrays. corner2center works both for coordinate meshes
   as well as data values defined on those meshes.
   
      >>> corner2center([1,3,5])
      array([ 2.,  4.])
      
      >>> corner2center([[1,3,5]])
      array([[ 2.,  4.]])
      
      >>> corner2center([[1],[3],[5]])
      array([[ 2.],
             [ 4.]])
      
      >>> corner2center([[1,3,5],[2,6,10]])
      array([[ 3.,  6.]])
      """
   cor = numpy.asarray(cor)
   shp = numpy.shape(cor)
   if   len(shp)==1 and len(cor)<2:
       raise ValueError('at least 2 elements required')
   elif len(shp)==1:
       cen = numpy.zeros(tuple([shp[0]-1,1]))
       cen = cor[ :-1] + cor[1:  ] 
       return cen/2.
   elif len(shp)==2 and shp[0]==1:
       cen = numpy.zeros(tuple([1,shp[1]-1]))
       cen[:,:] = cor[:, :-1] + cor[:,1:  ] 
       return cen/2.       
   elif len(shp)==2 and shp[1]==1:
       cen = numpy.zeros(tuple([shp[0]-1,1]))
       cen[:,:] = cor[ :-1,:] + cor[1:  ,:] 
       return cen/2.           
   elif len(shp)==2:       
       cen = numpy.zeros(tuple(numpy.array(shp)-1))
       cen[:,:] = cor[ :-1, :-1] + \
                  cor[ :-1,1:  ] + \
                  cor[1:  ,1:  ] + \
                  cor[1:  , :-1] 
   elif len(shp)>3:
       raise NotImplementedError('only 1D and 2D arrays implemented, only intervals and pixels, no voxels')
                  
   return cen/4.

if __name__ == '__main__':
    import doctest
    doctest.testmod()