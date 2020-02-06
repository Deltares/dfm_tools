# WIP python version of Matlab counterpart patch2tri.m

__version__ = "$HeadURL: https://svn.oss.deltares.nl/repos/openearthtools/trunk/python/OpenEarthTools/openearthtools/io/netcdf/netCDF4_tutorial_grid_lat_lon_curvilinear.py $" + "$Revision: 8907 $"

import numpy
import matplotlib.mlab
find = matplotlib.mlab.find
sum  = numpy.sum

#for i in [1]:
def patch2tri(corx,cory,mapn):
   """[triangles,map3,ntyp]=patch2tri(x,y,mapn) 
   Map grid with triangles, quadrilaterals, pentagons and hexagons to a mesh with 
   triangles only. triangles is a (?,3) array with 0-based integer pointers into 
   vectors (x,y), ready for use with matplotlib.tri.Triangulation() and 
   mayavi.mlab.triangular_mesh(). mapn is a (?,6) array analogous to triangles, 
   padded to width 6 to accomodate hexagons. For triangles, quadralatals and 
   pentagons respectively the last 3,2,1 elements are -1. map3 is  an array that 
   duplicates center values into triangle center values.
   
   Example:  6--7--8
   	     | \|/ |
   	     3--4--5
   	     |  |  |
   	     0--1--2
   
      # ipython --pylab
      import numpy as np
      import patch2tri
      from mayavi import mlab
      
      corx =  [0,1,2,0,1,2,0,1,2] 
      cory =  [0,0,0,1,1,1,2,2,2]
      corz =  [0,1,2,5,4,3,6,7,8]
      mapn = [[0,1,4,3,-1,-1],[1,2,5,4,-1,-1],[3,4,6,-1,-1,-1],[4,7,6,-1,-1,-1],[4,5,8,-1,-1,-1],[4,8,7,-1,-1,-1]]
      tri  = patch2tri.patch2tri(corx,cory,mapn)[0] # tri
      mlab.triangular_mesh(corx,cory,corz,tri,scalars=corz)
   
   NOTE: mapn, triangles and map3 are 0-based (!) index arrays into x and y (!)
   NOTE: Currently implemented only for mapn containing triangles or quadrilaterals
   
   >>> patch2tri([0,1,2,0,1,2,0,1,2],[0,0,0,1,1,1,2,2,2],[[0,1,4,3,-1,-1],[1,2,5,4,-1,-1],[3,4,6,-1,-1,-1],[4,7,6,-1,-1,-1],[4,5,8,-1,-1,-1],[4,8,7,-1,-1,-1]])[0] # tri
   array([[ 3.,  4.,  6.],
          [ 4.,  7.,  6.],
          [ 4.,  5.,  8.],
          [ 4.,  8.,  7.],
          [ 0.,  1.,  4.],
          [ 1.,  2.,  5.],
          [ 0.,  4.,  3.],
          [ 1.,  5.,  4.]])
   >>> patch2tri([0,1,2,0,1,2,0,1,2],[0,0,0,1,1,1,2,2,2],[[0,1,4,3,-1,-1],[1,2,5,4,-1,-1],[3,4,6,-1,-1,-1],[4,7,6,-1,-1,-1],[4,5,8,-1,-1,-1],[4,8,7,-1,-1,-1]])[1] # map3
   array([[ 2.],
          [ 3.],
          [ 4.],
          [ 5.],
          [ 0.],
          [ 1.],
          [ 0.],
          [ 1.]])
   >>> patch2tri([0,1,2,0,1,2,0,1,2],[0,0,0,1,1,1,2,2,2],[[0,1,4,3,-1,-1],[1,2,5,4,-1,-1],[3,4,6,-1,-1,-1],[4,7,6,-1,-1,-1],[4,5,8,-1,-1,-1],[4,8,7,-1,-1,-1]])[2] # ntyp
   array([0, 0, 4, 2, 0, 0])
   >>> patch2tri([0,2,3.5,4.5,5.5,6.5,5.5,4.5,3.5,2,0,0],[1,1,1,2,2,3,4,4,3,3,5,3],[[9,10,11,-1,-1,-1],[0,1,9,11,-1,-1],[1,2,3,8,9,-1],[3,4,5,6,7,8]])[0]
   array([[  9.,  10.,  11.],
          [  0.,   1.,   9.],
          [  0.,   9.,  11.],
          [  9.,   1.,   2.],
          [  9.,   2.,   8.],
          [  2.,   3.,   8.],
          [  8.,   3.,   7.],
          [  3.,   4.,   7.],
          [  7.,   4.,   6.],
          [  4.,   5.,   6.]])
   >>> patch2tri([0,2,3.5,4.5,5.5,6.5,5.5,4.5,3.5,2,0,0],[1,1,1,2,2,3,4,4,3,3,5,3],[[9,10,11,-1,-1,-1],[0,1,9,11,-1,-1],[1,2,3,8,9,-1],[3,4,5,6,7,8]])[1]
   array([[ 0.],
          [ 1.],
          [ 1.],
          [ 2.],
          [ 2.],
          [ 2.],
          [ 3.],
          [ 3.],
          [ 3.],
          [ 3.]])
   >>> patch2tri([0,2,3.5,4.5,5.5,6.5,5.5,4.5,3.5,2,0,0],[1,1,1,2,2,3,4,4,3,3,5,3],[[9,10,11,-1,-1,-1],[0,1,9,11,-1,-1],[1,2,3,8,9,-1],[3,4,5,6,7,8]])[2]   
   array([0, 0, 1, 1, 1, 1])
           
   See also: matplotlib.tri.Triangulation
             mayavi.mlab.triangular_mesh    
             http://oss.deltares.nl/web/delft3d/d-flow-flexible-mesh/-/message_boards/view_message/220151
             https://groups.google.com/forum/#!topic/ugrid-interoperability/IIFRDFsKgJM
   
      """
   # 1st doctest
   #corx =  [0,1,2,0,1,2,0,1,2] 
   #cory =  [0,0,0,1,1,1,2,2,2]
   #z    =  [0,1,2,5,4,3,6,7,8]
   #mapn = [[0,1,4,3,-1,-1],[1,2,5,4,-1,-1],[3,4,6,-1,-1,-1],[4,7,6,-1,-1,-1],[4,5,8,-1,-1,-1],[4,8,7,-1,-1,-1]]
   
   corx = numpy.asarray(corx)
   cory = numpy.asarray(cory)
   mapn = numpy.asarray(mapn) # -1 # go from 1-based to 0-based
   
   # determine topology of patch objects
   nface      = numpy.asarray(sum(mapn > -1,1)) # map is 0-based now
   
   ntyp = numpy.asarray([0,0,sum(nface==3),sum(nface==4),sum(nface==5),sum(nface==6)])
   # quadrangle become 2 triangles each (4-2)
   # pentagon   become 3 triangles each (5-2)
   # hexahon    become 4 triangles each (6-2)
   
# WIP                  
# ---------------
###ntri = sum(ntyp*[0,0,1,2,3,4])
   ntri = sum(ntyp*[0,0,1,2,0,0]) # WIP
   tri  = numpy.zeros([ntri,3])-1 # 0-based internally !
   map3 = numpy.zeros([ntri,1])-1 # 0-based internally !
      
   # 3: re-use existing triangles
   ind = find([nface==3])
   tri[0:len(ind),:]=mapn[ind,0:3]
   map3[0:len(ind)] = ind.reshape(len(ind),1)
    
   # 4: quadrilaterals: 2 triangles each
   #    delaunay used to fails for perfect square (R2009a and lower), so we gotta do it ourselves
   msk = ([nface==4])
   n   = ntyp[2]
   
   # Split quad into two triangles directly: cut 'along' shortest diagonal.
   mapnset = mapn[msk]
   ismaindiag = numpy.array(((corx[mapnset[:,0]] - corx[mapnset[:,2]])**2 \
                           + (cory[mapnset[:,0]] - cory[mapnset[:,2]])**2 \
                         ) > (corx[mapnset[:,1]] - corx[mapnset[:,3]])**2 \
                           + (cory[mapnset[:,1]] - cory[mapnset[:,3]])**2)
   ind = find(msk)

# Quads where 2--4 is shortest diagonal: tri's are (1-based) 1-2-4 and 2-3-4
   ntri1 =  sum(ismaindiag.nonzero())
   if ntri1 > 0:
      i1 = ind[find(ismaindiag)]
      mapnset = mapn[i1,:]
      tri [n        :n+  ntri1,    0:3] = mapnset[:,[0,1,3]]
      tri [n  +ntri1:n+2*ntri1,    0:3] = mapnset[:,[1,2,3]]
      map3[n        :n+  ntri1]         = i1.reshape(len(i1),1)
      map3[n  +ntri1:n+2*ntri1]         = i1.reshape(len(i1),1)
      n =            n+2*ntri1

# Quads where 1--3 is shortest diagonal: tri's are 1-2-3 and 1-3-4
   ntri2 = len(ismaindiag) - ntri1
   if ntri2 > 0:
      i2 = ind[find(ismaindiag==False)]
      mapnset = mapn[i2,:]
      tri [n        :n+  ntri2,    0:3] = mapnset[:,[0,1,2]]
      tri [n  +ntri2:n+2*ntri2,    0:3] = mapnset[:,[0,2,3]]
      map3[n        :n+  ntri2]         = i2.reshape(len(i2),1)
      map3[n  +ntri2:n+2*ntri2]         = i2.reshape(len(i2),1)
      n =            n+2*ntri2
   
# WIP                  
# ---------------
# 
# # 5: pentagons: 3 triangles each
# 
# ind = find(nface==5);
# if ~(n== sum(ntyp.*[0 0 1 2 0 0]));error('error after tri + quad');end
# [tri,n,map3] = nface2tri(map,corx,cory,tri,5,n,ind,OPT.debug,'pentagon',map3,OPT.quiet);
# 
# #plot(tri(:,1));hold on;plot(map3,'r');xlim([0 ntri]);pausedisp;clf
# 
# # 6
# 
# ind = find(nface==6);
# if ~(n== sum(ntyp.*[0 0 1 2 3 0]));error('error after tri + quad + pent');end
# [tri,n,map3] = nface2tri(map,corx,cory,tri,6,n,ind,OPT.debug,'hexagon',map3,OPT.quiet);
# 
# #plot(tri(:,1));hold on;plot(map3,'r');xlim([0 ntri]);pausedisp;clf
# 
# # out

   tri  = tri  #+1 # go back from 0-based to 1-based
   map3 = map3 #+1 # go back from 0-based to 1-based

   return tri,map3,ntyp
   
# WIP
# ---------------
# # genericish subsidiary for quad-, pent- and hexagons
# 
# def function = nface2tri(Map,X,Y,tri,type,n,ind,debug,txt,map3,quiet)
# 
#     order = type-2;
# 
#     for i=1:length(ind)
#     
#         pointers = Map(ind(i),1:type);
#         x        = X(pointers);
#         y        = Y(pointers);
#         trilocal = delaunay(x,y); % sometimes fails, and does not always yield 3 triangles
#         
#         if size(trilocal,1) > order
#            if ~quiet
#                warning([txt,' is not divided onto ',num2str(order),' triangles but ',num2str(size(trilocal,1)),': triangle(s) ingnored']);
#            end
#            trilocal = trilocal(1:order,:);
#            if debug
#               plot(x,y,'c-o','linewidth',10)
#               pausedisp
#            end
#         elseif size(trilocal,1) < order
#            if ~quiet
#                warning([txt,' is not divided onto ',num2str(order),' triangles but ',num2str(size(trilocal,1)),': triangle(s) duplicated']);
#            end
#            
#            trilocal = repmat(trilocal,[order 1]);
#            
#            if debug
#               plot(x,y,'c-o','linewidth',10);
#               pausedisp
#            end
#         end
#         
#         tri (n+1:n+order,:) = pointers(trilocal);
#         map3(n+1:n+order,:) = ind(i);
#         n                   = n + order;
#         
#     end       
#     
#     yield tri,map3,n 
    

if __name__ == '__main__':
    import doctest
    doctest.testmod()    