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

Created on Fri Feb 14 11:23:12 2020

@author: veenstra
"""

class UGrid:
    """Unstructured grid"""
    def __init__(self, mesh2d_node_x, mesh2d_node_y, mesh2d_face_nodes, verts, mesh2d_node_z=None, edge_verts=None):
        self.mesh2d_node_x = mesh2d_node_x
        self.mesh2d_node_y = mesh2d_node_y
        self.mesh2d_face_nodes = mesh2d_face_nodes
        self.verts = verts
        self.mesh2d_node_z = mesh2d_node_z
        #if mesh2d_node_z is not None:
        #    self.mesh2d_node_z = mesh2d_node_z
        #else:
        #    self.mesh2d_node_z = np.zeros(self.mesh2d_node_x.shape)
        self.edge_verts=edge_verts #can be none?
    @staticmethod
    def fromfile(file_nc):
        import numpy as np
        from netCDF4 import Dataset
        from dfm_tools.get_nc_helpers import get_varname_fromnc, ghostcell_filter
        
        def nodexyfaces2verts(node_x,node_y, faces):
            quatrangles = faces-1 #convert 1-based indexing of cell numbering in ugrid to 0-based indexing
            quatrangles_filled = quatrangles.filled(int(-999)) #necessary since fill value was changed from -999 to -2147483647 and this cannot be handled for some reason
            #https://stackoverflow.com/questions/49640311/matplotlib-unstructered-quadrilaterals-instead-of-triangles
            #https://stackoverflow.com/questions/52202014/how-can-i-plot-2d-fem-results-using-matplotlib
            yz = np.c_[node_x,node_y]
            verts= yz[quatrangles_filled]
            verts[quatrangles.mask==True,:] = np.nan #remove all masked values by making them nan
            return verts
        
        def nodexyfaces2edgeverts(node_x,node_y,edge_nodes, face_x,face_y,edge_faces):
            quatranglese = edge_nodes-1 #convert 1-based indexing of cell numbering in ugrid to 0-based indexing
            quatranglesf = edge_faces-1 #convert 1-based indexing of cell numbering in ugrid to 0-based indexing
            #https://stackoverflow.com/questions/49640311/matplotlib-unstructered-quadrilaterals-instead-of-triangles
            #https://stackoverflow.com/questions/52202014/how-can-i-plot-2d-fem-results-using-matplotlib
            yze = np.c_[node_x,node_y]
            yzf = np.c_[face_x,face_y]
            yzf = np.concatenate([yzf,[[np.nan,np.nan]]]) #add dummy row to index -1, since edge_faces has one 0-value when it is a model border edge. This gets turned into -1 towards quatranglesf (last index), which must be nan eventually
            vertse= yze[quatranglese]
            vertsf= yzf[quatranglesf]
            vertse[quatranglese.mask==True,:] = np.nan #remove all masked values by making them nan
            vertsf[quatranglesf.mask==True,:] = np.nan #remove all masked values by making them nan
            verts_raw = np.concatenate([vertse,vertsf],axis=1)
            edge_verts = verts_raw[:,[0,2,1,3],:] #set ordering in (counter)clockwise direction instead of updown-leftright. second column will sometimes contain nans, then the face is not plotted.
            return edge_verts
        
        data_nc = Dataset(file_nc)
        
        varn_mesh2d_node_x = get_varname_fromnc(data_nc,'mesh2d_node_x',vardim='var')
        varn_mesh2d_node_y = get_varname_fromnc(data_nc,'mesh2d_node_y',vardim='var')
        if varn_mesh2d_node_x is None or varn_mesh2d_node_y is None:
            raise Exception('file does not contain variables "mesh2d_node_x" and "mesh2d_node_y" or similar, are you sure this is an unstructured grid?')
        mesh2d_node_x = data_nc.variables[varn_mesh2d_node_x][:]
        mesh2d_node_y = data_nc.variables[varn_mesh2d_node_y][:]
        
        varn_mesh2d_node_z = get_varname_fromnc(data_nc,'mesh2d_node_z',vardim='var')
        if varn_mesh2d_node_z is not None: # node_z variable is present
            mesh2d_node_z = data_nc.variables[varn_mesh2d_node_z][:]
        else:
            mesh2d_node_z = None

        varn_mesh2d_face_nodes = get_varname_fromnc(data_nc,'mesh2d_face_nodes',vardim='var')
        if varn_mesh2d_face_nodes is not None: # node_z variable is present
            mesh2d_face_nodes = data_nc.variables[varn_mesh2d_face_nodes][:, :]
        else:
            raise Exception('ERROR: provided file does not contain a variable mesh2d_face_nodes or similar:\n%s\nPlease do one of the following:\n- plot grid from *_map.nc file\n- import and export the grid with RGFGRID\n- import and save the gridd "with cellfinfo" from interacter'%(file_nc))
        verts = nodexyfaces2verts(mesh2d_node_x, mesh2d_node_y, mesh2d_face_nodes) #xy coordinates of face nodes
        
        #varn_mesh2d_edge_x = get_varname_fromnc(data_nc,'mesh2d_edge_x',vardim='var')
        varn_mesh2d_edge_faces = get_varname_fromnc(data_nc,'mesh2d_edge_faces',vardim='var')
        if varn_mesh2d_edge_faces is not None: # mesh2d_edge_x (and mesh2d_edge_y) variable is present
            #get coordinates and mapping of edge start/end node
            #mesh2d_edge_x = data_nc.variables[varn_mesh2d_edge_x][:]
            #mesh2d_edge_y = data_nc.variables[get_varname_fromnc(data_nc,'mesh2d_edge_y',vardim='var')][:]
            mesh2d_edge_nodes = data_nc.variables[get_varname_fromnc(data_nc,'mesh2d_edge_nodes',vardim='var')][:]
            #get center coordinates and mapping of two bordering faces
            mesh2d_face_x = data_nc.variables[get_varname_fromnc(data_nc,'mesh2d_face_x',vardim='var')][:]
            mesh2d_face_y = data_nc.variables[get_varname_fromnc(data_nc,'mesh2d_face_y',vardim='var')][:]
            mesh2d_edge_faces = data_nc.variables[get_varname_fromnc(data_nc,'mesh2d_edge_faces',vardim='var')][:]
            #combine to edge_verts
            edge_verts = nodexyfaces2edgeverts(mesh2d_node_x, mesh2d_node_y, mesh2d_edge_nodes, mesh2d_face_x, mesh2d_face_y, mesh2d_edge_faces) #xy coordinates of face nodes
            
        else:
            edge_verts = None
            
        #remove ghost cells from faces and verts
        nonghost_bool = ghostcell_filter(file_nc)
        if nonghost_bool is not None:
            mesh2d_face_nodes = mesh2d_face_nodes[nonghost_bool]
            verts = verts[nonghost_bool]
            if 0: # remove edges from partition boundaries if there are partitions
                part_edges_removebool = (mesh2d_edge_faces==0).any(axis=1) # Array is 1 based indexed, 0 means missing # & (np.in1d(mesh2d_edge_faces[:,0],ghost_removeids-1) | np.in1d(mesh2d_edge_faces[:,1],ghost_removeids-1))
                edge_verts = edge_verts[~part_edges_removebool]

        data_nc.close()
        ugrid = UGrid(mesh2d_node_x, mesh2d_node_y, mesh2d_face_nodes, verts, mesh2d_node_z=mesh2d_node_z, edge_verts=edge_verts)
        return ugrid


    def polygon_intersect(self, line_array, optimize_dist=False, calcdist_fromlatlon=False):
        import numpy as np
        from matplotlib.path import Path
        
        from dfm_tools.testutils import try_importmodule
        try_importmodule(modulename='shapely')
        import shapely
        from shapely.geometry import LineString, Polygon, MultiLineString, Point
        from dfm_tools.get_nc import calc_dist_pythagoras, calc_dist_haversine

        print('defining celinlinebox')
        
        line_section = LineString(line_array)
        
        verts_xmax = np.nanmax(self.verts[:,:,0].data,axis=1)
        verts_xmin = np.nanmin(self.verts[:,:,0].data,axis=1)
        verts_ymax = np.nanmax(self.verts[:,:,1].data,axis=1)
        verts_ymin = np.nanmin(self.verts[:,:,1].data,axis=1)
        
        if not optimize_dist:
            cellinlinebox_all_bool = (((np.min(line_array[:,0]) <= verts_xmax) &
                                       (np.max(line_array[:,0]) >= verts_xmin)) &
                                      ((np.min(line_array[:,1]) <= verts_ymax) & 
                                       (np.max(line_array[:,1]) >= verts_ymin))
                                      )
        elif type(optimize_dist) in [int,float]: #not properly tested and documented
            #calculate angles wrt x axis
            angles_wrtx = []
            nlinecoords = line_array.shape[0]
            for iL in range(nlinecoords-1):
                dx = line_array[iL+1,0] - line_array[iL,0]
                dy = line_array[iL+1,1] - line_array[iL,1]
                angles_wrtx.append(np.rad2deg(np.arctan2(dy,dx)))
            angles_toprev = np.concatenate([[90],np.diff(angles_wrtx),[90]])
            angles_wrtx_ext = np.concatenate([[angles_wrtx[0]-90],np.array(angles_wrtx),[angles_wrtx[-1]+90]])
            angtot_wrtx = angles_wrtx_ext[:-1] + 0.5*(180+angles_toprev)
            #distance over xy-axis from original points
            dxynewpoints = optimize_dist * np.array([np.cos(np.deg2rad(angtot_wrtx)),np.sin(np.deg2rad(angtot_wrtx))]).T
            newpoints1 = line_array+dxynewpoints
            newpoints2 = line_array-dxynewpoints
            pol_inpol = np.concatenate([newpoints1, np.flip(newpoints2,axis=0)])
            pol_inpol_path = Path(pol_inpol)
            bool_all = []
            for iC in range(self.verts.shape[1]):
                test = pol_inpol_path.contains_points(self.verts[:,iC,:])
                bool_all.append(test)
            test_all = np.array(bool_all)
            cellinlinebox_all_bool = (test_all==True).any(axis=0)
        else:
            raise Exception('ERROR: invalid type for optimize_dist argument')
        
        #intersect_coords = np.empty((0,2,2))
        intersect_coords = np.empty((0,4))
        intersect_gridnos = np.empty((0),dtype=int) #has to be numbers, since a boolean is differently ordered
        verts_inlinebox = self.verts[cellinlinebox_all_bool,:,:]
        verts_inlinebox_nos = np.where(cellinlinebox_all_bool)[0]
        print('finding crossing flow links (can take a while if linebox over xy covers a lot of cells, %i of %i cells are being processed)'%(cellinlinebox_all_bool.sum(),len(cellinlinebox_all_bool)))
        
        for iP, pol_data in enumerate(verts_inlinebox):
            pol_shp = Polygon(pol_data[~np.isnan(pol_data).all(axis=1)])
            intersect_result = pol_shp.intersection(line_section)
            if isinstance(intersect_result,shapely.geometry.multilinestring.MultiLineString): #in the rare case that a cell (pol_shp) is crossed by multiple parts of the line
                intersect_result_multi = intersect_result
            elif isinstance(intersect_result,shapely.geometry.linestring.LineString): #if one linepart trough cell (ex/including node), make multilinestring anyway
                if intersect_result.coords == []: #when the line does not cross this cell, intersect_results.coords is an empty linestring and this cell can be skipped (continue makes forloop continue with next in line without finishing the rest of the steps for this instance)
                    continue
                elif len(intersect_result.coords.xy[0]) == 0: #for newer cartopy versions, when line does not cross this cell, intersect_result.coords.xy is (array('d'), array('d')), and both arrays in tuple have len 0.
                    continue
                intersect_result_multi = MultiLineString([intersect_result])
                
            for iLL, intesect_result_one in enumerate(intersect_result_multi.geoms): #loop over multilinestrings, will mostly only contain one linestring. Will be two if the line crosses a cell more than once.
                intersection_line = intesect_result_one.coords
                intline_xyshape = np.array(intersection_line.xy).shape
                #print('len(intersection_line.xy): %s'%([intline_xyshape]))
                for numlinepart_incell in range(1,intline_xyshape[1]): #is mostly 1, but more if there is a linebreakpoint in this cell (then there are two or more lineparts)
                    intersect_gridnos = np.append(intersect_gridnos,verts_inlinebox_nos[iP])
                    #intersect_coords = np.concatenate([intersect_coords,np.array(intersection_line.xy)[np.newaxis,:,numlinepart_incell-1:numlinepart_incell+1]],axis=0)
                    intersect_coords = np.concatenate([intersect_coords,np.array(intersection_line.xy).T[numlinepart_incell-1:numlinepart_incell+1].flatten()[np.newaxis]])
        
        if intersect_coords.shape[0] != len(intersect_gridnos):
            raise Exception('something went wrong, intersect_coords.shape[0] and len(intersect_gridnos) are not equal')
        
        import pandas as pd
        intersect_pd = pd.DataFrame(intersect_coords,index=intersect_gridnos,columns=['x1','y1','x2','y2'])
        intersect_pd.index.name = 'gridnumber'
        
        #TODO up to here could come from meshkernelpy
        
        print('calculating distance for all crossed cells, from first point of line (should not take long, but if it does, optimisation is needed)')
        nlinecoords = line_array.shape[0]
        nlinedims = len(line_array.shape)
        ncrosscellparts = len(intersect_pd)
        if nlinecoords<2 or nlinedims != 2:
            raise Exception('ERROR: line_array should at least contain two xy points [[x,y],[x,y]]')
        
        #calculate distance between celledge-linepart crossing (is zero when line iL crosses cell)
        distperline_tostart = np.zeros((ncrosscellparts,nlinecoords-1))
        distperline_tostop = np.zeros((ncrosscellparts,nlinecoords-1))
        linepart_length = np.zeros((nlinecoords))
        for iL in range(nlinecoords-1):
            #calculate length of lineparts
            line_section_part = LineString(line_array[iL:iL+2,:])
            if calcdist_fromlatlon:
                linepart_length[iL+1] = calc_dist_haversine(line_array[iL,0],line_array[iL+1,0],line_array[iL,1],line_array[iL+1,1])
            else:
                linepart_length[iL+1] = line_section_part.length
        
            #get distance between all lineparts and point (later used to calculate distance from beginpoint of closest linepart)
            for iP in range(ncrosscellparts):
                distperline_tostart[iP,iL] = line_section_part.distance(Point(intersect_coords[:,0][iP],intersect_coords[:,1][iP]))
                distperline_tostop[iP,iL] = line_section_part.distance(Point(intersect_coords[:,2][iP],intersect_coords[:,3][iP]))
        linepart_lengthcum = np.cumsum(linepart_length)
        cross_points_closestlineid = np.argmin(np.maximum(distperline_tostart,distperline_tostop),axis=1)
        intersect_pd['closestlineid'] = cross_points_closestlineid
        print('finished calculating distance for all crossed cells, from first point of line')
        
        if not calcdist_fromlatlon:
            crs_dist_starts = calc_dist_pythagoras(line_array[cross_points_closestlineid,0], intersect_coords[:,0], line_array[cross_points_closestlineid,1], intersect_coords[:,1]) + linepart_lengthcum[cross_points_closestlineid]
            crs_dist_stops = calc_dist_pythagoras(line_array[cross_points_closestlineid,0], intersect_coords[:,2], line_array[cross_points_closestlineid,1], intersect_coords[:,3]) + linepart_lengthcum[cross_points_closestlineid]
        else:
            crs_dist_starts = calc_dist_haversine(line_array[cross_points_closestlineid,0], intersect_coords[:,0], line_array[cross_points_closestlineid,1], intersect_coords[:,1]) + linepart_lengthcum[cross_points_closestlineid]
            crs_dist_stops = calc_dist_haversine(line_array[cross_points_closestlineid,0], intersect_coords[:,2], line_array[cross_points_closestlineid,1], intersect_coords[:,3]) + linepart_lengthcum[cross_points_closestlineid]
        intersect_pd['crs_dist_starts'] = crs_dist_starts
        intersect_pd['crs_dist_stops'] = crs_dist_stops
        intersect_pd['linepartlen'] = crs_dist_stops-crs_dist_starts
        intersect_pd = intersect_pd.sort_values('crs_dist_starts')
        
        #dimensions (gridnos, xy, firstsecond)
        print('done finding crossing flow links: %i of %i'%(len(intersect_gridnos),len(cellinlinebox_all_bool)))
        return intersect_pd
    



