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
            #https://stackoverflow.com/questions/49640311/matplotlib-unstructered-quadrilaterals-instead-of-triangles
            #https://stackoverflow.com/questions/52202014/how-can-i-plot-2d-fem-results-using-matplotlib
            yz = np.c_[node_x,node_y]
            verts= yz[quatrangles]
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
        nonghost_ids = ghostcell_filter(file_nc)
        if nonghost_ids is not None:
            mesh2d_face_nodes = mesh2d_face_nodes[nonghost_ids]
            verts = verts[nonghost_ids]
        
        data_nc.close()
        ugrid = UGrid(mesh2d_node_x, mesh2d_node_y, mesh2d_face_nodes, verts, mesh2d_node_z=mesh2d_node_z, edge_verts=edge_verts)
        return ugrid

    def polygon_intersect(self, line_array, optimize_dist=None):
        import numpy as np
        from matplotlib.path import Path
        
        from dfm_tools.testutils import try_importmodule
        try_importmodule(modulename='shapely')
        from shapely.geometry import LineString, Polygon
            
        print('finding crossing flow links (can take a while if linebox over xy covers a lot of cells)')
        #allpol = []
        intersect_gridnos = np.empty((0),dtype=int)
        intersect_coords1 = np.empty((0,2))
        intersect_coords2 = np.empty((0,2))
        line_section = LineString(line_array)
        
        verts_xmax = np.nanmax(self.verts[:,:,0].data,axis=1)
        verts_xmin = np.nanmin(self.verts[:,:,0].data,axis=1)
        verts_ymax = np.nanmax(self.verts[:,:,1].data,axis=1)
        verts_ymin = np.nanmin(self.verts[:,:,1].data,axis=1)
        
        if optimize_dist is None:
            cellinlinebox_all_bool = (((np.min(line_array[:,0]) < verts_xmax) &
                                       (np.max(line_array[:,0]) > verts_xmin)) &
                                      ((np.min(line_array[:,1]) < verts_ymax) & 
                                       (np.max(line_array[:,1]) > verts_ymin))
                                      )
        elif type(optimize_dist) is int or type(optimize_dist) is float:
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
        
        
        verts_inlinebox = self.verts[cellinlinebox_all_bool,:,:]
        verts_inlinebox_nos = np.where(cellinlinebox_all_bool)[0]
        
        for iP, pol_data in enumerate(verts_inlinebox):
            pol_data_nonan = pol_data[~np.isnan(pol_data).all(axis=1)]
            pol_shp = Polygon(pol_data_nonan)
            try:
                intersection_line = pol_shp.intersection(line_section).coords
            except: #in the rare case that a cell (pol_shp) is crossed by multiple parts of the line
                intersection_line = pol_shp.intersection(line_section)[0].coords
            
            if intersection_line != []:
                intersect_gridnos = np.concatenate([intersect_gridnos,[verts_inlinebox_nos[iP]]])
                intersect_coords1=np.concatenate([intersect_coords1,np.array([list(intersection_line)[0]])])
                intersect_coords2=np.concatenate([intersect_coords2,np.array([list(intersection_line)[1]])])
                #all_intersect.append(list(intersection_line))
            #print(iP)

        intersect_coords = np.stack([intersect_coords1,intersect_coords2], axis=2)
        #dimensions (gridnos, xy, firstsecond)
        #allpol_multi = MultiPolygon(allpol)
        #intersection_line = allpol_multi.intersection(line_section).coords #does not work
        print('done finding crossing flow links')
        return intersect_gridnos, intersect_coords
    



