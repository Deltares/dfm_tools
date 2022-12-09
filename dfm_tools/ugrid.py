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
along with this program.  if not, see <http://www.gnu.org/licenses/>.

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

import warnings


class UGrid: #TODO: remove UGrid class, not needed anymore with xugrid
    """Unstructured grid"""
    def __init__(self, mesh2d_node_x, mesh2d_node_y, mesh2d_face_nodes, verts, mesh2d_node_z=None, edge_verts=None):
        warnings.warn(DeprecationWarning('dfm_tools.ugrid.UGrid() will be deprecated, since there is an xarray alternative for multidomain FM files (xugrid). Open your file like this and use xarray sel/isel (example in postprocessing notebook):\n    data_xr_mapmerged = dfmt.open_partitioned_dataset(file_nc_map).\nIf you need the verts for some reason, use: ugrid_all_verts = dfmt.get_ugrid_verts(data_frommap_merged)'))
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


    def polygon_intersect(self, line_array, optimize_dist=False, calcdist_fromlatlon=False): #TODO: remove this code
        raise DeprecationWarning('ugrid.polygon_intersect() is deprecated. Cross sections are now computed with xugrid/xarray objects as input like this (example in postprocessing notebook):\n    data_xr_mapmerged = dfmt.open_partitioned_dataset(file_nc_map)\n    intersect_pd = dfmt.polygon_intersect(data_xr_mapmerged, line_array)\n    crs_verts, crs_plotdata = dfmt.get_xzcoords_onintersection(data_xr_mapmerged, varname="mesh2d_sa1", intersect_pd=intersect_pd, timestep=3)\n    fig, ax = plt.subplots()\n    pc = dfmt.plot_netmapdata(crs_verts, values=crs_plotdata, ax=ax)')
        
    
