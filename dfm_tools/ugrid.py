# -*- coding: utf-8 -*-
"""
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
        from dfm_tools.get_nc_helpers import get_varname_mapnc, ghostcell_filter
        
        def nodexyfaces2verts(node_x,node_y, faces):
            quatrangles = faces-1 #convert 1-based indexing of cell numbering in ugrid to 0-based indexing
            #https://stackoverflow.com/questions/49640311/matplotlib-unstructered-quadrilaterals-instead-of-triangles
            #https://stackoverflow.com/questions/52202014/how-can-i-plot-2d-fem-results-using-matplotlib
            yz = np.c_[node_x,node_y]
            verts= yz[quatrangles]
            verts[quatrangles.mask==True,:] = np.nan #remove all masked values by making them nan
            return verts
        
        data_nc = Dataset(file_nc)

        mesh2d_node_x = data_nc.variables[get_varname_mapnc(data_nc,'mesh2d_node_x')][:]
        mesh2d_node_y = data_nc.variables[get_varname_mapnc(data_nc,'mesh2d_node_y')][:]
        varn_mesh2d_node_z = get_varname_mapnc(data_nc,'mesh2d_node_z')
        if varn_mesh2d_node_z is not None: # node_z variable is present
            mesh2d_node_z = data_nc.variables[varn_mesh2d_node_z][:]
        else:
            mesh2d_node_z = None
        varn_mesh2d_face_nodes = get_varname_mapnc(data_nc,'mesh2d_face_nodes')
        if varn_mesh2d_face_nodes is not None: # node_z variable is present
            mesh2d_face_nodes = data_nc.variables[varn_mesh2d_face_nodes][:, :]
        else:
            raise Exception('ERROR: provided file does not contain a variable mesh2d_face_nodes or similar:\n%s\nPlease do one of the following:\n- plot grid from *_map.nc file\n- import and export the grid with RGFGRID\n- import and save the gridd "with cellfinfo" from interacter'%(file_nc))
        verts = nodexyfaces2verts(mesh2d_node_x, mesh2d_node_y, mesh2d_face_nodes) #xy coordinates of face nodes
        
        varn_mesh2d_edge_x = get_varname_mapnc(data_nc,'mesh2d_edge_x')
        if varn_mesh2d_edge_x is not None: # mesh2d_edge_x (and mesh2d_edge_y) variable is present
            mesh2d_edge_x = data_nc.variables[varn_mesh2d_edge_x][:]
            mesh2d_edge_y = data_nc.variables[get_varname_mapnc(data_nc,'mesh2d_edge_y')][:]
            mesh2d_edge_nodes = data_nc.variables[get_varname_mapnc(data_nc,'mesh2d_edge_nodes')][:]
            edge_verts = nodexyfaces2verts(mesh2d_edge_x, mesh2d_edge_y, mesh2d_edge_nodes) #xy coordinates of face nodes
        else:
            edge_verts = None
            
        #remove ghost cells from faces and verts
        ghostcells_bool, nonghost_ids = ghostcell_filter(file_nc)
        if ghostcells_bool:
            mesh2d_face_nodes = mesh2d_face_nodes[nonghost_ids]
            verts = verts[nonghost_ids]
        
        data_nc.close()
        ugrid = UGrid(mesh2d_node_x, mesh2d_node_y, mesh2d_face_nodes, verts, mesh2d_node_z=mesh2d_node_z, edge_verts=edge_verts)
        return ugrid

    def polygon_intersect(self, line_array):
        import numpy as np
        try:
            from shapely.geometry import Polygon, LineString
        except:
            raise Exception('ERROR: cannot ')
            
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
        
        cellinlinebox_all_bool = (((np.min(line_array[:,0]) < verts_xmax) &
                                   (np.max(line_array[:,0]) > verts_xmin)) &
                                  ((np.min(line_array[:,1]) < verts_ymax) & 
                                   (np.max(line_array[:,1]) > verts_ymin))
                                  )
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
    
