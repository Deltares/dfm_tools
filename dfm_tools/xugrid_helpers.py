# -*- coding: utf-8 -*-
"""
Created on Thu Jun 22 17:02:16 2023

@author: veenstra
"""

import numpy as np
import xugrid as xu


def curvilinear_to_UgridDataset(ds):
    """
    This is a first version of a function that creates a xugrid UgridDataset from a curvilinear dataset like CMCC. Curvilinear means in this case 2D lat/lon variables and i/j indexing. The CMCC dataset does contain vertices, which is essential for conversion to ugrid.
    """
    vertices_longitude = ds.vertices_longitude.to_numpy()
    vertices_longitude = vertices_longitude.reshape(-1,vertices_longitude.shape[-1])
    vertices_latitude = ds.vertices_latitude.to_numpy()
    vertices_latitude = vertices_latitude.reshape(-1,vertices_latitude.shape[-1])

    #convert from 0to360 to -180 to 180
    vertices_longitude = (vertices_longitude+180)%360 - 180 #TODO: check if periodic cell filter still works properly
    
    # face_xy = np.stack([longitude,latitude],axis=-1)
    # face_coords_x, face_coords_y = face_xy.T
    #a,b = np.unique(face_xy,axis=0,return_index=True) #TODO: there are non_unique face_xy values, inconvenient
    face_xy_vertices = np.stack([vertices_longitude,vertices_latitude],axis=-1)
    face_xy_vertices_flat = face_xy_vertices.reshape(-1,2)
    uniq,inv = np.unique(face_xy_vertices_flat, axis=0, return_inverse=True)
    #len(uniq) = 104926 >> amount of unique node coords
    #uniq.max() = 359.9654541015625 >> node_coords_xy
    #len(inv) = 422816 >> is length of face_xy_vertices.reshape(-1,2)
    #inv.max() = 104925 >> node numbers
    node_coords_x, node_coords_y = uniq.T
    
    face_node_connectivity = inv.reshape(face_xy_vertices.shape[:2]) #fnc.max() = 104925
    
    #remove all faces that have only 1 unique node (does not result in a valid grid) #TODO: not used yet
    fnc_all_duplicates = (face_node_connectivity.T==face_node_connectivity[:,0]).all(axis=0)
    
    #create bool of cells with duplicate nodes (some have 1 unique node, some 3, all these are dropped) #TODO: support also triangles
    fnc_closed = np.c_[face_node_connectivity,face_node_connectivity[:,0]]
    fnc_has_duplicates = (np.diff(fnc_closed,axis=1)==0).any(axis=1)
    
    #create boolean of perodic cells (going around the back) #TODO: put function in xugrid to convert periodic to non-periodic grid: https://github.com/Deltares/xugrid/issues/63
    face_node_coordinates = uniq[face_node_connectivity]
    fn_coords_diff = (face_node_coordinates[...,0].max(axis=1)-face_node_coordinates[...,0].min(axis=1))
    bool_periodic_cells = (fn_coords_diff>180)
    
    #only keep cells that are not periodic and have 4 unique nodes
    bool_combined = ~fnc_has_duplicates & ~bool_periodic_cells
    print(f'WARNING: dropping {fnc_has_duplicates.sum()} faces with duplicate nodes ({fnc_all_duplicates.sum()} with one unique node), dropping {bool_periodic_cells.sum()} periodic cells')
    face_node_connectivity = face_node_connectivity[bool_combined]
    
    grid = xu.Ugrid2d(node_x=node_coords_x,
                      node_y=node_coords_y,
                      face_node_connectivity=face_node_connectivity,
                      fill_value=-1,
                      )
    # fig, ax = plt.subplots()
    # grid.plot(ax=ax)
    
    face_dim = grid.face_dimension
    ds_stacked = ds.stack({face_dim:('i','j')}).sel({face_dim:bool_combined}) #TODO: lev/time bnds are dropped, avoid this. maybe stack initial dataset since it would also simplify the rest of the function a bit
    ds_stacked = ds_stacked.drop_vars(['i','j','mesh2d_nFaces']) #TODO: solve "DeprecationWarning: Deleting a single level of a MultiIndex is deprecated", solved by removing mesh2d_nFaces variable
    uds = xu.UgridDataset(ds_stacked,grids=[grid])
    return uds

