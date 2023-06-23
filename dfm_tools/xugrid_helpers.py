# -*- coding: utf-8 -*-
"""
Created on Thu Jun 22 17:02:16 2023

@author: veenstra
"""

import numpy as np
import xugrid as xu
import xarray as xr


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


def delft3d4_findnanval(data_nc_XZ,data_nc_YZ):
    values, counts = np.unique(data_nc_XZ, return_counts=True)
    X_nanval = values[np.argmax(counts)]
    values, counts = np.unique(data_nc_YZ, return_counts=True)
    Y_nanval = values[np.argmax(counts)]
    if X_nanval!=Y_nanval:
        XY_nanval = None
    else:
        XY_nanval = X_nanval
    return XY_nanval


def open_dataset_delft3d4(file_nc):
    
    ds = xr.open_dataset(file_nc,chunks={'time':1}) #TODO: move chunks/kwargs to input arguments
    
    #average U1/V1 values to M/N
    for varn in ['U1','V1']:#ds.data_vars:
        var_attrs = ds[varn].attrs
        ds_var = ds[varn]
        if varn=='U1':
            ds_var = ds_var.where(ds.KFU,0)
            ds_var = (ds_var + ds_var.shift(MC=1))/2 #TODO: or MC=-1
            ds_var = ds_var.rename({'MC':'M'})
            var_attrs['long_name'] = var_attrs['long_name'].replace('U-point','zeta point')
        else:
            ds_var = ds_var.where(ds.KFV,0)
            ds_var = (ds_var + ds_var.shift(NC=1))/2 #TODO: or NC=-1
            ds_var = ds_var.rename({'NC':'N'})
            var_attrs['long_name'] = var_attrs['long_name'].replace('V-point','zeta point')
        var_attrs['location'] = 'face'
        ds[varn] = ds_var.assign_attrs(var_attrs)
    
    #compute ux/uy/umag/udir #TODO: add attrs to variables
    ALFAS_rad = np.deg2rad(ds.ALFAS)
    vel_x = ds.U1*np.cos(ALFAS_rad) - ds.V1*np.sin(ALFAS_rad)
    vel_y = ds.U1*np.sin(ALFAS_rad) + ds.V1*np.cos(ALFAS_rad)
    ds['ux'] = vel_x
    ds['uy'] = vel_y
    ds['umag'] = np.sqrt(vel_x**2 + vel_y**2)
    ds['udir'] = np.rad2deg(np.arctan2(vel_y, vel_x))%360
    
    mn_slice = slice(1,None)
    ds = ds.isel(M=mn_slice,N=mn_slice) #cut off first values of M/N (centers), since they are fillvalues and should have different size than MC/NC (corners)
    
    #find and set nans in XZ/YZ arrays, these are ignored in xugrid but still nice to mask
    data_nc_XZ = ds.XZ
    data_nc_YZ = ds.YZ
    XY_nanval = delft3d4_findnanval(data_nc_XZ,data_nc_YZ)
    if XY_nanval is not None:
        mask_XY = (data_nc_XZ==XY_nanval) & (data_nc_YZ==XY_nanval)
        ds['XZ'] = data_nc_XZ.where(~mask_XY)
        ds['YZ'] = data_nc_YZ.where(~mask_XY)

    #find and set nans in XCOR/YCOR arrays
    data_nc_XCOR = ds.XCOR
    data_nc_YCOR = ds.YCOR
    XY_nanval = delft3d4_findnanval(data_nc_XCOR,data_nc_YCOR) #-999.999 in kivu and 0.0 in curvedbend
    if XY_nanval is not None:
        mask_XYCOR = (data_nc_XCOR==XY_nanval) & (data_nc_YCOR==XY_nanval)
        ds['XCOR'] = data_nc_XCOR.where(~mask_XYCOR)
        ds['YCOR'] = data_nc_YCOR.where(~mask_XYCOR)

    #convert to ugrid
    node_coords_x = ds.XCOR.to_numpy().ravel()
    node_coords_y = ds.YCOR.to_numpy().ravel()
    xcor_shape = ds.XCOR.shape
    xcor_nvals = xcor_shape[0] * xcor_shape[1]
    
    #remove weird outlier values in kivu model
    node_coords_x[node_coords_x<-1000] = np.nan
    node_coords_y[node_coords_y<-1000] = np.nan
    
    #find nodes with nan coords
    if not (np.isnan(node_coords_x) == np.isnan(node_coords_y)).all():
        raise Exception('node_coords_xy do not have nans in same location')
    nan_nodes_bool = np.isnan(node_coords_x)
    node_coords_x = node_coords_x[~nan_nodes_bool]
    node_coords_y = node_coords_y[~nan_nodes_bool]
    
    node_idx_square = -np.ones(xcor_nvals,dtype=int)
    node_idx_nonans = np.arange((~nan_nodes_bool).sum())
    node_idx_square[~nan_nodes_bool] = node_idx_nonans
    node_idx = node_idx_square.reshape(xcor_shape)
    face_node_connectivity = np.stack([node_idx[1:,:-1].ravel(), #ll
                                       node_idx[1:,1:].ravel(), #lr
                                       node_idx[:-1,1:].ravel(), #ur
                                       node_idx[:-1,:-1].ravel(), #ul
                                       ],axis=1)
    
    add_triangles = False #seems to be all zeros, so remove from code?
    if add_triangles: # also include faces with 3 nodes?
        faces_w3nodes_bool = (face_node_connectivity!=-1).sum(axis=1)==3
        fnc_3nodes_sorted = np.sort(face_node_connectivity[faces_w3nodes_bool],axis=1)[:,::-1]
        face_node_connectivity[faces_w3nodes_bool] = fnc_3nodes_sorted
        keep_faces_bool = (face_node_connectivity!=-1).sum(axis=1)>=3
    else:
        keep_faces_bool = (face_node_connectivity!=-1).sum(axis=1)==4
    
    face_node_connectivity = face_node_connectivity[keep_faces_bool]
    
    grid = xu.Ugrid2d(node_x=node_coords_x,
                      node_y=node_coords_y,
                      face_node_connectivity=face_node_connectivity,
                      fill_value=-1,
                      )
    
    face_dim = grid.face_dimension
    ds_stacked = ds.stack({face_dim:('M','N')}).sel({face_dim:keep_faces_bool})
    ds_stacked = ds_stacked.drop_vars(['M','N','mesh2d_nFaces'])
    uds = xu.UgridDataset(ds_stacked,grids=[grid]) 
    
    uds = uds.drop_vars(['XCOR','YCOR'])#,'KCU','KCV','KFU','KFV','DP0','DPU0','DPV0']) #TODO: #drop additional vars with MC/NC (automate)
    
    return uds
