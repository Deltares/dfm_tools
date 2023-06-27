# -*- coding: utf-8 -*-
"""
Created on Thu Jun 22 17:02:16 2023

@author: veenstra
"""

import numpy as np
import xugrid as xu
import xarray as xr
import datetime as dt
import pandas as pd
from dfm_tools.xarray_helpers import file_to_list


def get_vertical_dimensions(uds): #TODO: maybe add layer_dimension and interface_dimension properties to xugrid?
    """
    get vertical_dimensions from grid_info of ugrid mapfile (this will fail for hisfiles). The info is stored in the layer_dimension and interface_dimension attribute of the mesh2d variable of the dataset (stored in uds.grid after reading with xugrid)
    
    processing cb_3d_map.nc
        >> found layer/interface dimensions in file: mesh2d_nLayers mesh2d_nInterfaces
    processing Grevelingen-FM_0*_map.nc
        >> found layer/interface dimensions in file: nmesh2d_layer nmesh2d_interface (these are updated in open_partitioned_dataset)
    processing DCSM-FM_0_5nm_0*_map.nc
        >> found layer/interface dimensions in file: mesh2d_nLayers mesh2d_nInterfaces
    processing MB_02_0*_map.nc
        >> found layer/interface dimensions in file: mesh2d_nLayers mesh2d_nInterfaces
    """
    
    if not hasattr(uds,'grid'): #early return in case of e.g. hisfile
        return None, None
        
    gridname = uds.grid.name
    grid_info = uds.grid.to_dataset()[gridname]
    if hasattr(grid_info,'layer_dimension'):
        return grid_info.layer_dimension, grid_info.interface_dimension
    else:
        return None, None


def remove_ghostcells(uds): #TODO: create JIRA issue: remove ghostcells from output (or make values in ghostcells the same as not-ghostcells, is now not the case for velocities, probably due to edge-to-center interpolation)
    """
    Dropping ghostcells if there is a domainno variable present and there is a domainno in the filename.
    Not using most-occurring domainno in var, since this is not a valid assumption for merged datasets and might be invalid for a very small partition.
    
    """
    gridname = uds.grid.name
    varn_domain = f'{gridname}_flowelem_domain'
    
    #check if dataset has domainno variable, return uds if not present
    if varn_domain not in uds.data_vars:
        print('[nodomainvar] ',end='')
        return uds
    
    #derive domainno from filename, return uds if not present
    fname = uds.encoding['source']
    if '_' not in fname: #safety escape in case there is no _ in the filename
        print('[nodomainfname] ',end='')
        return uds
    fname_splitted = fname.split('_')
    part_domainno_fromfname = fname_splitted[-2] #this is not valid for rstfiles (date follows after partnumber), but they cannot be read with xugrid anyway since they are mapformat=1
    if not part_domainno_fromfname.isnumeric() or len(part_domainno_fromfname)!=4:
        print('[nodomainfname] ',end='')
        return uds
    
    #drop ghostcells
    part_domainno_fromfname = int(part_domainno_fromfname)
    da_domainno = uds[varn_domain]
    idx = np.flatnonzero(da_domainno == part_domainno_fromfname)
    uds = uds.isel({uds.grid.face_dimension:idx})
    return uds


def remove_periodic_cells(uds): #TODO: implement proper fix: https://github.com/Deltares/xugrid/issues/63
    """
    For global models with grids that go "around the back". Temporary fix to drop all faces that are larger than grid_extent/2 (eg 360/2=180 degrees in case of GTSM)
    
    """
    face_node_x = uds.grid.face_node_coordinates[:,:,0]
    grid_extent = uds.grid.bounds[2] - uds.grid.bounds[0]
    face_node_maxdx = np.nanmax(face_node_x,axis=1) - np.nanmin(face_node_x,axis=1)
    bool_face = face_node_maxdx < grid_extent/2
    if bool_face.all(): #early return for when no cells have to be removed (might increase performance)
        return uds
    print(f'>> removing {(~bool_face).sum()} periodic cells from dataset: ',end='')
    dtstart = dt.datetime.now()
    uds = uds.sel({uds.grid.face_dimension:bool_face})
    print(f'{(dt.datetime.now()-dtstart).total_seconds():.2f} sec')
    return uds


def open_partitioned_dataset(file_nc, remove_ghost=True, **kwargs): 
    """
    using xugrid to read and merge partitions, with some additional features (remaning old layerdim, timings, set zcc/zw as data_vars)

    Parameters
    ----------
    file_nc : TYPE
        DESCRIPTION.
    kawrgs : TYPE, optional
        arguments that are passed to xr.open_dataset. The chunks argument is set if not provided
        chunks={'time':1} increases performance significantly upon reading, but causes memory overloads when performing sum/mean/etc actions over time dimension (in that case 100/200 is better). The default is {'time':1}.

    Raises
    ------
    Exception
        DESCRIPTION.

    Returns
    -------
    ds_merged_xu : TYPE
        DESCRIPTION.
    
    file_nc = 'p:\\1204257-dcsmzuno\\2006-2012\\3D-DCSM-FM\\A18b_ntsu1\\DFM_OUTPUT_DCSM-FM_0_5nm\\DCSM-FM_0_5nm_0*_map.nc' #3D DCSM
    file_nc = 'p:\\archivedprojects\\11206813-006-kpp2021_rmm-2d\\C_Work\\31_RMM_FMmodel\\computations\\model_setup\\run_207\\results\\RMM_dflowfm_0*_map.nc' #RMM 2D
    file_nc = 'p:\\1230882-emodnet_hrsm\\GTSMv5.0\\runs\\reference_GTSMv4.1_wiCA_2.20.06_mapformat4\\output\\gtsm_model_0*_map.nc' #GTSM 2D
    file_nc = 'p:\\11208053-005-kpp2022-rmm3d\\C_Work\\01_saltiMarlein\\RMM_2019_computations_02\\computations\\theo_03\\DFM_OUTPUT_RMM_dflowfm_2019\\RMM_dflowfm_2019_0*_map.nc' #RMM 3D
    file_nc = 'p:\\archivedprojects\\11203379-005-mwra-updated-bem\\03_model\\02_final\\A72_ntsu0_kzlb2\\DFM_OUTPUT_MB_02\\MB_02_0*_map.nc'
    Timings (xu.open_dataset/xu.merge_partitions):
        - DCSM 3D 20 partitions  367 timesteps: 231.5/ 4.5 sec (decode_times=False: 229.0 sec)
        - RMM  2D  8 partitions  421 timesteps:  55.4/ 4.4 sec (decode_times=False:  56.6 sec)
        - GTSM 2D  8 partitions  746 timesteps:  71.8/30.0 sec (decode_times=False: 204.8 sec)
        - RMM  3D 40 partitions  146 timesteps: 168.8/ 6.3 sec (decode_times=False: 158.4 sec)
        - MWRA 3D 20 partitions 2551 timesteps:  74.4/ 3.4 sec (decode_times=False:  79.0 sec)
    
    """
    #TODO: FM-mapfiles contain wgs84/projected_coordinate_system variables. xugrid has .crs property, projected_coordinate_system/wgs84 should be updated to be crs so it will be automatically handled? >> make dflowfm issue (and https://github.com/Deltares/xugrid/issues/42)
    #TODO: add support for multiple grids via keyword? GTSM+riv grid also only contains only one grid, so no testcase available
    #TODO: speed up open_dataset https://github.com/Deltares/dfm_tools/issues/225 (also remove_ghost)
    
    if 'chunks' not in kwargs:
        kwargs['chunks'] = {'time':1}
    
    dtstart_all = dt.datetime.now()
    file_nc_list = file_to_list(file_nc)
    
    print(f'>> xu.open_dataset() with {len(file_nc_list)} partition(s): ',end='')
    dtstart = dt.datetime.now()
    partitions = []
    for iF, file_nc_one in enumerate(file_nc_list):
        print(iF+1,end=' ')
        ds = xr.open_dataset(file_nc_one, **kwargs)
        if 'nFlowElem' in ds.dims and 'nNetElem' in ds.dims: #for mapformat1 mapfiles: merge different face dimensions (rename nFlowElem to nNetElem) to make sure the dataset topology is correct
            print('[mapformat1] ',end='')
            ds = ds.rename({'nFlowElem':'nNetElem'})
        uds = xu.core.wrap.UgridDataset(ds)
        if remove_ghost: #TODO: this makes it way slower (at least for GTSM, although merging seems faster), but is necessary since values on overlapping cells are not always identical (eg in case of Venice ucmag)
            uds = remove_ghostcells(uds)
        partitions.append(uds)
    print(': ',end='')
    print(f'{(dt.datetime.now()-dtstart).total_seconds():.2f} sec')
    
    if len(partitions) == 1: #do not merge in case of 1 partition
        return partitions[0]
    
    print(f'>> xu.merge_partitions() with {len(file_nc_list)} partition(s): ',end='')
    dtstart = dt.datetime.now()
    ds_merged_xu = xu.merge_partitions(partitions)
    print(f'{(dt.datetime.now()-dtstart).total_seconds():.2f} sec')
    
    #print variables that are dropped in merging procedure. Often only ['mesh2d_face_x_bnd', 'mesh2d_face_y_bnd'], which can be derived by combining node_coordinates (mesh2d_node_x mesh2d_node_y) and face_node_connectivity (mesh2d_face_nodes). >> can be removed from FM-mapfiles (email of 16-1-2023)
    varlist_onepart = list(partitions[0].variables.keys())
    varlist_merged = list(ds_merged_xu.variables.keys())
    varlist_dropped_bool = ~pd.Series(varlist_onepart).isin(varlist_merged)
    varlist_dropped = pd.Series(varlist_onepart).loc[varlist_dropped_bool]
    if varlist_dropped_bool.any():
        print(f'>> some variables dropped with merging of partitions: {varlist_dropped.tolist()}')
    
    print(f'>> dfmt.open_partitioned_dataset() total: {(dt.datetime.now()-dtstart_all).total_seconds():.2f} sec')
    return ds_merged_xu


def open_dataset_curvilinear(file_nc,
                             varn_vert_lon='vertices_longitude', #'grid_x'
                             varn_vert_lat='vertices_latitude', #'grid_y'
                             ij_dims=['i','j'], #['N','M']
                             convert_360to180=False,
                             **kwargs):
    """
    This is a first version of a function that creates a xugrid UgridDataset from a curvilinear dataset like CMCC. Curvilinear means in this case 2D lat/lon variables and i/j indexing. The CMCC dataset does contain vertices, which is essential for conversion to ugrid.
    It also works for WAQUA files that are converted with getdata
    """
    
    if 'chunks' not in kwargs:
        kwargs['chunks'] = {'time':1}
    
    ds = xr.open_mfdataset(file_nc, **kwargs)
    
    print('>> getting vertices from ds: ',end='') #long (but does not reflect in 
    dtstart = dt.datetime.now()
    vertices_longitude = ds[varn_vert_lon].to_numpy()
    vertices_longitude = vertices_longitude.reshape(-1,vertices_longitude.shape[-1])
    vertices_latitude = ds[varn_vert_lat].to_numpy()
    vertices_latitude = vertices_latitude.reshape(-1,vertices_latitude.shape[-1])
    print(f'{(dt.datetime.now()-dtstart).total_seconds():.2f} sec')
    
    #convert from 0to360 to -180 to 180
    if convert_360to180:
        vertices_longitude = (vertices_longitude+180)%360 - 180 #TODO: check if periodic cell filter still works properly after doing this
    
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
    
    #remove all faces that have only 1 unique node (does not result in a valid grid) #TODO: not used yet except for print
    fnc_all_duplicates = (face_node_connectivity.T==face_node_connectivity[:,0]).all(axis=0)
    
    #create bool of cells with duplicate nodes (some have 1 unique node, some 3, all these are dropped) #TODO: support also triangles?
    fnc_closed = np.c_[face_node_connectivity,face_node_connectivity[:,0]]
    fnc_has_duplicates = (np.diff(fnc_closed,axis=1)==0).any(axis=1)
    
    #only keep cells that have 4 unique nodes
    bool_combined = ~fnc_has_duplicates
    print(f'WARNING: dropping {fnc_has_duplicates.sum()} faces with duplicate nodes ({fnc_all_duplicates.sum()} with one unique node)')
    face_node_connectivity = face_node_connectivity[bool_combined]
    
    grid = xu.Ugrid2d(node_x=node_coords_x,
                      node_y=node_coords_y,
                      face_node_connectivity=face_node_connectivity,
                      fill_value=-1,
                      )
    
    print('>> stacking ds i/j coordinates: ',end='') #fast
    dtstart = dt.datetime.now()
    face_dim = grid.face_dimension
    ds_stacked = ds.stack({face_dim:ij_dims}).sel({face_dim:bool_combined}) #TODO: lev/time bnds are dropped, avoid this. maybe stack initial dataset since it would also simplify the rest of the function a bit
    ds_stacked = ds_stacked.drop_vars(ij_dims+['mesh2d_nFaces']) #TODO: solve "DeprecationWarning: Deleting a single level of a MultiIndex is deprecated", solved by removing mesh2d_nFaces variable
    print(f'{(dt.datetime.now()-dtstart).total_seconds():.2f} sec')
    
    print('>> init uds: ',end='') #long
    dtstart = dt.datetime.now()
    uds = xu.UgridDataset(ds_stacked,grids=[grid])
    print(f'{(dt.datetime.now()-dtstart).total_seconds():.2f} sec')
    
    #remove faces that link to node coordinates that are nan (occurs in waqua models)
    bool_faces_wnannodes = np.isnan(uds.grid.face_node_coordinates[:,:,0]).any(axis=1)
    if bool_faces_wnannodes.any():
        print(f'>> drop {bool_faces_wnannodes.sum()} faces with nan nodecoordinates from uds: ',end='') #long
        dtstart = dt.datetime.now()
        uds = uds.sel({face_dim:~bool_faces_wnannodes})
        print(f'{(dt.datetime.now()-dtstart).total_seconds():.2f} sec')
    
    return uds


def get_delft3d4_nanmask(x,y):
    bool_0 = (x==0) & (y==0)
    bool_1 = (x==-999) & (y==-999)
    bool_2 = (x==-999.999) & (y==-999.999)
    bool_mask = bool_0 | bool_1 | bool_2
    return bool_mask


def open_dataset_delft3d4(file_nc, **kwargs):
    
    if 'chunks' not in kwargs:
        kwargs['chunks'] = {'time':1}
    
    ds = xr.open_dataset(file_nc, **kwargs)
    
    if ('U1' in ds.data_vars) and ('V1' in ds.data_vars):
        #mask u and v separately with 0 to avoid high velocities (cannot be nan, since (nan+value)/2= nan instead of value=2
        mask_u1 = (ds.U1==-999) | (ds.U1==-999.999)
        mask_v1 = (ds.V1==-999) | (ds.V1==-999.999)
        u1_mn = ds.U1.where(~mask_u1,0)
        v1_mn = ds.V1.where(~mask_v1,0)
        
        #create combined uv mask (have to rename dimensions)
        mask_u1_mn = mask_u1.rename({'MC':'M'})
        mask_v1_mn = mask_v1.rename({'NC':'N'})
        mask_uv1_mn = mask_u1_mn & mask_v1_mn

        #average U1/V1 values to M/N
        u1_mn = (u1_mn + u1_mn.shift(MC=1))/2 #TODO: or MC=-1
        u1_mn = u1_mn.rename({'MC':'M'})
        u1_mn = u1_mn.where(~mask_uv1_mn,np.nan) #replace temporary zeros with nan
        v1_mn = (v1_mn + v1_mn.shift(NC=1))/2 #TODO: or NC=-1
        v1_mn = v1_mn.rename({'NC':'N'})
        v1_mn = v1_mn.where(~mask_uv1_mn,np.nan) #replace temporary zeros with nan
        ds = ds.drop_vars(['U1','V1']) #to avoid creating large chunks, alternative is to overwrite the vars with the MN-averaged vars, but it requires passing and updating of attrs
        
        #compute ux/uy/umag/udir #TODO: add attrs to variables
        alfas_rad = np.deg2rad(ds.ALFAS)
        vel_x = u1_mn*np.cos(alfas_rad) - v1_mn*np.sin(alfas_rad)
        vel_y = u1_mn*np.sin(alfas_rad) + v1_mn*np.cos(alfas_rad)
        ds['ux'] = vel_x
        ds['uy'] = vel_y
        ds['umag'] = np.sqrt(vel_x**2 + vel_y**2)
        ds['udir'] = np.rad2deg(np.arctan2(vel_y, vel_x))%360
    
    mn_slice = slice(1,None)
    ds = ds.isel(M=mn_slice,N=mn_slice) #cut off first values of M/N (centers), since they are fillvalues and should have different size than MC/NC (corners)
    
    #find and set nans in XCOR/YCOR arrays
    mask_xy = get_delft3d4_nanmask(ds.XCOR,ds.YCOR) #-999.999 in kivu and 0.0 in curvedbend, both in westernscheldt
    ds['XCOR'] = ds.XCOR.where(~mask_xy)
    ds['YCOR'] = ds.YCOR.where(~mask_xy)

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
    
    uds = uds.drop_vars(['XCOR','YCOR'])
    uds = uds.drop_dims(['MC','NC']) #clean up dataset by dropping corner dims (drops also variabes with U/V masks and U/V/C bedlevel)
    
    return uds
