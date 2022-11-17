# -*- coding: utf-8 -*-
"""
Created on Mon Sep 26 13:41:09 2022

@author: veenstra
"""

import os
import numpy as np
import xugrid as xu
import xarray as xr
import matplotlib.pyplot as plt
plt.close('all')
import datetime as dt
import dfm_tools as dfmt
import glob
import pandas as pd

dtstart_all = dt.datetime.now()

#dir_model = r"c:\tmp\xugrid\Grevelingen_run01"
dir_testdata = r'c:\DATA\dfm_tools_testdata'


if 0: #REFERENCE DFM_TOOLS
    
    file_nc = os.path.join(dir_testdata,r'DFM_3D_z_Grevelingen\computations\run01','DFM_OUTPUT_Grevelingen-FM','Grevelingen-FM_0000_map.nc')
    
    timestep = 3
    layno = 33 #35 is top
    
    #get ugrid data
    ugrid_all = dfmt.get_netdata(file_nc=file_nc)
    
    print('plot grid and values from mapdata (salinity on layer, 3dim, on cell centers)')
    data_frommap = dfmt.get_ncmodeldata(file_nc=file_nc, varname='mesh2d_s1', timestep=timestep)#, layer=layno)
    fig, ax = plt.subplots()
    pc = dfmt.plot_netmapdata(ugrid_all.verts, values=data_frommap, ax=None, linewidth=0.5, cmap="jet", edgecolor='face')
    cbar = fig.colorbar(pc, ax=ax)
    cbar.set_label('%s [%s]'%(data_frommap.var_ncattrs['long_name'], data_frommap.var_ncattrs['units']))
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    


if 1: #method huite, adjusted
    
    def open_partitioned_dataset(file_nc, only_faces=False, chunks={'time':1}): #chunks={'time':1} increases performance significantly
        
        if isinstance(file_nc,list):
            file_nc_list = file_nc
        else:
            file_nc_list = glob.glob(file_nc)
        if len(file_nc_list)==0:
            raise Exception('file(s) not found, empty file_nc_list')
        
        partitions = [xu.open_dataset(file_nc_one,chunks=chunks) for file_nc_one in file_nc_list]
        if 'idomain' in partitions[0].data_vars: #if netfile, rename to mapfile varnames
            rename_dict = {'idomain':'mesh2d_flowelem_domain',
                           'iglobal_s':'mesh2d_flowelem_globalnr',
                           'nNetElemMaxNode':'max_nmesh2d_face_nodes'} #TODO: extend with others
            partitions = [part.rename(rename_dict) for part in partitions]
        
        all_indices = []
        all_faces = []
        all_nodes_x = []
        all_nodes_y = []
        accumulator = 0
        for i, part in enumerate(partitions):
            # For ghost nodes, keep the values of the lower numbered partition.
            grid = part.ugrid.grid
            idx = np.flatnonzero(part['mesh2d_flowelem_domain'] >= i) #TODO: wonder if domains are now correct, since I see different domain numbers in the domain plot of 1 partition
            faces = grid.face_node_connectivity[idx]
            faces[faces != grid.fill_value] += accumulator
            accumulator += grid.n_node
            all_indices.append(idx)
            all_faces.append(faces)
            all_nodes_x.append(grid.node_x)
            all_nodes_y.append(grid.node_y)
        node_x = np.concatenate(all_nodes_x)
        node_y = np.concatenate(all_nodes_y)
        node_xy = np.column_stack([node_x, node_y])
        merged_nodes, inverse = np.unique(node_xy, return_inverse=True, axis=0)
        n_face_total = sum(len(faces) for faces in all_faces)
        n_max_node = max(faces.shape[1] for faces in all_faces)
        merged_faces = np.full((n_face_total, n_max_node), -1, dtype=np.intp)
        start = 0 #TODO: use the globalnumbers instead?
        for faces in all_faces:
            n_face, n_max_node = faces.shape
            end = start + n_face
            merged_faces[start:end, :n_max_node] = faces
            start = end
        isnode = merged_faces != -1
        faces_flat = merged_faces[isnode]
        renumbered = inverse[faces_flat]
        merged_faces[isnode] = renumbered
        merged_grid = xu.Ugrid2d(
            node_x=merged_nodes[:, 0],
            node_y=merged_nodes[:, 1],
            fill_value=-1,
            face_node_connectivity=merged_faces,
        )
        facedim = partitions[0].ugrid.grid.face_dimension
        nodedim = partitions[0].ugrid.grid.node_dimension
        edgedim = partitions[0].ugrid.grid.edge_dimension
        #define list of variables per dimension

        ds_face_list = []
        ds_node_list = []
        ds_edge_list = []
        ds_rest_list = []
        dtstart = dt.datetime.now()
        for idx, uds in zip(all_indices, partitions):
            
            #TODO: commented part does not work yet since we need exceptions for eg 'max_nmesh2d_face_nodes'
            # uds_face = dfmt.Dataset_varswithdim(uds,dimname=facedim)
            # uds_node = dfmt.Dataset_varswithdim(uds,dimname=nodedim)
            # uds_edge = dfmt.Dataset_varswithdim(uds,dimname=edgedim)
            face_variables = []        
            node_variables = []        
            edge_variables = []        
            for varname in uds.variables.keys(): #TODO: do this only for one partition (or is that unsafe?)
                if 'max_nmesh2d_face_nodes' in uds[varname].dims: #TODO: not possible to concatenate this dim yet (size varies per partition) #therefore, vars mesh2d_face_x_bnd and mesh2d_face_y_bnd cannot be included currently
                    continue
                if facedim in uds[varname].dims:
                    face_variables.append(varname)
                if nodedim in uds[varname].dims:
                    node_variables.append(varname)
                if edgedim in uds[varname].dims:
                    edge_variables.append(varname)
            ds_face = uds.ugrid.obj[face_variables] #TODO: why/how does this work?
            if not only_faces:
                ds_node = uds.ugrid.obj[node_variables]
                ds_edge = uds.ugrid.obj[edge_variables]
                ds_rest = uds.drop_dims([facedim,nodedim,edgedim]) #TODO: this is a ugrid dataset, cannot be concatenated (drop_dims is not available under uds.ugrid) #contains 4/6 dropped variables
            
            ds_face_list.append(ds_face.isel({facedim: idx}))
            if not only_faces:
                ds_node_list.append(ds_node)#.isel({nodedim: idx})) #TODO: add ghostcell removal for nodes and edges
                ds_edge_list.append(ds_edge)#.isel({edgedim: idx}))
                ds_rest_list.append(ds_rest)#.isel({edgedim: idx}))
        print(f'>>time .sel/xr.append: {(dt.datetime.now()-dtstart).total_seconds():.2f} sec')
        dtstart = dt.datetime.now()
        ds_face_concat = xr.concat(ds_face_list, dim=facedim)
        if not only_faces:
            ds_node_concat = xr.concat(ds_node_list, dim=nodedim)
            ds_edge_concat = xr.concat(ds_edge_list, dim=edgedim)
            #ds_rest_concat = xr.concat(ds_rest_list, dim=edgedim)
        print(f'>>time xr.concat: {(dt.datetime.now()-dtstart).total_seconds():.2f} sec')
        
        if only_faces:
            ds_merged = ds_face_concat
        else:
            ds_merged = xr.merge([ds_face_concat,ds_node_concat,ds_edge_concat])#,ds_rest_concat])
            
        varlist_onepart = list(partitions[0].variables.keys())
        varlist_merged = list(ds_merged.variables.keys())
        varlist_dropped_bool = ~pd.Series(varlist_onepart).isin(varlist_merged)
        varlist_dropped = pd.Series(varlist_onepart).loc[varlist_dropped_bool]
        if varlist_dropped_bool.any():
            print(f'WARNING: some variables dropped with merging of partitions:\n{varlist_dropped}')#\nOne partition with some of these:\n{ds_rest}')
        
        ds_merged = ds_merged.rename({facedim: merged_grid.face_dimension}) #TODO: why is this necessary?
        return xu.UgridDataset(ds_merged, grids=[merged_grid])

    mode = 'map_partitioned' #'net' 'map_single' 'map_partitioned'
    if mode=='net':
        file_nc = os.path.join(dir_testdata,r'DFM_3D_z_Grevelingen\computations\run01','Grevelingen_FM_grid_20190603_*_net.nc')
    elif mode=='map_partitioned':
        file_nc = os.path.join(dir_testdata,r'DFM_3D_z_Grevelingen\computations\run01','DFM_OUTPUT_Grevelingen-FM','Grevelingen-FM_*_map.nc')
        file_nc = [os.path.join(dir_testdata,r'DFM_3D_z_Grevelingen\computations\run01','DFM_OUTPUT_Grevelingen-FM',f'Grevelingen-FM_{i:04d}_map.nc') for i in range(3)] #works also with one file or list of some of the partion files
    elif mode=='map_single': #TODO: make it work for non-partitioned files also (skip many steps)
        file_nc = os.path.join(dir_testdata,r'DFM_sigma_curved_bend\DFM_OUTPUT_cb_3d\cb_3d_map.nc')
        
    chunks = {'time':1}
    merged = open_partitioned_dataset(file_nc,only_faces=False,chunks=chunks)
    
    fig,ax = plt.subplots()
    merged['mesh2d_flowelem_domain'].ugrid.plot()#edgecolor='face')
    
    fig,ax = plt.subplots()
    merged['mesh2d_flowelem_globalnr'].ugrid.plot()
    
    if mode!='net':
        fig,ax = plt.subplots()
        pc = merged['mesh2d_s1'].isel(time=3).ugrid.plot(vmin=-0.3,vmax=0.3,edgecolor='face') #TODO: geen sel op time geeft onduidelijke error
        fig,ax = plt.subplots()
        pc = merged['mesh2d_sa1'].isel(time=3,nmesh2d_layer=34).ugrid.plot(vmin=28.7,vmax=30,edgecolor='face') #TODO: have to sel nmesh2d_layer, otherwise error (with xarray you get a histogram, possible to give same behaviour?)
    
    


time_passed_all = (dt.datetime.now()-dtstart_all).total_seconds()
print(f'>>time passed: {time_passed_all:.2f} sec')

