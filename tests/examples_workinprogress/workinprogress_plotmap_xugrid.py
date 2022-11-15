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

dtstart = dt.datetime.now()

#dir_model = r"c:\tmp\xugrid\Grevelingen_run01"
dir_model = r'c:\DATA\dfm_tools_testdata\DFM_3D_z_Grevelingen\computations\run01'


if 0: #REFERENCE DFM_TOOLS
    
    file_nc = os.path.join(dir_model,'DFM_OUTPUT_Grevelingen-FM','Grevelingen-FM_0000_map.nc')
    
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
    


if 0: #only merging, goes wrong with numbering
    
    files_map = [os.path.join(dir_model,'DFM_OUTPUT_Grevelingen-FM',f"Grevelingen-FM_{i:04d}_map.nc") for i in range(8)]
    data_xr1 = xr.open_dataset(files_map[0])
    data_xr2 = xr.open_dataset(files_map[1])
    
    facedim, nodedim, edgedim = 'nmesh2d_face', 'nmesh2d_node', 'nmesh2d_edge'
    data_xr1_face = dfmt.Dataset_varswithdim(data_xr1,dimname=facedim)
    data_xr2_face = dfmt.Dataset_varswithdim(data_xr2,dimname=facedim)
    data_xr1_node = dfmt.Dataset_varswithdim(data_xr1,dimname=nodedim)
    data_xr2_node = dfmt.Dataset_varswithdim(data_xr2,dimname=nodedim)
    data_xr1_edge = dfmt.Dataset_varswithdim(data_xr1,dimname=edgedim)
    data_xr2_edge = dfmt.Dataset_varswithdim(data_xr2,dimname=edgedim)
    
    data_xr_merged_face = xr.concat([data_xr1_face,data_xr2_face],dim=facedim)
    data_xr_merged_node = xr.concat([data_xr1_node,data_xr2_node],dim=nodedim)
    data_xr_merged_edge = xr.concat([data_xr1_edge,data_xr2_edge],dim=edgedim)
    
    data_xr_merged = xr.merge([data_xr_merged_face,data_xr_merged_node,data_xr_merged_edge])
    data_xr_merged['mesh2d'] = data_xr1.mesh2d
    
    data_xu_merged = xu.UgridDataset(data_xr_merged)
    
    data_xu_merged['mesh2d_sa1'].isel(time=2,nmesh2d_layer=35).ugrid.plot()


if 1: #method huite, slow for multiple variables
    mode = 'map' #'net' 'map'
    
    if mode=='net':
        partitions = [xu.open_dataset(os.path.join(dir_model,f"Grevelingen_FM_grid_20190603_{i:04d}_net.nc")) for i in range(8)] #netfiles
        varname_idomain = 'idomain'
        varname_globalnr = 'iglobal_s'
    else:
        partitions = [xu.open_dataset(os.path.join(dir_model,'DFM_OUTPUT_Grevelingen-FM',f"Grevelingen-FM_{i:04d}_map.nc")) for i in range(8)] #mapfiles
        varname_idomain = 'mesh2d_flowelem_domain'
        varname_globalnr = 'mesh2d_flowelem_globalnr'
    
    idomains = [uds[varname_idomain] for uds in partitions]
    globalids = [uds[varname_globalnr] for uds in partitions]
    
    # %%
    def merge(partitions, face_variables=None):
        
        all_indices = []
        all_faces = []
        all_nodes_x = []
        all_nodes_y = []
        accumulator = 0
        for i, part in enumerate(partitions):
            # For ghost nodes, keep the values of the lower numbered partition.
            grid = part.ugrid.grid
            idx = np.flatnonzero(part[varname_idomain] >= i)
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
        start = 0
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
        selection = []
        for idx, uds in zip(all_indices, partitions):
            if face_variables is None: #make list of all face_variables (will be concatenated)
                face_variables = []
                for varname in uds.variables.keys():
                    if 'max_nmesh2d_face_nodes' in uds[varname].dims: #TODO: not possible to concatenate this dim yet (varies per partition)
                        continue
                    if facedim in uds[varname].dims:
                        face_variables.append(varname)
            ds = uds.ugrid.obj[face_variables]
            selection.append(ds.isel({facedim: idx}))
        merged_data = xr.concat(selection, dim=facedim)
        merged_data = merged_data.rename({facedim: merged_grid.face_dimension})
        return xu.UgridDataset(merged_data, grids=[merged_grid])
    # %%
    
    face_variables = None#[varname_idomain, varname_globalnr, 'mesh2d_s1']
    
    merged = merge(partitions, face_variables=face_variables)
    
    # fig,ax = plt.subplots()
    # merged[varname_idomain].ugrid.plot()
    
    # fig,ax = plt.subplots()
    # merged[varname_globalnr].ugrid.plot()
    
    fig,ax = plt.subplots()
    pc = merged['mesh2d_s1'].isel(time=3).ugrid.plot(vmin=-0.3,vmax=0.3,edgecolor='face') #TODO: geen sel op time geeft onduidelijke error
    
    


time_passed = (dt.datetime.now()-dtstart).total_seconds()
print(f'>>time passed: {time_passed:.2f} sec')


""" DELETE: mainly looping over files and plotting separately, no merging
# xugrid alternative for plotting a grid (use random face property and use facecolor='none')
file_net = r'c:\DATA\dfm_tools_testdata\DFM_3D_z_Grevelingen\computations\run01\Grevelingen_FM_grid_20190603_net.nc'
data_xru = xugrid.open_dataset(file_net)
fig,(ax1) = plt.subplots()
pc = data_xru.mesh2d_face_x.ugrid.plot(ax=ax1, facecolor='none', edgecolor='grey', linewidth=0.5, alpha=0.5, add_colorbar=False)


file_nc = r'c:\DATA\dfm_tools_testdata\DFM_3D_z_Grevelingen\computations\run01\DFM_OUTPUT_Grevelingen-FM\Grevelingen-FM_00*_map.nc'
file_list = glob.glob(file_nc)

clim_min,clim_max = np.nan,np.nan
list_pc = []

fig,(ax1) = plt.subplots()
for iF, file_nc in enumerate(file_list):#[:2]):
    data_xr = xugrid.open_dataset(file_nc)
    
    data_dom = data_xr['mesh2d_flowelem_domain']
    bool_nonghost = data_dom==iF
    data_dom_nonghost = data_dom.sel(nmesh2d_face=bool_nonghost)
    data_wl = data_xr['mesh2d_sa1'].isel(time=3,nmesh2d_layer=-1)#.sel(nmesh2d_face=bool_nonghost)
    #data_wl = data_xr.get('mesh2d_s1').isel(time=0)
    
    pc = data_wl.ugrid.plot(ax=ax1,cmap='jet',add_colorbar=False)#,edgecolor='face')#,vmin=-1,vmax=1)
    list_pc.append(pc)
    clim_min = np.nanmin([clim_min,data_wl.min()])
    clim_max = np.nanmin([clim_max,data_wl.max()])
    
    #data_dom_nonghost.ugrid.plot(ax=ax1,cmap='jet',vmin=0,vmax=len(file_list))

for pc in list_pc:
    pc.set_clim(clim_min,clim_max)
    #pc.set_clim(27,29.1)
    
ax1.set_xlim(46500,71000)
ax1.set_ylim(408000,426000)
fig.colorbar(pc,ax=ax1)
ax1.set_aspect('equal')
"""
