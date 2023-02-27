# -*- coding: utf-8 -*-
"""
Created on Tue Jun 15 15:32:15 2021

@author: buckman
"""

import os
import numpy as np
import dfm_tools as dfmt
import xarray as xr
import xugrid as xu

file_nc_ug = os.path.join(r'c:\DATA\dfm_tools_testdata','DFM_3D_z_Grevelingen','computations','run01','DFM_OUTPUT_Grevelingen-FM','Grevelingen-FM_0*_map.nc')
file_nc = r'p:\11206811-002-kpp-veerse-meer\model\runs_2011-2012\VM_WQ_3D_run9_c\DFM_OUTPUT_VM_WQ_3D\VM_WQ_3D_0000_20130101_000000_rst.nc'
dir_output = '.' #r'p:\11206811-002-kpp-veerse-meer\model\initial_conditions\VM_WQ_3D_run9c'

if not os.path.exists(dir_output):
    os.makedirs(dir_output)


mesh2d_dict_mf1 = {'cf_role': 'mesh_topology',
                 'long_name': 'Topology data of 2D network',
                 'topology_dimension': 2,
                 'node_coordinates': 'NetNode_x NetNode_y', #'mesh2d_node_x mesh2d_node_y', #TODO: NetNode_x/NetNode_y but missing in rst
                 'node_dimension': 'nNetNode', #'nmesh2d_node', #TODO: nNetNode missing in rst
                 'max_face_nodes_dimension': 'nNetElemMaxNode', # 'max_nmesh2d_face_nodes',
                 'edge_node_connectivity': 'NetLink', #'mesh2d_edge_nodes',
                 #'edge_dimension': 'nFlowLink', #'nmesh2d_edge',
                 'edge_dimension': 'nNetLink', #'nmesh2d_edge',
                 'edge_coordinates': 'FlowLink_xu FlowLink_yu', #'mesh2d_edge_x mesh2d_edge_y',
                 'face_node_connectivity': 'NetElemNode', #'mesh2d_face_nodes',
                 'face_dimension': 'nNetElem', #'nmesh2d_face', #TODO: rst (but also present in map)
                 #'face_dimension': 'nFlowElem', #'nmesh2d_face', TODO: or nNetElem
                 'edge_face_connectivity': 'NetElemLink', #'mesh2d_edge_faces',
                 #'edge_face_connectivity': 'FlowLink', #'mesh2d_edge_faces',
                 'face_coordinates': 'FlowElem_xzw FlowElem_yzw', #'mesh2d_face_x mesh2d_face_y', #TODO: rst, but also map
                 #'face_coordinates': 'FlowElem_xcc FlowElem_ycc', #'mesh2d_face_x mesh2d_face_y', #TODO: present in map, but is bounds?
                 'layer_dimension': 'laydim', #'nmesh2d_layer',
                 'interface_dimension': 'wdim', #'nmesh2d_interface',
                 'vertical_dimensions': 'laydim: wdim (padding: none)'} #'nmesh2d_layer: nmesh2d_interface (padding: none)'}
#nNetElem

ds_ug = xr.open_dataset(file_nc_ug.replace('0*','0000'))
print(ds_ug.mesh2d)

ds_rst = xr.open_dataset(file_nc)
ds_rst['mesh2d'] = xr.DataArray(data=int()).assign_attrs(mesh2d_dict_mf1)
#uds_rst = xu.core.wrap.UgridDataset(ds_rst)

ds_map = xr.open_dataset(r'p:\1230882-emodnet_hrsm\GTSMv5.0\runs\GM43_2000m_eu0900m_ITfac5p5_wx\output\gtsm_model_0000_map.nc')
ds_map['mesh2d'] = xr.DataArray(data=int()).assign_attrs(mesh2d_dict_mf1)
ds_map.NetLink
uds_map = xu.core.wrap.UgridDataset(ds_map)

breakit

uds = dfmt.open_partitioned_dataset(file_nc)

# define map variables for VM
subs = ['DetCS1','DetNS1','DetPS1','DetSiS1']
x_coords = dfmt.get_ncmodeldata(file_nc=file_nc, varname='FlowElem_xzw', multipart=True)
y_coords = dfmt.get_ncmodeldata(file_nc=file_nc, varname='FlowElem_yzw', multipart=True)

#data_frommap_merged = dfmt.open_partitioned_dataset(file_nc.replace('_0000_','_0*_')) #TODO: make starred default, but not supported by older code #TODO: 0000 not at end of filename and no mesh2d variable in the file

for sub in subs:
    #get data to plot
    data_frommap = dfmt.get_ncmodeldata(file_nc=file_nc, varname=sub, timestep=0, layer=None, multipart=True)
    xyz = np.stack([x_coords,y_coords,data_frommap.transpose().reshape(-1)])
    
    xyz_out = xyz.transpose()
    np.savetxt(os.path.join(sub+'.xyz'),xyz_out,delimiter='\t')
