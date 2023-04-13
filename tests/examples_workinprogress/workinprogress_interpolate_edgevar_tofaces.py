# -*- coding: utf-8 -*-
"""
Created on Thu Apr 13 14:04:56 2023

@author: veenstra
"""

import os
import dfm_tools as dfmt

dir_testinput = r'c:\DATA\dfm_tools_testdata'

file_nc = os.path.join(dir_testinput,'DFM_sigma_curved_bend\\DFM_OUTPUT_cb_3d\\cb_3d_map.nc') #sigmalayer
file_nc = os.path.join(dir_testinput,'DFM_3D_z_Grevelingen','computations','run01','DFM_OUTPUT_Grevelingen-FM','Grevelingen-FM_0*_map.nc') #zlayer
# file_nc = r'p:\dflowfm\maintenance\JIRA\05000-05999\05477\c103_ws_3d_fourier\DFM_OUTPUT_westerscheldt01_0subst\westerscheldt01_0subst_map.nc' #zsigma model without fullgrid output but with new ocean_sigma_z_coordinate variable
# file_nc = r'p:\archivedprojects\11206813-006-kpp2021_rmm-2d\C_Work\31_RMM_FMmodel\computations\model_setup\run_207\results\RMM_dflowfm_0*_map.nc' #2D model
# file_nc = r'p:\1204257-dcsmzuno\2006-2012\3D-DCSM-FM\A18b_ntsu1\DFM_OUTPUT_DCSM-FM_0_5nm\DCSM-FM_0_5nm_0*_map.nc' #fullgrid
# file_nc = r'p:\archivedprojects\11203379-005-mwra-updated-bem\03_model\02_final\A72_ntsu0_kzlb2\DFM_OUTPUT_MB_02\MB_02_0*_map.nc'

uds = dfmt.open_partitioned_dataset(file_nc.replace('0*','0000'))

uds_edges = dfmt.Dataset_varswithdim(uds, uds.grid.edge_dimension)

if 'mesh2d_vicwwu' in uds.data_vars:
    varn_onedges = 'mesh2d_vicwwu'
elif 'mesh2d_czu' in uds.data_vars:
    varn_onedges = 'mesh2d_czu'
else:
    varn_onedges = 'mesh2d_edge_type' #if all else fails, interpolate this one

print('construct indexer')
varn_fnc = uds.grid.to_dataset().mesh2d.attrs['face_node_connectivity']
dimn_maxfn = uds.grid.to_dataset().mesh2d.attrs['max_face_nodes_dimension']
dimn_interface = uds.grid.to_dataset().mesh2d.attrs['interface_dimension']

data_fnc = uds.grid.to_dataset()[varn_fnc]

if hasattr(data_fnc,'_FillValue'):
    data_fnc_validbool = data_fnc!=data_fnc.attrs['_FillValue']
else:
    data_fnc_validbool = None

if hasattr(data_fnc,'start_index'):
    if data_fnc.attrs['start_index'] != 0:
        data_fnc = data_fnc - data_fnc.attrs['start_index'] #TODO: this drops attrs, avoid this

#TODO: interpolation is slow for many timesteps, so maybe use .sel() on time dimension first
print('interpolation with indexer: step 1')
edgevar_tofaces_onint_step1 = uds[varn_onedges].isel({uds.grid.edge_dimension:data_fnc}) #TODO: fails for cb_3d_map.nc, westernscheldt
print('interpolation with indexer: step 2')
edgevar_tofaces_onint_step2 = edgevar_tofaces_onint_step1.where(data_fnc_validbool) #replace all values for fillvalue edges (-1) with nan
print('interpolation with indexer: step 3')
edgevar_tofaces_onint = edgevar_tofaces_onint_step2.mean(dim=dimn_maxfn)
print('interpolation with indexer: done')

if dimn_interface in edgevar_tofaces_onint.dims:
    #TODO: interpolate interfaces to centers (first we need zcc variable also, but this is 3D so inconvenient)
    #TODO: easier might be to shift variable, add that to original variable, drop top/bottom values, change dimname from interfaces to centers
    #edgevar_tofaces = edgevar_tofaces_onint.interp()#{dimn_interface:})
    pass
else:
    edgevar_tofaces = edgevar_tofaces_onint
#TODO: add inverse distance weighing, below example is from faces to edges, so make other way round
"""
edge_coords = grid.edge_coordinates
edge_faces = grid.edge_face_connectivity
boundary = (edge_faces[:, 1] == -1)
edge_faces[boundary, 1] = edge_faces[boundary, 0]
face_coords = grid.face_coordinates[edge_faces]
distance = np.linalg.norm(edge_coords[:, np.newaxis, :] - face_coords, axis=-1)
weights = distance / distance.sum(axis=1)[:, np.newaxis]
values = (face_values[edge_faces] * weights).sum(axis=1)
"""