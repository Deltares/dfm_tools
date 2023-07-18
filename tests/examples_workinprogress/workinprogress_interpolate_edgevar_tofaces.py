# -*- coding: utf-8 -*-
"""
Created on Thu Apr 13 14:04:56 2023

@author: veenstra
"""

import dfm_tools as dfmt
import matplotlib.pyplot as plt
plt.close('all')
import xarray as xr

# file_nc = dfmt.data.fm_curvedbend_map(return_filepath=True) #sigmalayer
file_nc = dfmt.data.fm_grevelingen_map(return_filepath=True) #zlayer
# file_nc = r'p:\dflowfm\maintenance\JIRA\05000-05999\05477\c103_ws_3d_fourier\DFM_OUTPUT_westerscheldt01_0subst\westerscheldt01_0subst_map.nc' #zsigma model without fullgrid output but with new ocean_sigma_z_coordinate variable
# file_nc = r'p:\archivedprojects\11206813-006-kpp2021_rmm-2d\C_Work\31_RMM_FMmodel\computations\model_setup\run_207\results\RMM_dflowfm_0*_map.nc' #2D model
# file_nc = r'p:\1204257-dcsmzuno\2006-2012\3D-DCSM-FM\A18b_ntsu1\DFM_OUTPUT_DCSM-FM_0_5nm\DCSM-FM_0_5nm_0*_map.nc' #fullgrid
# file_nc = r'p:\archivedprojects\11203379-005-mwra-updated-bem\03_model\02_final\A72_ntsu0_kzlb2\DFM_OUTPUT_MB_02\MB_02_0*_map.nc'

uds = dfmt.open_partitioned_dataset(file_nc.replace('0*','0000'))
if 'RMM_dflowfm' in file_nc:
    uds = uds.isel(time=[-1])

uds_edges = dfmt.Dataset_varswithdim(uds, uds.grid.edge_dimension)

#TODO: can also be done for all edge variables, re-add to original dataset as *onfaces variables?
if 'mesh2d_vicwwu' in uds.data_vars:
    varn_onedges = 'mesh2d_vicwwu'
elif 'mesh2d_czu' in uds.data_vars:
    varn_onedges = 'mesh2d_czu'
else:
    varn_onedges = 'mesh2d_edge_type' #if all else fails, interpolate this one

fig,ax = plt.subplots()
uds[varn_onedges].isel(time=-1,nmesh2d_interface=-2,mesh2d_nInterfaces=-2,missing_dims='ignore').ugrid.plot()

mesh2d_var = uds.grid.to_dataset().mesh2d

print('construct indexer')
dimn_faces = uds.grid.face_dimension
dimn_maxfn = mesh2d_var.attrs['max_face_nodes_dimension']
dimn_layer, dimn_interface = dfmt.get_vertical_dimensions(uds)
dimn_edges = uds.grid.edge_dimension
fill_value = uds.grid.fill_value

data_fec = xr.DataArray(uds.grid.face_edge_connectivity,dims=(dimn_faces,dimn_maxfn)) # (8355, 4)
data_fec_validbool = data_fec!=fill_value
data_fec = data_fec.where(data_fec_validbool,-1).astype(int)
data_fec = data_fec.compute()

#TODO: interpolation is slow for many timesteps, so maybe use .sel() on time dimension first
print('interpolation with indexer: step 1 (for each face, select all corresponding edge values)')
edgevar_tofaces_onint_step1 = uds[varn_onedges].isel({dimn_edges:data_fec})
print('interpolation with indexer: step 2 (replace nonexistent edges with nan)')
if hasattr(data_fec,'_FillValue'):
    edgevar_tofaces_onint_step2 = edgevar_tofaces_onint_step1.where(data_fec_validbool) #replace all values for fillvalue edges (-1) with nan
else:
    edgevar_tofaces_onint_step2 = edgevar_tofaces_onint_step1
print('interpolation with indexer: step 3 (average edge values per face)')
edgevar_tofaces_onint = edgevar_tofaces_onint_step2.mean(dim=dimn_maxfn)
print('interpolation with indexer: done')


if dimn_interface in edgevar_tofaces_onint.dims:
    print('average from interfaces to layers (so in z-direction) in case of a 3D model')
    #select all top interfaces and all bottom interfaces, sum, divide by two (same as average)
    edgevar_tofaces_topint = edgevar_tofaces_onint.isel({dimn_interface:slice(1,None)})
    edgevar_tofaces_botint = edgevar_tofaces_onint.isel({dimn_interface:slice(None,-1)})
    edgevar_tofaces = (edgevar_tofaces_topint + edgevar_tofaces_botint)/2
    #rename int to lay dimension and re-assign variable attributes
    edgevar_tofaces = edgevar_tofaces.rename({dimn_interface:dimn_layer}).assign_attrs(edgevar_tofaces_onint_step1.attrs)
else:
    edgevar_tofaces = edgevar_tofaces_onint.assign_attrs(edgevar_tofaces_onint_step1.attrs)

fig,ax = plt.subplots()
edgevar_tofaces.isel(time=-1,nmesh2d_layer=-1,mesh2d_nLayers=-1,missing_dims='ignore').ugrid.plot()
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