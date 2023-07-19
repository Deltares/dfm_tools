# -*- coding: utf-8 -*-
"""
Created on Thu Apr 13 14:04:56 2023

@author: veenstra
"""

import dfm_tools as dfmt
import matplotlib.pyplot as plt
plt.close('all')


file_nc = dfmt.data.fm_curvedbend_map(return_filepath=True) #sigmalayer
file_nc = dfmt.data.fm_grevelingen_map(return_filepath=True) #zlayer
# file_nc = r'p:\dflowfm\maintenance\JIRA\05000-05999\05477\c103_ws_3d_fourier\DFM_OUTPUT_westerscheldt01_0subst\westerscheldt01_0subst_map.nc' #zsigma model without fullgrid output but with new ocean_sigma_z_coordinate variable
# file_nc = r'p:\archivedprojects\11206813-006-kpp2021_rmm-2d\C_Work\31_RMM_FMmodel\computations\model_setup\run_207\results\RMM_dflowfm_0*_map.nc' #2D model
# file_nc = r'p:\1204257-dcsmzuno\2006-2012\3D-DCSM-FM\A18b_ntsu1\DFM_OUTPUT_DCSM-FM_0_5nm\DCSM-FM_0_5nm_0*_map.nc' #fullgrid
# file_nc = r'p:\archivedprojects\11203379-005-mwra-updated-bem\03_model\02_final\A72_ntsu0_kzlb2\DFM_OUTPUT_MB_02\MB_02_0*_map.nc'

uds = dfmt.open_partitioned_dataset(file_nc.replace('0*','0002')) 

if 'RMM_dflowfm' in file_nc: #interpolation is slow for many timesteps, so subsetting first (alternatively use different chunks argument)
    uds = uds.isel(time=slice(-1,None))
#uds_edges = dfmt.Dataset_varswithdim(uds, uds.grid.edge_dimension)

#TODO: can also be done for all edge variables, re-add to original dataset as *onfaces variables?
if 'mesh2d_vicwwu' in uds.data_vars:
    varn_onedges = 'mesh2d_vicwwu'
elif 'mesh2d_czu' in uds.data_vars:
    varn_onedges = 'mesh2d_czu'
else:
    varn_onedges = 'mesh2d_edge_type' #if all else fails, interpolate this one
uda_edge = uds[varn_onedges]

fig,ax = plt.subplots()
uda_edge.isel(time=-1,nmesh2d_interface=-2,mesh2d_nInterfaces=-2,missing_dims='ignore').ugrid.plot()

uda_face = dfmt.uda_edges_to_faces(uda_edge)
uds[f'{varn_onedges}_onfaces'] = uda_face

fig,ax = plt.subplots()
uda_face.isel(time=-1,nmesh2d_layer=-1,mesh2d_nLayers=-1,missing_dims='ignore').ugrid.plot()

