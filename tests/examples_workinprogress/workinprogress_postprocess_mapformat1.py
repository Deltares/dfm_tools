# -*- coding: utf-8 -*-
"""
Created on Tue Feb 28 12:59:52 2023

@author: veenstra
"""

import dfm_tools as dfmt
import xarray as xr
import xugrid as xu
import matplotlib.pyplot as plt
plt.close('all')


file_nc_list = [#r'p:\11206811-002-kpp-veerse-meer\model\runs_2011-2012\VM_WQ_3D_run9_c\DFM_OUTPUT_VM_WQ_3D\VM_WQ_3D_0000_20130101_000000_rst.nc', #mf1_rstfile (without topology var). TODO: fails since no grids are present (no variable with cf_role:mesh_topology attr, which can be reconstructed but also no node_coordinates present in dataset)
                r'p:\1230882-emodnet_hrsm\GTSMv5.0\runs\GM43_2000m_eu0900m_ITfac5p5_wx\gtsm_200s_2000m_eu0900m_ca2000m_v4_0000_net.nc', #mf1_netfile
                r'p:\1230882-emodnet_hrsm\GTSMv5.0\runs\GM43_2000m_eu0900m_ITfac5p5_wx\output\gtsm_model_0000_map.nc', #mf1_mapfile
                ]

for file_nc in file_nc_list:
        
    ds = xr.open_dataset(file_nc)
    #uds = xu.core.wrap.UgridDataset(ds) #does not work for mf1_mapfile since two facedims are present in the dataset, this is catched in dfmt.open_partitioned_dataset()
    uds = dfmt.open_partitioned_dataset(file_nc)
    print(uds.grid.face_dimension)
    print(uds.grid.node_dimension)
    print(uds.grid.edge_dimension)
    
    fig,ax = plt.subplots()
    if 's1' in uds.data_vars:
        uds.s1.isel(time=-1).ugrid.plot(vmin=-2,vmax=2,edgecolor='face')
    else: #netfile
        uds.grid.plot()
    #fig.savefig()