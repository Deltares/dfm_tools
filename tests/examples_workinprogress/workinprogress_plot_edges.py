# -*- coding: utf-8 -*-
"""
Created on Fri Dec  9 11:07:16 2022

@author: veenstra
"""

import os
import matplotlib.pyplot as plt
plt.close('all')
import dfm_tools as dfmt

dir_testinput = r'c:\DATA\dfm_tools_testdata'
dir_output = '.'


file_nc_list = [os.path.join(dir_testinput,'DFM_sigma_curved_bend\\DFM_OUTPUT_cb_3d\\cb_3d_map.nc'), #sigmalayer
                os.path.join(dir_testinput,'DFM_3D_z_Grevelingen','computations','run01','DFM_OUTPUT_Grevelingen-FM','Grevelingen-FM_0000_map.nc'), #zlayer
                r'p:\1204257-dcsmzuno\2006-2012\3D-DCSM-FM\A18b_ntsu1\DFM_OUTPUT_DCSM-FM_0_5nm\DCSM-FM_0_5nm_0000_map.nc', #fullgrid
                r'p:\11206813-006-kpp2021_rmm-2d\C_Work\31_RMM_FMmodel\computations\model_setup\run_207\results\RMM_dflowfm_0000_map.nc', #2D model
                #r'p:\archivedprojects\11203379-005-mwra-updated-bem\03_model\02_final\A72_ntsu0_kzlb2\DFM_OUTPUT_MB_02\MB_02_0000_map.nc',
                ]

for file_nc in file_nc_list:
    print('processing %s'%(os.path.basename(file_nc)))
    basename = os.path.basename(file_nc).replace('.','')
    
    data_frommap_merged = dfmt.open_partitioned_dataset(file_nc.replace('_0000_','_0*_'))
    
    fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,figsize=(12,8))
    pc = data_frommap_merged.mesh2d_edge_x.ugrid.plot(ax=ax1,cmap='jet')
    pc = data_frommap_merged.mesh2d_edge_y.ugrid.plot(ax=ax2,cmap='jet')
    if 'mesh2d_flowelem_domain' in data_frommap_merged.data_vars:
        pc = data_frommap_merged.mesh2d_flowelem_domain.ugrid.plot(ax=ax3,cmap='jet',edgecolor='face')
    if 'mesh2d_flowelem_globalnr' in data_frommap_merged.data_vars:
        pc = data_frommap_merged.mesh2d_flowelem_globalnr.ugrid.plot(ax=ax4,cmap='jet',edgecolor='face')
    fig.tight_layout()
    fig.savefig(os.path.join(dir_output,f'{basename}_plot_edges'))
