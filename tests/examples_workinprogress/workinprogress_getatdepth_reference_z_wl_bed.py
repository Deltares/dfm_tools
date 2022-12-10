# -*- coding: utf-8 -*-
"""
Created on Fri Dec  9 20:39:55 2022

@author: veenstra
"""


import os
import matplotlib.pyplot as plt
plt.close('all')
import numpy as np
import contextily as ctx
import datetime as dt
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
    if 'cb_3d_map' in file_nc:
        clim_sal = None
    elif 'Grevelingen' in file_nc:
        clim_sal = [28,30.2]
    elif 'DCSM-FM_0_5nm' in file_nc:
        clim_sal = [25,36]
        
    data_frommap_merged = dfmt.open_partitioned_dataset(file_nc.replace('_0000_','_0*_')) #TODO: make starred default, but not yet supported by older code
    #reconstruct_zw_zcc_fromz(data_frommap_merged.isel(mesh2d_nFaces=[10000,4000,4005],time=[0,1]).load())
    
    #get ugrid data, vars informatin and grid units (latter from bedlevel coordinates)
    vars_pd = dfmt.get_ncvarproperties(file_nc=file_nc)
    timestep = 3

    
    print('plot grid and values from mapdata (salinity on layer, 3dim, on cell centers) >> on fixed depth') #depth w.r.t. z0/waterlevel/bedlevel
    data_frommap_timesel = data_frommap_merged.isel(time=timestep) #select data for all layers
    data_frommap_timesel_ondepth = dfmt.get_mapdata_atdepth(data_xr_map=data_frommap_timesel, depth=-4, reference='z0')#, zlayer_z0_interp=True)
    fig, ax = plt.subplots()
    pc = data_frommap_timesel_ondepth['mesh2d_sa1'].ugrid.plot(edgecolor='face',cmap='jet')
    pc.set_clim(clim_sal)
    ax.set_aspect('equal')
    fig.tight_layout()


    