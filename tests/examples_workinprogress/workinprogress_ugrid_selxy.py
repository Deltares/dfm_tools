# -*- coding: utf-8 -*-
"""
Created on Mon Dec 12 07:45:11 2022

@author: veenstra
"""


import os
import matplotlib.pyplot as plt
plt.close('all')
import dfm_tools as dfmt
import xugrid as xu

dir_testinput = r'c:\DATA\dfm_tools_testdata'
dir_output = '.'

file_nc_list = [#os.path.join(dir_testinput,'DFM_sigma_curved_bend\\DFM_OUTPUT_cb_3d\\cb_3d_map.nc'), #sigmalayer
                os.path.join(dir_testinput,'DFM_3D_z_Grevelingen','computations','run01','DFM_OUTPUT_Grevelingen-FM','Grevelingen-FM_0000_map.nc'), #zlayer
                #r'p:\1204257-dcsmzuno\2006-2012\3D-DCSM-FM\A18b_ntsu1\DFM_OUTPUT_DCSM-FM_0_5nm\DCSM-FM_0_5nm_0000_map.nc', #fullgrid
                #r'p:\11206813-006-kpp2021_rmm-2d\C_Work\31_RMM_FMmodel\computations\model_setup\run_207\results\RMM_dflowfm_0000_map.nc', #2D model
                #r'p:\archivedprojects\11203379-005-mwra-updated-bem\03_model\02_final\A72_ntsu0_kzlb2\DFM_OUTPUT_MB_02\MB_02_0000_map.nc',
                ]

for file_nc in file_nc_list:
    print('processing %s'%(os.path.basename(file_nc)))
    basename = os.path.basename(file_nc).replace('.','')
    
    if 'cb_3d_map' in file_nc:
        timestep = 72
        layno = 5
        xslice = slice(750,2500)
        val_ylim = None
        clim_bl = None
        clim_sal = None
        crs = None
        file_nc_fou = None
    elif 'Grevelingen' in file_nc:
        timestep = 3
        layno = 33 #35 is top
        xslice = slice(45000,72000)#slice(50000,62000)
        yslice = slice(416000,421000)
        val_ylim = [-25,5]
        clim_bl = None
        clim_sal = [28,30.2]
        crs = "EPSG:28992"
        file_nc_fou = None
    elif 'DCSM-FM_0_5nm' in file_nc:
        timestep = 365
        layno = 45
        xslice = slice(-3,1)
        val_ylim = [-600,1]
        clim_bl = [-500,0]
        clim_sal = [25,36]
        crs = "EPSG:4326"
        file_nc_fou = None
    elif 'RMM_dflowfm' in file_nc:
        timestep = 365 #50
        layno = None
        xslice = slice(None,None)
        val_ylim = None
        clim_bl = [-10,10]
        clim_sal = None
        crs = "EPSG:28992"
    else:
        raise Exception('ERROR: no settings provided for this mapfile')
    
    file_nc_list = [file_nc.replace('_0000_',f'_{x:04d}_') for x in [2,3]]#][0,1,4,5,6,7]]
    data_frommap_merged = dfmt.open_partitioned_dataset(file_nc_list)#file_nc.replace('_0000_','_0*_'))
    #data_frommap_merged = xu.open_dataset(file_nc.replace('_0000_','_0002_'),chunks={'time':1})
    
    data_frommap_merged = data_frommap_merged.ugrid.sel(y=yslice) #TODO xugrid: if sel() results in empty, plotting crashes with "ImportError: Plotting of arrays of cftime.datetime objects or arrays indexed by cftime.datetime objects requires the optional `nc-time-axis` (v1.2.0 or later) package.". If using xugrid directly "AttributeError: 'DataArray' object has no attribute 'ugrid'"
    
    print('plot grid and bedlevel (constantvalue, 1 dim)')
    #get bedlevel and create plot with ugrid and cross section line
    fig, ax_input = plt.subplots()
    pc = data_frommap_merged['mesh2d_flowelem_domain'].ugrid.plot(edgecolor='grey')#,cmap='jet')
    pc.set_clim(clim_bl)
    ax_input.set_aspect('equal')
    line, = ax_input.plot([], [],'o-') # empty line
    #linebuilder = dfmt.LineBuilder(line) #this makes it possible to interactively click a line in the bedlevel figure. Use linebuilder.line_array as alternative line_array
    fig.tight_layout()
    
    
if 1:#'Grevelingen' in file_nc:
    
    uds1 = xu.open_dataset(file_nc.replace('_0000_','_0000_'),chunks={'time':1})
    fig, ax = plt.subplots(1,1)
    pc = uds1['mesh2d_flowelem_domain'].ugrid.plot(ax=ax, edgecolor='r',alpha=0.6,linewidth=2)#,cmap='jet')
    ax.set_aspect('equal')
    
    uds2 = xu.open_dataset(file_nc.replace('_0000_','_0002_'),chunks={'time':1})
    fig, ax = plt.subplots(1,1)
    pc = uds2['mesh2d_flowelem_domain'].ugrid.plot(ax=ax, edgecolor='r',alpha=0.6,linewidth=2)#,cmap='jet')
    ax.set_aspect('equal')
    ax.set_xlim(50200,51300)
    ax.set_ylim(420000,420900)



