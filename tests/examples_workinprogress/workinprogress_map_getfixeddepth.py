# -*- coding: utf-8 -*-
"""
Created on Fri Dec  9 11:51:29 2022

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

file_nc_list = [#os.path.join(dir_testinput,'DFM_sigma_curved_bend\\DFM_OUTPUT_cb_3d\\cb_3d_map.nc'), #sigmalayer
                #os.path.join(dir_testinput,'DFM_3D_z_Grevelingen','computations','run01','DFM_OUTPUT_Grevelingen-FM','Grevelingen-FM_0000_map.nc'), #zlayer
                r'p:\1204257-dcsmzuno\2006-2012\3D-DCSM-FM\A18b_ntsu1\DFM_OUTPUT_DCSM-FM_0_5nm\DCSM-FM_0_5nm_0000_map.nc', #fullgrid
                #r'p:\11206813-006-kpp2021_rmm-2d\C_Work\31_RMM_FMmodel\computations\model_setup\run_207\results\RMM_dflowfm_0000_map.nc', #2D model
                #r'p:\archivedprojects\11203379-005-mwra-updated-bem\03_model\02_final\A72_ntsu0_kzlb2\DFM_OUTPUT_MB_02\MB_02_0000_map.nc',
                ]
timestep = 3

for file_nc in file_nc_list:
    print('processing %s'%(os.path.basename(file_nc)))
    basename = os.path.basename(file_nc).replace('.','')
    
    
    data_frommap_merged = dfmt.open_partitioned_dataset(file_nc.replace('_0000_','_0*_')) #TODO: make starred default, but not supported by older code
    """
    #example for how subsetting works
    data_xr_grid = data_frommap_merged.ugrid.grid.to_dataset()
    face_nos = data_xr_grid.mesh2d_face_nodes.isel(mesh2d_nMax_face_nodes=0).load()
    print(face_nos)
    face_nnodecoords_x = data_xr_grid.mesh2d_node_x.isel(mesh2d_nNodes=face_nos)
    """
    #get ugrid data, vars informatin and grid units (latter from bedlevel coordinates)
    vars_pd = dfmt.get_ncvarproperties(file_nc=file_nc)
    
    print('plot grid and values from mapdata (salinity on layer, 3dim, on cell centers) >> on fixed depth (currently only for z-layer models)')
    data_frommap_sel = data_frommap_merged['mesh2d_sa1'].isel(time=timestep) #select data for all layers
    
    #TODO: retrieve on z-depth is not possible yet for fullgrid/zsigma (mesh2d_flowelem_zw) or sigma (mesh2d_layer_sigma)
    fixed_depth = -50
    clim_sal = [32,36] #should be taken care of automatically, but max salinity >200 for some reason

    
    fig, ax = plt.subplots()
    if 'mesh2d_flowelem_zcc' in data_frommap_sel.coords: #TODO: maybe do it a bit more exact with mesh2d_flowelem_zw (mesh2d_nInterfaces) instead
        bool_valid = data_frommap_merged.isel(time=timestep).mesh2d_flowelem_zw.min(dim='mesh2d_nInterfaces') <= fixed_depth #TODO: C:\Users\veenstra\Anaconda3\envs\dfm_tools_env\lib\site-packages\dask\array\reductions.py:640: RuntimeWarning: All-NaN slice encountered. return np.nanmax(x_chunk, axis=axis, keepdims=keepdims)
        bool_mindist = data_frommap_sel.nmesh2d_layer==abs(data_frommap_sel.mesh2d_flowelem_zcc - fixed_depth).argmin(dim='nmesh2d_layer').load()
        from dask.diagnostics import ProgressBar
        print('performing isel on fixed depth for zsigma/fullgrid model')
        with ProgressBar():
            #data_frommap_ondepth = data_frommap_sel.isel(nmesh2d_layer=isel_array) #TODO: this results in <xarray.DataArray 'mesh2d_sa1' (mesh2d_nFaces: 289713, nmesh2d_layer: 289713)> (should not have layer dimension) #Huite commented that .where(boolean) works here, not .isel(ints): https://stackoverflow.com/questions/66027917/python-xarray-vectorized-indexing/66126022#66126022
            data_frommap_ondepth = data_frommap_sel.where(bool_valid&bool_mindist).max(dim='nmesh2d_layer') #set all layers but one to nan, than reduce (take max value in this case, but is arbitrary)
        if len(data_frommap_ondepth.shape)>1:
            raise Exception('should have 1 dimension')
        
    elif 'mesh2d_layer_z' in data_frommap_sel.coords:
        data_frommap_sel = data_frommap_sel.set_index({'nmesh2d_layer':'mesh2d_layer_z'}) #set depth as index on layers
        data_frommap_ondepth = data_frommap_sel.interp(nmesh2d_layer=fixed_depth) #interpolate to fixed z-depth
    elif 'mesh2d_layer_sigma' in data_frommap_sel.coords:
        raise Exception('mesh2d_layer_sigma not yet supported')
    else:
        print('WARNING: depth variable not found, returning the same dataarray')
        data_frommap_ondepth = data_frommap_sel
    
    pc = data_frommap_ondepth.ugrid.plot(edgecolor='face',cmap='jet')
    pc.set_clim(clim_sal)
    ax.set_aspect('equal')
    fig.tight_layout()
    fig.savefig(os.path.join(dir_output,f'{basename}_mesh2d_sa1_onfixeddepth'))