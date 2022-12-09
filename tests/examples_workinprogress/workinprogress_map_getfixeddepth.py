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

file_nc_list = [os.path.join(dir_testinput,'DFM_sigma_curved_bend\\DFM_OUTPUT_cb_3d\\cb_3d_map.nc'), #sigmalayer
                #os.path.join(dir_testinput,'DFM_3D_z_Grevelingen','computations','run01','DFM_OUTPUT_Grevelingen-FM','Grevelingen-FM_0000_map.nc'), #zlayer
                #r'p:\1204257-dcsmzuno\2006-2012\3D-DCSM-FM\A18b_ntsu1\DFM_OUTPUT_DCSM-FM_0_5nm\DCSM-FM_0_5nm_0000_map.nc', #fullgrid
                #r'p:\11206813-006-kpp2021_rmm-2d\C_Work\31_RMM_FMmodel\computations\model_setup\run_207\results\RMM_dflowfm_0000_map.nc', #2D model
                #r'p:\archivedprojects\11203379-005-mwra-updated-bem\03_model\02_final\A72_ntsu0_kzlb2\DFM_OUTPUT_MB_02\MB_02_0000_map.nc',
                ]
timestep = 3
varname = 'mesh2d_sa1'
fixed_depth = -50
clim_sal = [32,36] #should be taken care of automatically, but max salinity >200 for some reason

for file_nc in file_nc_list:
    print('processing %s'%(os.path.basename(file_nc)))
    basename = os.path.basename(file_nc).replace('.','')
    
    
    data_frommap_merged = dfmt.open_partitioned_dataset(file_nc.replace('_0000_','_0*_')) #TODO: make starred default, but not supported by older code
    
    dtstart_all = dt.datetime.now()
    
    #get ugrid data, vars informatin and grid units (latter from bedlevel coordinates)
    vars_pd = dfmt.get_ncvarproperties(file_nc=file_nc)
    
    print('plot grid and values from mapdata (salinity on layer, 3dim, on cell centers) >> on fixed depth (currently only for z-layer models)')
    data_frommap_timesel = data_frommap_merged.isel(time=timestep) #select data for all layers
    
    def get_mapdata_onfixedepth(data_xr_map, z, varname=None):
        """
        data_xr_map:
            has to be dataset, otherwise mesh2d_flowelem_zw etc are not available (interface z values)
            in case of zsigma layers (or fullgrid), it cannot contain a time dimension
        #TODO: retrieve on z-depth is not possible yet for sigma (mesh2d_layer_sigma)
        """
        
        print('>> subsetting data on fixed depth: ',end='')
        dtstart = dt.datetime.now()
        
        if not 'nmesh2d_layer' in data_xr_map.dims: #TODO: maybe raise exception?
            print(f'{(dt.datetime.now()-dtstart).total_seconds():.2f} sec')
            print('WARNING: depth variable not found, returning input Dataset')
            return data_xr_map
        
        if varname is not None:
            data_xr_map_var = data_xr_map[varname]
        else:
            data_xr_map_var = data_xr_map
        #from dask.diagnostics import ProgressBar
        #with ProgressBar():
        if 'mesh2d_flowelem_zcc' in data_xr_map.coords:
            print('[zsigma/fullgrid model]',end='')
            if 'time' in data_xr_map.dims:
                raise Exception('ERROR: supply dataset with only one time for models with zsigma/sigma layers')
            bool_valid = data_xr_map.mesh2d_flowelem_zw.min(dim='mesh2d_nInterfaces') <= z #TODO suppress warning: C:\Users\veenstra\Anaconda3\envs\dfm_tools_env\lib\site-packages\dask\array\reductions.py:640: RuntimeWarning: All-NaN slice encountered. return np.nanmax(x_chunk, axis=axis, keepdims=keepdims)
            bool_mindist = data_xr_map.nmesh2d_layer==abs(data_xr_map.mesh2d_flowelem_zcc - z).argmin(dim='nmesh2d_layer').load()
            print('performing .where() on fixed depth for zsigma/fullgrid model')
            data_frommap_ondepth = data_xr_map_var.where(bool_valid&bool_mindist).max(dim='nmesh2d_layer') #set all layers but one to nan, followed by an arbitrary reduce (max in this case)
        elif 'mesh2d_layer_z' in data_xr_map.coords: #TODO: would be better to take interfaces into account also (it currently interpolates between z-center values)
            print('[z model]',end='')
            dict_layer_num2z = {'nmesh2d_layer':'mesh2d_layer_z'} #TODO: also transfer attrs from z to lay
            data_xr_map_var = data_xr_map_var.set_index(dict_layer_num2z).rename(dict_layer_num2z) #set depth as index on layers, to be able to interp to depths instead of layernumbers
            data_frommap_ondepth = data_xr_map_var.interp(mesh2d_layer_z=z,kwargs=dict(bounds_error=True)) #interpolate to fixed z-depth
        elif 'mesh2d_layer_sigma' in data_xr_map.coords:
            print('[sigma model]',end='')
            raise Exception('get_mapdata_onfixedepth does not support sigma-layer models yet (mesh2d_layer_sigma), let us know if you need it')
            if 'time' in data_xr_map.dims:
                raise Exception('ERROR: supply dataset with only one time for models with zsigma/sigma layers')
        else:
            raise Exception('layers present, but unknown layertype')
        print(f'{(dt.datetime.now()-dtstart).total_seconds():.2f} sec')
        
        return data_frommap_ondepth
    
    if 1: #subset variable after onfixeddepth function
        data_frommap_timesel_ondepth = get_mapdata_onfixedepth(data_xr_map=data_frommap_timesel, z=fixed_depth)
        data_frommap_timesel_ondepth = data_frommap_timesel_ondepth[varname]
    else:
        data_frommap_timesel_ondepth = get_mapdata_onfixedepth(data_xr_map=data_frommap_timesel, z=fixed_depth, varname=varname)
    
    fig, ax = plt.subplots()
    pc = data_frommap_timesel_ondepth.ugrid.plot(edgecolor='face',cmap='jet')
    pc.set_clim(clim_sal)
    ax.set_aspect('equal')
    fig.tight_layout()
    #fig.savefig(os.path.join(dir_output,f'{basename}_mesh2d_sa1_onfixeddepth'))
    print(f'script runtime (excl opening mapfile): {(dt.datetime.now()-dtstart_all).total_seconds():.2f} sec')
