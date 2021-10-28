# -*- coding: utf-8 -*-
"""
Created on Wed Oct 27 14:42:52 2021

@author: veenstra
"""

import os
import matplotlib.pyplot as plt
plt.close('all')

from dfm_tools.get_nc import get_netdata, get_ncmodeldata, plot_netmapdata
from dfm_tools.get_nc_helpers import get_ncvardimlist
dir_testinput = r'c:\DATA\dfm_tools_testdata'
dir_output = '.'

file_nc_list = [os.path.join(dir_testinput,r'DFM_sigma_curved_bend\DFM_OUTPUT_cb_3d\cb_3d_map.nc'),
                #r'p:\11205258-006-kpp2020_rmm-g6\C_Work\08_RMM_FMmodel\computations\run_180\DFM_OUTPUT_RMM_dflowfm\RMM_dflowfm_0000_map.nc',
                os.path.join(dir_testinput,r'DFM_3D_z_Grevelingen\computations\run01\DFM_OUTPUT_Grevelingen-FM\Grevelingen-FM_0000_map.nc'),
                ]

for file_nc in file_nc_list:
    if 'cb_3d_map' in file_nc:
        timestep = 3
        layer = 5
        clim_bl = None
        clim_wl = [-0.5,1]
        clim_sal = None
        clim_tem = None
    elif 'Grevelingen-FM_0000_map' in file_nc:
        timestep = 3
        layer = 33
        clim_bl = None
        clim_wl = [-0.5,1]
        clim_sal = [28,30.2]
        clim_tem = [4,10]
    elif 'RMM_dflowfm_0000_map' in file_nc:
        timestep = 50
        layer = None
        clim_bl = [-10,10]
        clim_wl = [-2,2]
        clim_sal = None
        clim_tem = None
    else:
        raise Exception('ERROR: no settings provided for this mapfile')
        
    
    vars_pd, dims_pd = get_ncvardimlist(file_nc=file_nc)

    #PLOT GRID
    print('plot only grid from mapdata')
    ugrid_all = get_netdata(file_nc=file_nc)#,multipart=False)
    fig, ax = plt.subplots()
    pc = plot_netmapdata(ugrid_all.verts, values=None, ax=None, linewidth=0.5, color="crimson", facecolor="None")
    ax.set_xlabel('x-direction')
    ax.set_ylabel('y-direction')
    ax.set_aspect('equal')
    fig.tight_layout()
    fig.savefig(os.path.join(dir_output,'%s_grid'%(os.path.basename(file_nc).replace('.',''))))


    #PLOT bedlevel
    if not 'cb_3d_map' in file_nc:
        print('plot grid and values from mapdata (constantvalue, 1 dim)')
        ugrid = get_netdata(file_nc=file_nc)#,multipart=False)
        #iT = 3 #for iT in range(10):
        data_frommap = get_ncmodeldata(file_nc=file_nc, varname='mesh2d_flowelem_bl')#, multipart=False)
        data_frommap_flat = data_frommap.flatten()
        fig, ax = plt.subplots()
        pc = plot_netmapdata(ugrid.verts, values=data_frommap_flat, ax=None, linewidth=0.5, cmap="jet")
        pc.set_clim(clim_bl)
        cbar = fig.colorbar(pc, ax=ax)
        cbar.set_label('%s [%s]'%(data_frommap.var_ncattrs['long_name'], data_frommap.var_ncattrs['units']))
        coordnames_xy = data_frommap.var_ncattrs['coordinates'].split()
        varcoords_x = get_ncmodeldata(file_nc=file_nc, varname=coordnames_xy[0])
        varcoords_y = get_ncmodeldata(file_nc=file_nc, varname=coordnames_xy[1])
        ax.set_xlabel('%s [%s]'%(varcoords_x.var_ncattrs['long_name'],varcoords_x.var_ncattrs['units']))
        ax.set_ylabel('%s [%s]'%(varcoords_y.var_ncattrs['long_name'],varcoords_y.var_ncattrs['units']))
        ax.set_aspect('equal')
        fig.tight_layout()
        fig.savefig(os.path.join(dir_output,'%s_mesh2d_flowelem_bl'%(os.path.basename(file_nc).replace('.',''))))
        
    
    #PLOT water level on map
    print('plot grid and values from mapdata (waterlevel, 2dim)')
    data_frommap = get_ncmodeldata(file_nc=file_nc, varname='mesh2d_s1', timestep=timestep)#, multipart=False)
    data_frommap_flat = data_frommap.flatten()
    fig, ax = plt.subplots()
    pc = plot_netmapdata(ugrid_all.verts, values=data_frommap_flat, ax=None, linewidth=0.5, cmap="jet")
    pc.set_clim(clim_wl)
    cbar = fig.colorbar(pc, ax=ax)
    cbar.set_label('%s [%s]'%(data_frommap.var_ncattrs['long_name'], data_frommap.var_ncattrs['units']))
    coordnames_xy = data_frommap.var_ncattrs['coordinates'].split()
    varcoords_x = get_ncmodeldata(file_nc=file_nc, varname=coordnames_xy[0])
    varcoords_y = get_ncmodeldata(file_nc=file_nc, varname=coordnames_xy[1])
    ax.set_xlabel('%s [%s]'%(varcoords_x.var_ncattrs['long_name'],varcoords_x.var_ncattrs['units']))
    ax.set_ylabel('%s [%s]'%(varcoords_y.var_ncattrs['long_name'],varcoords_y.var_ncattrs['units']))
    ax.set_aspect('equal')
    fig.tight_layout()
    fig.savefig(os.path.join(dir_output,'%s_mesh2d_s1'%(os.path.basename(file_nc).replace('.',''))))

    #PLOT var layer on map
    if 'RMM_dflowfm_0000_map' in file_nc:
        print('plot grid and values from mapdata (wind x velocity on cell edges)')
        data_frommap = get_ncmodeldata(file_nc=file_nc, varname='mesh2d_windxu', timestep=timestep, layer=layer)#, multipart=False)
        data_frommap_flat = data_frommap.flatten()
        fig, ax = plt.subplots()
        pc = plot_netmapdata(ugrid_all.edge_verts, values=data_frommap_flat, ax=None, linewidth=0.5, cmap="jet", edgecolor='face')
        #pc.set_clim(0,5)
        cbar = fig.colorbar(pc, ax=ax)
        cbar.set_label('%s [%s]'%(data_frommap.var_ncattrs['long_name'], data_frommap.var_ncattrs['units']))
        coordnames_xy = data_frommap.var_ncattrs['coordinates'].split()
        varcoords_x = get_ncmodeldata(file_nc=file_nc, varname=coordnames_xy[0])
        varcoords_y = get_ncmodeldata(file_nc=file_nc, varname=coordnames_xy[1])
        ax.set_xlabel('%s [%s]'%(varcoords_x.var_ncattrs['long_name'],varcoords_x.var_ncattrs['units']))
        ax.set_ylabel('%s [%s]'%(varcoords_y.var_ncattrs['long_name'],varcoords_y.var_ncattrs['units']))
        ax.set_aspect('equal')
        fig.tight_layout()
        fig.savefig(os.path.join(dir_output,'%s_mesh2d_windxu_edges'%(os.path.basename(file_nc).replace('.',''))))
    else:
        print('plot grid and values from mapdata (salinity on layer, 3dim)')
        data_frommap = get_ncmodeldata(file_nc=file_nc, varname='mesh2d_sa1', timestep=timestep, layer=layer)#, multipart=False)
        data_frommap_flat = data_frommap.flatten()
        fig, ax = plt.subplots()
        pc = plot_netmapdata(ugrid_all.verts, values=data_frommap_flat, ax=None, linewidth=0.5, cmap="jet")
        pc.set_clim(clim_sal)
        cbar = fig.colorbar(pc, ax=ax)
        cbar.set_label('%s [%s]'%(data_frommap.var_ncattrs['long_name'], data_frommap.var_ncattrs['units']))
        coordnames_xy = data_frommap.var_ncattrs['coordinates'].split()
        varcoords_x = get_ncmodeldata(file_nc=file_nc, varname=coordnames_xy[0])
        varcoords_y = get_ncmodeldata(file_nc=file_nc, varname=coordnames_xy[1])
        ax.set_xlabel('%s [%s]'%(varcoords_x.var_ncattrs['long_name'],varcoords_x.var_ncattrs['units']))
        ax.set_ylabel('%s [%s]'%(varcoords_y.var_ncattrs['long_name'],varcoords_y.var_ncattrs['units']))
        ax.set_aspect('equal')
        fig.tight_layout()
        fig.savefig(os.path.join(dir_output,'%s_mesh2d_sa1'%(os.path.basename(file_nc).replace('.',''))))
    
        print('plot grid and values from mapdata (temperature on layer, 3dim)')
        data_frommap = get_ncmodeldata(file_nc=file_nc, varname='mesh2d_tem1', timestep=timestep, layer=layer)#, multipart=False)
        data_frommap_flat = data_frommap.flatten()
        fig, ax = plt.subplots()
        pc = plot_netmapdata(ugrid_all.verts, values=data_frommap_flat, ax=None, linewidth=0.5, cmap="jet")
        pc.set_clim(clim_tem)
        cbar = fig.colorbar(pc, ax=ax)
        cbar.set_label('%s [%s]'%(data_frommap.var_ncattrs['long_name'], data_frommap.var_ncattrs['units']))
        coordnames_xy = data_frommap.var_ncattrs['coordinates'].split()
        varcoords_x = get_ncmodeldata(file_nc=file_nc, varname=coordnames_xy[0])
        varcoords_y = get_ncmodeldata(file_nc=file_nc, varname=coordnames_xy[1])
        ax.set_xlabel('%s [%s]'%(varcoords_x.var_ncattrs['long_name'],varcoords_x.var_ncattrs['units']))
        ax.set_ylabel('%s [%s]'%(varcoords_y.var_ncattrs['long_name'],varcoords_y.var_ncattrs['units']))
        ax.set_aspect('equal')
        fig.tight_layout()
        fig.savefig(os.path.join(dir_output,'%s_mesh2d_tem1'%(os.path.basename(file_nc).replace('.',''))))
  
    