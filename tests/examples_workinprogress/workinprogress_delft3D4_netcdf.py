# -*- coding: utf-8 -*-
"""
Created on Wed Oct 27 21:54:57 2021

@author: veenstra

to get delft3D to write netCDF output instead of .dat files, add these lines to your mdf:
    FlNcdf= #maphis#
    ncFormat=4

"""

import os
import matplotlib.pyplot as plt
plt.close('all')
import dfm_tools as dfmt

dir_testinput = r'c:\DATA\dfm_tools_testdata'
dir_output = '.'

file_nc_list = [r'p:\archivedprojects\1220688-lake-kivu\3_modelling\1_FLOW\7_heatfluxinhis\062_netcdf\trim-thiery_002_coarse.nc',
                #os.path.join(dir_testinput,'D3D_3D_sigma_curved_bend_nc\\trim-cb2-sal-added-3d.nc'),
                ]

for file_nc in file_nc_list:
    
    if 'trim-cb2-sal' in file_nc:
        name = 'curvedbend'
        timestep = 4
        res = 100
        scale = 25
        figsize = (6,5)
        def maybe_add_coastlines(ax):
            ax.set_aspect('equal')
    else:
        name = 'kivu'
        timestep = 10
        res = 1/70
        scale = 3
        figsize = (6,7)
        def maybe_add_coastlines(ax):
            dfmt.plot_coastlines(ax=ax,res='h')
    
    uds = dfmt.open_dataset_delft3d4(file_nc)
    uds_sel = uds.isel(time=timestep,KMAXOUT_RESTR=-2)
    uds_raster = dfmt.rasterize_ugrid(uds_sel,resolution=res)
    
    fig,ax = plt.subplots(figsize=figsize)
    uds.grid.plot(ax=ax, color='b',linewidth=0.2)
    maybe_add_coastlines(ax)
    fig.savefig(os.path.join(dir_output,f'{name}_mesh'))
    
    fig,ax = plt.subplots(figsize=figsize)
    uds_sel.DPS0.ugrid.plot(ax=ax, center=False, cmap='jet')
    maybe_add_coastlines(ax)
    fig.savefig(os.path.join(dir_output,f'{name}_bedlevel'))
    
    fig,ax = plt.subplots(figsize=figsize)
    uds_sel.umag.ugrid.plot(ax=ax, center=False, cmap='jet')
    pc = ax.quiver(uds_raster.x, uds_raster.y, uds_raster.ux, uds_raster.uy,
              scale=scale,color='w',width=0.005, cmap='jet')
    maybe_add_coastlines(ax)
    fig.savefig(os.path.join(dir_output,f'{name}_umag'))
    
    fig,ax = plt.subplots(figsize=figsize)
    uds.grid.plot(ax=ax, color='b',linewidth=0.2, zorder=0, alpha=0.6)
    pc = ax.quiver(uds_raster.x, uds_raster.y, uds_raster.ux, uds_raster.uy, uds_raster.umag,
              scale=scale,color='w',width=0.005, cmap='jet')
    fig.colorbar(pc,ax=ax)
    maybe_add_coastlines(ax)
    fig.savefig(os.path.join(dir_output,f'{name}_umag_quiv'))

