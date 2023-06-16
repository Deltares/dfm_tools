# -*- coding: utf-8 -*-
"""
Created on Wed Oct 27 22:03:01 2021

@author: veenstra

"""

import os
import numpy as np
import matplotlib.pyplot as plt
plt.close('all')
import xarray as xr
import dfm_tools as dfmt


dir_testinput = r'c:\DATA\dfm_tools_testdata'
dir_output = '.'

file_nc_list = [r'p:\archivedprojects\1220688-lake-kivu\2_data\COSMO\COSMOCLM_2012_out02_merged_4Wouter.nc', #COSMO
                r'p:\11202255-sfincs\Testbed\Original_tests\01_Implementation\08_restartfile\sfincs_map.nc', #SFINCS
                ]

for file_nc in file_nc_list:
    basename = os.path.basename(file_nc).replace('.','')
    
    data_xr = xr.open_dataset(file_nc)
    if 'COSMOCLM' in file_nc: #add ldb for Kivu
        name_wl_x,name_wl_y,name_wl = None,None,None
        name_uv_x,name_uv_y,name_u,name_v = 'lon','lat','U_10M','V_10M'
        thinning = 2
        data_xr = data_xr.drop(['height_10m','height_2m']) #gives cleaner figure title
        clim_wl = [0,0.15]
        clim_uv = [0,5]
        add_ldb = True
        fig_args = dict(nrows=1,ncols=3,figsize=(14,6))
    elif 'sfincs' in file_nc:
        name_wl_x,name_wl_y,name_wl = 'x','y','zs'
        name_uv_x,name_uv_y,name_u,name_v = 'edge_x','edge_y','u','v' #TODO: uv is now on centers, so use new mapfile (and simplify/lessen xyvarnames)
        thinning = 5
        data_xr = data_xr.set_coords(['x','y','edge_x','edge_y']) #set coordinates for sfincs. TODO: request x/y/etc as coords in sfincs mapfile https://github.com/Deltares/SFINCS/issues/10
        clim_wl = [0,0.15]
        clim_uv = [0,0.6]
        add_ldb = False
        fig_args = dict(nrows=3,ncols=1,figsize=(10,8))
    
    
    if name_wl is not None:
        data_fromnc_zs = data_xr[name_wl]
        fig, axs = plt.subplots(**fig_args)
        for iT, timestep in enumerate([0,1,10]):
            ax=axs[iT]
            pc = data_fromnc_zs.isel(time=timestep).plot.pcolormesh(x=name_wl_x,y=name_wl_y,cmap='jet', ax=ax)
            pc.set_clim(clim_wl)
            ax.set_aspect('equal')
        fig.tight_layout()
        plt.savefig(os.path.join(dir_output,f'{basename}_waterlevel'))
        
    
    data_uv = data_xr[[name_u,name_v]]
    data_u = data_xr[name_u]
    data_v = data_xr[name_v]
    
    fig, axs = plt.subplots(**fig_args)
    for iT, timestep in enumerate([0,1,10]):
        ax=axs[iT]
        coordlist = list(data_u.coords)
        coordlist.remove('time')
        data_u_tsel = data_u.isel(time=timestep)
        data_v_tsel = data_v.isel(time=timestep)
        #data_uv_tsel = data_uv.isel(time=timestep) #TODO: also possible to quiver a uv dataset directly
        data_magn = np.sqrt(data_u_tsel**2 + data_v_tsel**2)
        pc = data_magn.plot.pcolormesh(x=name_uv_x,y=name_uv_y, cmap='jet',ax=ax)
        #pc = data_uv_tsel.plot.quiver(x='lon',y='lat',u='U_10M',v='V_10M', cmap='jet',ax=ax)
        pc.set_clim(clim_uv)
        if add_ldb:
            dfmt.plot_coastlines(ax=ax)
        data_u_thin = data_u_tsel.loc[::thinning,::thinning]
        data_v_thin = data_v_tsel.loc[::thinning,::thinning]
        ax.quiver(data_u_thin[name_uv_x], data_u_thin[name_uv_y], data_u_thin, data_v_thin, #raises warnings for some reason
                  color='w')#,scale=50,width=0.008), cmap='jet')
    fig.tight_layout()
    plt.savefig(os.path.join(dir_output,f'{basename}_magn_pcolorquiver'))


    #curved vector plot #TODO: workinprogress
    dist = 0.106
    reg_x_vec = np.arange(data_xr[name_uv_x].min(),data_xr[name_uv_x].max()+dist,dist)
    reg_y_vec = np.arange(data_xr[name_uv_y].min(),data_xr[name_uv_y].max()+dist,dist)
    reg_grid_X,reg_grid_Y = np.meshgrid(reg_x_vec,reg_y_vec)
    from scipy.interpolate import griddata
    
    fig, axs = plt.subplots(**fig_args,sharex=True,sharey=True)
    for iT, timestep in enumerate([0,1,10]):
        ax=axs[iT]
        data_u_tsel = data_u.isel(time=timestep)
        data_v_tsel = data_v.isel(time=timestep)
        lonvals_flat = data_u_tsel[name_uv_x].to_numpy().flatten()
        latvals_flat = data_u_tsel[name_uv_y].to_numpy().flatten()
        U = griddata((lonvals_flat,latvals_flat),data_u_tsel.to_numpy().flatten(),(reg_grid_X,reg_grid_Y),method='nearest') #TODO: this is probably easier with xarray
        V = griddata((lonvals_flat,latvals_flat),data_v_tsel.to_numpy().flatten(),(reg_grid_X,reg_grid_Y),method='nearest')
        speed = np.sqrt(U*U + V*V)
        #quiv_curved = dfmt.velovect(ax,reg_grid_X,reg_grid_Y,U,V, arrowstyle='fancy', density=5, grains=25, color=speed)#, cmap='jet') #TODO: requested adding curved-quiver to matplotlib https://github.com/matplotlib/matplotlib/issues/20038
        #quiv_curved.lines.set_clim(0,5)
        quiv_curved = ax.streamplot(reg_x_vec,reg_y_vec,U,V, arrowstyle='-|>', integration_direction='forward',broken_streamlines=False, color=speed, density=1)
        cbar = fig.colorbar(quiv_curved.lines, ax=ax)
        cbar.set_label('velocity magnitude (%s)'%(data_v.attrs['units']))
        ax.set_title(data_v_tsel.time.to_pandas())
        if add_ldb:
            dfmt.plot_coastlines(ax=ax)
    fig.tight_layout()
    plt.savefig(os.path.join(dir_output,f'{basename}_magn_curvedquiver'))


