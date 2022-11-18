# -*- coding: utf-8 -*-
"""
Created on Mon Sep 26 13:41:09 2022

@author: veenstra
"""

import os
import matplotlib.pyplot as plt
plt.close('all')
import datetime as dt
import dfm_tools as dfmt

dtstart_all = dt.datetime.now()

#dir_model = r"c:\tmp\xugrid\Grevelingen_run01"
dir_testdata = r'c:\DATA\dfm_tools_testdata'


if 0: #REFERENCE DFM_TOOLS
    
    file_nc = os.path.join(dir_testdata,r'DFM_3D_z_Grevelingen\computations\run01','DFM_OUTPUT_Grevelingen-FM','Grevelingen-FM_0000_map.nc')
    crs = "EPSG:28992"
    
    timestep = 3
    layno = 33 #35 is top
    
    #get ugrid data
    ugrid_all = dfmt.get_netdata(file_nc=file_nc)
    
    print('plot grid and values from mapdata (salinity on layer, 3dim, on cell centers)')
    data_frommap = dfmt.get_ncmodeldata(file_nc=file_nc, varname='mesh2d_s1', timestep=timestep)#, layer=layno)
    fig, ax = plt.subplots()
    pc = dfmt.plot_netmapdata(ugrid_all.verts, values=data_frommap, ax=None, linewidth=0.5, cmap="jet", edgecolor='face')
    cbar = fig.colorbar(pc, ax=ax)
    cbar.set_label('%s [%s]'%(data_frommap.var_ncattrs['long_name'], data_frommap.var_ncattrs['units']))
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    


if 1: #xugrid

    mode = 'map_partitioned' #'net' 'map_single' 'map_partitioned'
    if mode=='net':
        file_nc = os.path.join(dir_testdata,r'DFM_3D_z_Grevelingen\computations\run01','Grevelingen_FM_grid_20190603_*_net.nc')
        crs = "EPSG:28992"
    elif mode=='map_partitioned':
        file_nc = os.path.join(dir_testdata,r'DFM_3D_z_Grevelingen\computations\run01','DFM_OUTPUT_Grevelingen-FM','Grevelingen-FM_*_map.nc')
        #file_nc = [os.path.join(dir_testdata,r'DFM_3D_z_Grevelingen\computations\run01','DFM_OUTPUT_Grevelingen-FM',f'Grevelingen-FM_{i:04d}_map.nc') for i in [1]]#,2,3]]#range(3)] #works also with one file or list of some of the partion files
        crs = "EPSG:28992"
        layer = 34
    elif mode=='map_single':
        file_nc = os.path.join(dir_testdata,r'DFM_sigma_curved_bend\DFM_OUTPUT_cb_3d\cb_3d_map.nc')
        crs = None
        layer = 5
    
    chunks = {'time':1}
    merged = dfmt.open_partitioned_dataset(file_nc,only_faces=False,chunks=chunks)
    
    #merged = merged.where(merged.mesh2d_face_y<420000) nan in case of false
    #merged = merged.ugrid.sel(y=slice(None,420000)) #TODO: this should work, HB will make an issue (new resease required)
    
    gridname = merged.ugrid.grid.name #'mesh2d'
    if mode!='map_single':
        fig,ax = plt.subplots()
        merged[f'{gridname}_flowelem_domain'].ugrid.plot()#edgecolor='face')
        fig,ax = plt.subplots()
        merged[f'{gridname}_flowelem_globalnr'].ugrid.plot()
    
    if mode!='net':
        fig,ax = plt.subplots()
        pc = merged[f'{gridname}_s1'].isel(time=3).ugrid.plot(vmin=-0.3,vmax=0.3,edgecolor='face')
        fig,ax = plt.subplots()
        pc = merged[f'{gridname}_sa1'].isel(time=3,nmesh2d_layer=layer).ugrid.plot(edgecolor='face')#,vmin=28.7,vmax=30) #TODO: have to sel time and nmesh2d_layer, otherwise error (with xarray you get a histogram, possible to give same behaviour?)
        # if crs is not None:
        #     source = None #ctx.providers.Esri.WorldImagery # ctx.providers.Stamen.Terrain (default), ctx.providers.CartoDB.Voyager, ctx.providers.NASAGIBS.ViirsEarthAtNight2012, ctx.providers.Stamen.Watercolor
        #     ctx.add_basemap(ax=ax, source=source, crs=crs, attribution=False)


time_passed_all = (dt.datetime.now()-dtstart_all).total_seconds()
print(f'>>time passed: {time_passed_all:.2f} sec')

