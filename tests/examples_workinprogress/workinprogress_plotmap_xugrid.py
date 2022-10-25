# -*- coding: utf-8 -*-
"""
Created on Mon Sep 26 13:41:09 2022

@author: veenstra
"""
import matplotlib.pyplot as plt
plt.close('all')
import xugrid
import xarray as xr
import glob
import numpy as np
from dfm_tools.get_nc_helpers import get_ncvardimlist


# xugrid alternative for plotting a grid (use random face property and use facecolor='none')
file_net = r'c:\DATA\dfm_tools_testdata\DFM_3D_z_Grevelingen\computations\run01\Grevelingen_FM_grid_20190603_net.nc'
data_xru = xugrid.open_dataset(file_net)
fig,(ax1) = plt.subplots()
pc = data_xru.mesh2d_face_x.ugrid.plot(ax=ax1, facecolor='none', edgecolor='grey', linewidth=0.5, alpha=0.5, add_colorbar=False)


file_nc = r'c:\DATA\dfm_tools_testdata\DFM_3D_z_Grevelingen\computations\run01\DFM_OUTPUT_Grevelingen-FM\Grevelingen-FM_00*_map.nc'
file_list = glob.glob(file_nc)

clim_min,clim_max = np.nan,np.nan
list_pc = []

fig,(ax1) = plt.subplots()
for iF, file_nc in enumerate(file_list):#[:2]):
    data_xr = xugrid.open_dataset(file_nc)
    
    vars_pd,dims_pd = get_ncvardimlist(file_nc)
    data_dom = data_xr['mesh2d_flowelem_domain']
    bool_nonghost = data_dom==iF
    data_dom_nonghost = data_dom.sel(nmesh2d_face=bool_nonghost)
    data_wl = data_xr['mesh2d_sa1'].isel(time=3,nmesh2d_layer=-1)#.sel(nmesh2d_face=bool_nonghost)
    #data_wl = data_xr.get('mesh2d_s1').isel(time=0)
    
    pc = data_wl.ugrid.plot(ax=ax1,cmap='jet',add_colorbar=False)#,edgecolor='face')#,vmin=-1,vmax=1)
    list_pc.append(pc)
    clim_min = np.nanmin([clim_min,data_wl.min()])
    clim_max = np.nanmin([clim_max,data_wl.max()])
    
    #data_dom_nonghost.ugrid.plot(ax=ax1,cmap='jet',vmin=0,vmax=len(file_list))

for pc in list_pc:
    pc.set_clim(clim_min,clim_max)
    #pc.set_clim(27,29.1)
    
ax1.set_xlim(46500,71000)
ax1.set_ylim(408000,426000)
fig.colorbar(pc,ax=ax1)
ax1.set_aspect('equal')

