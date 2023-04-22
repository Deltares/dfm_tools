# -*- coding: utf-8 -*-
"""
Created on Wed Oct 27 14:35:12 2021

@author: veenstra
"""

import os
import matplotlib.pyplot as plt
plt.close('all')
import numpy as np
import xarray as xr
import dfm_tools as dfmt
import contextily as ctx
import cartopy.crs as ccrs #install cartopy with `conda install cartopy -c conda-forge`
import cartopy.feature as cf #install cartopy with `conda install cartopy -c conda-forge`

dir_testinput = r'c:\DATA\dfm_tools_testdata'
dir_output = '.'

#HIRLAM
file_nc = r'p:\1204257-dcsmzuno\2014\data\meteo\HIRLAM72_2018\h72_201803.nc' #TODO: xarray MissingDimensionsError
data_xr = xr.open_mfdataset(file_nc,drop_variables=['x','y'],preprocess=dfmt.preprocess_hirlam)
timestep = 0
coarsefac = 2 #coarsen dataset for more performance, but not necessary

data_u = data_xr['eastward_wind'].isel(time=timestep)
data_v = data_xr['northward_wind'].isel(time=timestep)
magn_attrs = {'standard_name':'velocity magnitude','units':data_u.attrs['units']}
magn = np.sqrt(data_u**2 + data_v**2).assign_attrs(magn_attrs)
magn = magn[::coarsefac,::coarsefac] #coarsening makes coordinate conversion faster

fig, ax = plt.subplots()
ax.pcolormesh(magn)
fig.savefig(os.path.join(dir_output,'cartopy_hirlam_raw'))

fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()})
pc = magn.plot(ax=ax,x='longitude',y='latitude',add_colorbar=False)
fig.savefig(os.path.join(dir_output,'cartopy_hirlam_aspect'))

fig, ax = plt.subplots(figsize=(9,5),subplot_kw={'projection': ccrs.PlateCarree()}) #provide axis projection on initialisation, cannot be edited later on
pc = magn.plot(ax=ax,x='longitude',y='latitude')
ax.coastlines(linewidth=1)
ax.add_feature(cf.BORDERS, linewidth=1, edgecolor='gray', facecolor='none')
ctx.add_basemap(ax=ax,source=ctx.providers.Esri.WorldStreetMap,crs='EPSG:4326', attribution=False)
ax.set_xticks(ax.get_xticks()[1:-1])
ax.set_yticks(ax.get_yticks()[1:-1])
fig.savefig(os.path.join(dir_output,'cartopy_hirlam_moreoptions'))

fig, ax = plt.subplots(figsize=(6,7),subplot_kw={'projection': ccrs.EuroPP()}) #provide axis projection on initialisation, cannot be edited later on
pc = magn.plot(ax=ax,x='longitude',y='latitude', transform=ccrs.PlateCarree(),add_colorbar=False)
ax.coastlines(linewidth=1)
ax.gridlines(draw_labels=True) #cannot use ax.get_xticks+ax.set_xticks since the data was transformed
fig.savefig(os.path.join(dir_output,'cartopy_hirlam_curvedgridlines'))

#GREVELINGEN
file_nc_map = os.path.join(dir_testinput,'DFM_3D_z_Grevelingen\\computations\\run01\\DFM_OUTPUT_Grevelingen-FM\\Grevelingen-FM_0*_map.nc')
data_frommap_merged = dfmt.open_partitioned_dataset(file_nc_map) #TODO: make starred default, but not supported by older code
fig, ax = plt.subplots(1,1, subplot_kw={'projection': ccrs.epsg(28992)}) #provide axis projection on initialisation, cannot be edited later on
pc = data_frommap_merged['mesh2d_flowelem_bl'].ugrid.plot(ax=ax, linewidth=0.5, cmap='jet', vmin=-40, vmax=10)
ctx.add_basemap(ax=ax,source=ctx.providers.Esri.WorldImagery,crs='EPSG:28992', attribution=False)
ax.coastlines(linewidth=1)
ax.set_xticks(ax.get_xticks()[1:-1])
ax.set_yticks(ax.get_yticks()[1:-1])
fig.savefig(os.path.join(dir_output,'cartopy_grevelingen_RD'))
