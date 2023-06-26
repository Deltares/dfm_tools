# -*- coding: utf-8 -*-
"""
Created on Wed Oct 27 14:35:12 2021

@author: veenstra
"""

import os
import sys
import matplotlib.pyplot as plt
plt.close('all')
import numpy as np
import xarray as xr
import dfm_tools as dfmt
import contextily as ctx

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

fig, ax = plt.subplots(figsize=(9,5))
pc = magn.plot(ax=ax,x='longitude',y='latitude')
ctx.add_basemap(ax=ax,source=ctx.providers.Esri.WorldStreetMap,crs='EPSG:4326', attribution=False)
dfmt.plot_coastlines(ax=ax,linewidth=1,res='l',min_area=1000)
fig.savefig(os.path.join(dir_output,'cartopy_hirlam_moreoptions'))

if 'cartopy' in sys.modules: #cartop is not a dfm_tools dependency, so only create this plot if it is installed
    import cartopy.crs as ccrs #install cartopy with `conda install cartopy -c conda-forge`
    fig, ax = plt.subplots(figsize=(6,7),subplot_kw={'projection': ccrs.EuroPP()}) #provide axis projection on initialisation, cannot be edited later on
    pc = magn.plot(ax=ax,x='longitude',y='latitude', transform=ccrs.PlateCarree(),add_colorbar=False)
    ax.coastlines(linewidth=1)
    ax.gridlines(draw_labels=True) #cannot use ax.get_xticks+ax.set_xticks since the data was transformed
    fig.savefig(os.path.join(dir_output,'cartopy_hirlam_curvedgridlines'))

#GREVELINGEN
file_nc_map = os.path.join(dir_testinput,'DFM_grevelingen_3D\\Grevelingen-FM_0*_map.nc')
data_frommap_merged = dfmt.open_partitioned_dataset(file_nc_map)
fig, ax = plt.subplots(1,1)
pc = data_frommap_merged['mesh2d_flowelem_bl'].ugrid.plot(ax=ax, linewidth=0.5, cmap='jet', vmin=-40, vmax=10)
ctx.add_basemap(ax=ax,source=ctx.providers.Esri.WorldImagery,crs='EPSG:28992', attribution=False)
dfmt.plot_coastlines(ax=ax,linewidth=1,res='h',crs='EPSG:28992',min_area=100)

fig.savefig(os.path.join(dir_output,'cartopy_grevelingen_RD'))
