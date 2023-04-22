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

def add_ticks(ax, nbins='auto'):
    """
    cartopy.mpl.geoaxes.GeoAxesSubplot does not have xticks/yticks by default, this function uses matplotlibs MaxNLocation to automatically create xticks and yticks
    """
    import matplotlib as mpl
    #check if xticklabels are different than xticks
    xticks = ax.get_xticks()
    xticklabels = np.array([x.get_text().replace('−','-') for x in ax.get_xticklabels()]).astype(float)
    if not (xticks==xticklabels).all():
        raise Exception('you are transforming coordinates, please use ax.gridlines(draw_labels=True) to get proper axis ticks+ticklabels')
    else:
        extent = ax.get_extent()
        mpl_al = mpl.ticker.MaxNLocator(nbins=nbins, prune='both') #https://github.com/matplotlib/matplotlib/blob/v3.7.1/lib/matplotlib/ticker.py#L1957-L2166
        xticks = mpl_al.tick_values(vmin=extent[0],vmax=extent[1])
        yticks = mpl_al.tick_values(vmin=extent[2],vmax=extent[3])
        ax.set_xticks(xticks)
        ax.set_yticks(yticks)

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
add_ticks(ax)
fig.savefig(os.path.join(dir_output,'cartopy_hirlam_moreoptions'))

fig, ax = plt.subplots(figsize=(6,7),subplot_kw={'projection': ccrs.EuroPP()}) #provide axis projection on initialisation, cannot be edited later on
pc = magn.plot(ax=ax,x='longitude',y='latitude', transform=ccrs.PlateCarree(),add_colorbar=False)
ax.coastlines(linewidth=1)
ax.gridlines(draw_labels=True) #cannot use add_ticks() since we transformed the data
fig.savefig(os.path.join(dir_output,'cartopy_hirlam_curvedgridlines'))

#GREVELINGEN
file_nc_map = os.path.join(dir_testinput,'DFM_3D_z_Grevelingen\\computations\\run01\\DFM_OUTPUT_Grevelingen-FM\\Grevelingen-FM_0*_map.nc')
data_frommap_merged = dfmt.open_partitioned_dataset(file_nc_map) #TODO: make starred default, but not supported by older code
fig, ax = plt.subplots(1,1, subplot_kw={'projection': ccrs.epsg(28992)}) #provide axis projection on initialisation, cannot be edited later on
pc = data_frommap_merged['mesh2d_flowelem_bl'].ugrid.plot(ax=ax, linewidth=0.5, cmap='jet', vmin=-40, vmax=10)
ctx.add_basemap(ax=ax,source=ctx.providers.Esri.WorldImagery,crs='EPSG:28992', attribution=False)
ax.coastlines(linewidth=1)
add_ticks(ax)
fig.savefig(os.path.join(dir_output,'cartopy_grevelingen_RD'))
