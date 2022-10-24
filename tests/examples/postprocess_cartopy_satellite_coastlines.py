# -*- coding: utf-8 -*-
"""
Created on Wed Oct 27 14:35:12 2021

@author: veenstra
"""
import os

import matplotlib.pyplot as plt
plt.close('all')
import numpy as np
import cartopy.crs as ccrs
import xarray as xr

from dfm_tools.get_nc import get_netdata, get_ncmodeldata, plot_netmapdata, plot_background
from dfm_tools.xarray_helpers import preprocess_hirlam

dir_testinput = r'c:\DATA\dfm_tools_testdata'
dir_output = '.'

#HIRLAM
file_nc = r'p:\1204257-dcsmzuno\2014\data\meteo\HIRLAM72_2018\h72_201803.nc' #TODO: xarray MissingDimensionsError
data_xr = xr.open_mfdataset(file_nc,drop_variables=['x','y'],preprocess=preprocess_hirlam)
timestep = 0
coarsefac = 2 #coarsen dataset for more performance, but not necessary

data_u = data_xr['eastward_wind'].isel(time=timestep)
data_v = data_xr['northward_wind'].isel(time=timestep)
magn = np.sqrt(data_u**2 + data_v**2)
magn = magn[::coarsefac,::coarsefac] #coarsening makes coordinate conversion faster

fig, ax = plt.subplots()
ax.pcolormesh(magn)
plt.savefig(os.path.join(dir_output,'cartopy_hirlam_raw'))

fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()})
ax.pcolormesh(magn.longitude,magn.latitude,magn)#, transform=ccrs.PlateCarree())
plt.savefig(os.path.join(dir_output,'cartopy_hirlam_aspect'))

fig, ax = plt.subplots(figsize=(9,5),subplot_kw={'projection': ccrs.PlateCarree()}) #provide axis projection on initialisation, cannot be edited later on
pc = ax.pcolormesh(magn.longitude,magn.latitude,magn)#, transform=ccrs.PlateCarree())
cbar = fig.colorbar(pc, ax=ax)
cbar.set_label('velocity magnitude (%s)'%(data_u.attrs['units']))
plot_background(ax=ax, resolution=1, google_style='street', features=['countries_highres'], linewidth=0.5, edgecolor='gray', facecolor='none', latlon_format=True)
plot_background(ax=ax, google_style=None, features=['coastlines_highres'], linewidth=0.5, latlon_format=True)
plt.savefig(os.path.join(dir_output,'cartopy_hirlam_moreoptions'))

fig, ax = plt.subplots(figsize=(6,7),subplot_kw={'projection': ccrs.EuroPP()}) #provide axis projection on initialisation, cannot be edited later on
pc = ax.pcolormesh(magn.longitude, magn.latitude, magn, transform=ccrs.PlateCarree())
plot_background(ax=ax, google_style=None, features=['coastlines_highres'], latlon_format=True, gridlines=True)
plt.savefig(os.path.join(dir_output,'cartopy_hirlam_curvedgridlines'))


#GREVELINGEN
file_nc_map = os.path.join(dir_testinput,'DFM_3D_z_Grevelingen\\computations\\run01\\DFM_OUTPUT_Grevelingen-FM\\Grevelingen-FM_0000_map.nc')
ugrid = get_netdata(file_nc=file_nc_map)
data_frommap_bl = get_ncmodeldata(file_nc=file_nc_map, varname='mesh2d_flowelem_bl')

fig, ax = plt.subplots(1,1, subplot_kw={'projection': ccrs.epsg(28992)}) #provide axis projection on initialisation, cannot be edited later on
pc = plot_netmapdata(ugrid.verts, values=data_frommap_bl, ax=ax, linewidth=0.5, cmap='jet')
plot_background(ax=ax, resolution=12, features=['coastlines_highres'], linewidth=0.5)
plt.savefig(os.path.join(dir_output,'cartopy_grevelingen_RD'))