# -*- coding: utf-8 -*-
"""
Created on Wed Oct 27 14:35:12 2021

@author: veenstra
"""
import os

import matplotlib.pyplot as plt
plt.close('all')
import numpy as np

from dfm_tools.testutils import try_importmodule
try_importmodule(modulename='cartopy') #check if cartopy was installed since it is an optional module, also happens in plot_cartopybasemap()
import cartopy.crs as ccrs

from dfm_tools.get_nc import get_netdata, get_ncmodeldata, plot_netmapdata, plot_background
from dfm_tools.get_nc_helpers import get_ncvardimlist

dir_testinput = r'c:\DATA\dfm_tools_testdata'
dir_output = '.'
"""
#HIRLAM
file_nc = r'p:\1204257-dcsmzuno\2014\data\meteo\HIRLAM72_2018\h72_201803.nc' #TODO: xarray MissingDimensionsError
vars_pd, dims_pd = get_ncvardimlist(file_nc=file_nc)

timestep = 0
mesh2d_node_x = get_ncmodeldata(file_nc=file_nc, varname='x')
mesh2d_node_y = get_ncmodeldata(file_nc=file_nc, varname='y')
mesh2d_node_x_sel = mesh2d_node_x[::2,::2]
mesh2d_node_y_sel = mesh2d_node_y[::2,::2]
data_v = get_ncmodeldata(file_nc=file_nc, varname='northward_wind',timestep=timestep)
data_u = get_ncmodeldata(file_nc=file_nc, varname='eastward_wind',timestep=timestep)
#airp = get_ncmodeldata(file_nc=file_nc, varname='air_pressure_fixed_height',timestep=0)[0,:,:]
magn = np.sqrt(data_u**2 + data_v**2)[0,::2,::2]

fig, ax = plt.subplots()
ax.pcolor(mesh2d_node_x_sel,mesh2d_node_y_sel,magn)
plt.savefig(os.path.join(dir_output,'cartopy_hirlam_raw'))

fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()})
ax.pcolor(mesh2d_node_x_sel,mesh2d_node_y_sel,magn)#, transform=ccrs.PlateCarree())
plt.savefig(os.path.join(dir_output,'cartopy_hirlam_aspect'))

fig, ax = plt.subplots(figsize=(9,5),subplot_kw={'projection': ccrs.PlateCarree()}) #provide axis projection on initialisation, cannot be edited later on
pc = ax.pcolor(mesh2d_node_x_sel,mesh2d_node_y_sel,magn)#, transform=ccrs.PlateCarree())
cbar = fig.colorbar(pc, ax=ax)
cbar.set_label('velocity magnitude (%s)'%(data_v.var_ncattrs['units']))
plot_background(ax=ax, resolution=1, google_style='street', features=['countries_highres'], linewidth=0.5, edgecolor='gray', facecolor='none', latlon_format=True)
plot_background(ax=ax, google_style=None, features=['coastlines_highres'], linewidth=0.5, latlon_format=True)
plt.savefig(os.path.join(dir_output,'cartopy_hirlam_moreoptions'))

fig, ax = plt.subplots(figsize=(6,7),subplot_kw={'projection': ccrs.EuroPP()}) #provide axis projection on initialisation, cannot be edited later on
pc = ax.pcolor(mesh2d_node_x_sel[:100,:100],mesh2d_node_y_sel[:100,:100],magn[:100,:100], transform=ccrs.PlateCarree()) #take subset of dataset to speed up coordinate transformation
plot_background(ax=ax, google_style=None, features=['coastlines_highres'], latlon_format=True, gridlines=True)
plt.savefig(os.path.join(dir_output,'cartopy_hirlam_curvedgridlines'))
"""
    
#GREVELINGEN
file_nc_map = os.path.join(dir_testinput,'DFM_3D_z_Grevelingen\\computations\\run01\\DFM_OUTPUT_Grevelingen-FM\\Grevelingen-FM_0000_map.nc')
ugrid = get_netdata(file_nc=file_nc_map)
data_frommap_bl = get_ncmodeldata(file_nc=file_nc_map, varname='mesh2d_flowelem_bl')

fig, ax = plt.subplots(1,1, subplot_kw={'projection': ccrs.epsg(28992)}) #provide axis projection on initialisation, cannot be edited later on
pc = plot_netmapdata(ugrid.verts, values=data_frommap_bl, ax=ax, linewidth=0.5, cmap='jet')
plot_background(ax=ax, resolution=12, features=['coastlines_highres'], linewidth=0.5)
plt.savefig(os.path.join(dir_output,'cartopy_grevelingen_RD'))