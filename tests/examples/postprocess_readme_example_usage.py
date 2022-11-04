# -*- coding: utf-8 -*-
"""
Created on Wed Oct 27 14:55:27 2021

@author: veenstra
"""

import os
import xarray as xr
import matplotlib.pyplot as plt
plt.close('all')
from dfm_tools.get_nc import get_netdata, get_ncmodeldata, plot_netmapdata
from dfm_tools.xarray_helpers import preprocess_hisnc

dir_testinput = os.path.join(r'n:\Deltabox\Bulletin\veenstra\info dfm_tools\test_input')
file_nc_map = os.path.join(dir_testinput,'DFM_sigma_curved_bend','DFM_OUTPUT_cb_3d','cb_3d_map.nc')
file_nc_his = os.path.join(dir_testinput,'DFM_sigma_curved_bend','DFM_OUTPUT_cb_3d','cb_3d_his.nc')

data_xr_his = xr.open_mfdataset(file_nc_his, preprocess=preprocess_hisnc)
stations_pd = data_xr_his['stations'].to_dataframe()

#retrieve his data and plot
fig, ax = plt.subplots(1,1,figsize=(10,5))
data_xr_his.waterlevel.plot.line(ax=ax, x='time')
ax.legend(data_xr_his.stations.to_series(),loc=1) #optional, to change legend location
fig.tight_layout()

#plot net/grid
ugrid_all = get_netdata(file_nc=file_nc_map)#,multipart=False)
fig, ax = plt.subplots()
pc = plot_netmapdata(ugrid_all.verts, values=None, ax=None, linewidth=0.5, color="crimson", facecolor="None")
ax.set_aspect('equal')

#plot water level on map
data_frommap_wl = get_ncmodeldata(file_nc=file_nc_map, varname='mesh2d_s1', timestep=3)#, multipart=False)
fig, ax = plt.subplots()
pc = plot_netmapdata(ugrid_all.verts, values=data_frommap_wl[0,:], ax=None, linewidth=0.5, cmap="jet")
pc.set_clim([-0.5,1])
fig.colorbar(pc, ax=ax)
ax.set_title('%s (%s)'%(data_frommap_wl.var_varname, data_frommap_wl.var_ncattrs['units']))
ax.set_aspect('equal')

#plot salinity on map
data_frommap_sal = get_ncmodeldata(file_nc=file_nc_map, varname='mesh2d_sa1', timestep=2, layer=5)#, multipart=False)
fig, ax = plt.subplots()
pc = plot_netmapdata(ugrid_all.verts, values=data_frommap_sal[0,:,0], ax=None, linewidth=0.5, cmap="jet")
fig.colorbar(pc, ax=ax)
ax.set_title('%s (%s)'%(data_frommap_sal.var_varname, data_frommap_sal.var_ncattrs['units']))
ax.set_aspect('equal')

#print contents of retrieved data withing data_frommap_sal variable
print_var = data_frommap_sal
print('++++++\nthe data in the variable %s is:\n%s\n'%(print_var.var_varname, print_var))
print('++++++\nthe time indices and times in the variable %s are:\n%s\n'%(print_var.var_varname, print_var.var_times))
print('++++++\nthe layer indices in the variable %s are:\n%s\n'%(print_var.var_varname, print_var.var_layers))
print('++++++\nthe shape of the variable %s is:\n%s\n'%(print_var.var_varname, print_var.shape))
print('++++++\nthe dimensions of the variable %s are (copied from netCDF variable):\n%s\n'%(print_var.var_varname, print_var.var_dimensions))
print('++++++\nthe netCDF variable where the data in variable %s comes from is:\n%s\n'%(print_var.var_varname, print_var.var_ncvarobject))
print('++++++\nsome example contents of this netCDF variable:')
print('\tthe dimension names of the netCDF variable %s are:\n\t\t%s'%(print_var.var_varname, print_var.var_dimensions))
print('\tthe shape of the netCDF variable %s is:\n\t\t%s'%(print_var.var_varname, print_var.var_shape))
print('\tthe units of the netCDF variable %s are:\n\t\t%s'%(print_var.var_varname, print_var.var_ncattrs['units']))
print('\tthe long_name of the netCDF variable %s is:\n\t\t%s'%(print_var.var_varname, print_var.var_ncattrs['long_name']))
print('\tthe standard_name of the netCDF variable %s is:\n\t\t%s'%(print_var.var_varname, print_var.var_ncattrs['standard_name']))

