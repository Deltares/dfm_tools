# -*- coding: utf-8 -*-
"""
Created on Wed Apr 26 12:25:27 2023

@author: veenstra
"""

import os
import xarray as xr
import matplotlib.pyplot as plt
plt.close('all')
import numpy as np
import dfm_tools as dfmt
import xugrid as xu
import datetime as dt
from dask.diagnostics import ProgressBar

ds_gebco = xr.open_dataset('p:\\metocean-data\\open\\GEBCO\\2022\\GEBCO_2022.nc')
#ds_gebco = ds_gebco.where(ds_gebco.elevation<100) #mask all land values
file_net = r'p:\1230882-emodnet_hrsm\global_tide_surge_model\trunk\gtsm4.1\step11_global_1p25eu_withcellinfo_net.nc'

uds = xu.open_dataset(file_net)

nnodes = uds.dims[uds.grid.node_dimension]

step = 5
print(f'interpolating GEBCO to {nnodes} nodes in 360/{step}={360/step} steps:')
dtstart = dt.datetime.now()
uds_part_list = []
for i in range(-180, 180+1, step):
    xslice = slice(i,i+step)
    print(xslice)
    
    uds_xsel = uds.ugrid.sel(x=xslice)
    x_sel, y_sel = uds_xsel.obj.NetNode_x, uds_xsel.obj.NetNode_y
    
    z_sel = ds_gebco.interp(lon=x_sel, lat=y_sel).reset_coords(['lat','lon']) #interpolates lon/lat gebcodata to mesh2d_nNodes dimension #TODO: if these come from xu_grid_uds (without ojb), the mesh2d_node_z var has no ugrid accessor since the dims are lat/lon instead of mesh2d_nNodes
    uds_xsel['NetNode_z'] = z_sel.elevation
    uds_part_list.append(uds_xsel)
    # if i==-170:
    #     breakit

print('merging slices')
uds_withz = xu.merge_partitions(uds_part_list)
print('save grid')
#uds_withz.NetNode_z.ugrid.plot()
print(f'{(dt.datetime.now()-dtstart).total_seconds():.2f} sec')


#xu_grid_uds['mesh2d_node_z'] = data_bathy_interp.clip(max=10)



