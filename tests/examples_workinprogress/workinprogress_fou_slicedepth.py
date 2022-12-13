# -*- coding: utf-8 -*-
"""
Created on Tue Dec 13 22:00:18 2022

@author: veenstra
"""

import os
import matplotlib.pyplot as plt
plt.close('all')
import dfm_tools as dfmt

dir_testinput = r'c:\DATA\dfm_tools_testdata'
dir_output = '.'

file_nc = r'p:\archivedprojects\11203379-005-mwra-updated-bem\03_model\02_final\A72_ntsu0_kzlb2\DFM_OUTPUT_MB_02_fou\MB_02_0*_fou.nc'

timestep = 10
layno = 45
val_ylim = [-600,1]
clim_bl = [-500,0]
clim_sal = [25,36]
crs = "EPSG:4326"
file_nc_fou = None

data_frommap_merged = dfmt.open_partitioned_dataset(file_nc)

vars_pd = dfmt.get_ncvarproperties(data_frommap_merged)

print('plot grid and values from mapdata (salinity on layer, 3dim, on cell centers), on layer')
fig, ax = plt.subplots()
if 'nmesh2d_layer' in data_frommap_merged['mesh2d_fourier001_mean'].dims: #use argument missing_dims='ignore' instead
    pc = data_frommap_merged['mesh2d_fourier001_mean'].isel(nmesh2d_layer=layno).ugrid.plot(edgecolor='face',cmap='jet')
else:
    pc = data_frommap_merged['mesh2d_fourier001_mean'].ugrid.plot(edgecolor='face',cmap='jet')
pc.set_clim(clim_sal)
ax.set_aspect('equal')
fig.tight_layout()


#TODO: the below part does not work yet, since mesh2d_s1 and mesh2d_bl vars are missing
print('plot grid and values from mapdata (salinity on layer, 3dim, on cell centers), on fixed depth(s)')
depths = [-1,-4]
data_frommap_timesel = data_frommap_merged #select data for all layers
data_frommap_timesel_atdepths = dfmt.get_Dataset_atdepths(data_xr=data_frommap_timesel, depths=depths, reference='z0') #depth w.r.t. z0/waterlevel/bedlevel
for dep in depths:
    fig, ax = plt.subplots()
    if 'depth_fromref' in data_frommap_timesel_atdepths.dims: #TODO: use missingdims=ignore so if-statement is not necessary
        pc = data_frommap_timesel_atdepths['mesh2d_fourier001_mean'].sel(depth_fromref=dep).ugrid.plot(edgecolor='face',cmap='jet')
    else:
        pc = data_frommap_timesel_atdepths['mesh2d_fourier001_mean'].ugrid.plot(edgecolor='face',cmap='jet')
    pc.set_clim(clim_sal)
    ax.set_aspect('equal')
    fig.tight_layout()
    