# -*- coding: utf-8 -*-
"""
Created on Fri Oct 28 06:40:03 2022

@author: veenstra
"""


import os
import matplotlib.pyplot as plt
plt.close('all')
import datetime as dt
import dfm_tools as dfmt

overwrite = True

dir_testinput = r'c:\DATA\dfm_tools_testdata'
file_nc = os.path.join(dir_testinput,r'DFM_3D_z_Grevelingen\computations\run01\DFM_OUTPUT_Grevelingen-FM\Grevelingen-FM_0000_map.nc')

dir_output = os.path.join(os.path.dirname(file_nc),'monthlymeans')
if not os.path.exists(dir_output):
    os.mkdir(dir_output)
postfix = '_'.join(os.path.basename(file_nc).split('_')[-2:])
file_out = os.path.join(dir_output,os.path.basename(file_nc).replace(postfix,f'monthlymean_{postfix}'))

print('opening dataset')
dtstart = dt.datetime.now()
data_frommap_merged = dfmt.open_partitioned_dataset(file_nc.replace('_0000_','_0*_'))

#for varname in data_frommap_merged.data_vars: #TODO: this should be avoided, can also be done with mask_and_scale argument, but is that desireable?
#    if (data_frommap_merged[varname].encoding['dtype']==int) & (data_frommap_merged[varname].dtype!=int):
#        breakit
#        print(f'converting {varname} back to int dtype')
#        data_frommap_merged[varname] = data_frommap_merged[varname].astype(int)
data_xr_var = data_frommap_merged['mesh2d_sa1']
time_passed = (dt.datetime.now()-dtstart).total_seconds()
print(f'>>time passed: {time_passed:.2f} sec')

print('computing monthly means')
dtstart = dt.datetime.now()
if True: #on unique month numbers
    #data_xr_monthmean = data_xr_var.groupby('time.month').mean()
    data_xr_monthmean = data_frommap_merged.groupby('time.month').mean() 
    data_xr_monthmean.rename({'month':'time'})
    #data_xr_monthmean['mesh2d'] = data_frommap_merged.mesh2d #TODO: mesh2d etc are dropped by xugrid, how to arrange this (and convert back to normal mapfile)
else: # on unique year+month combinations
    data_xr_monthmean = data_frommap_merged.copy()
    #data_xr_monthmean = data_xr_var.copy()
    data_xr_monthmean = data_xr_monthmean.resample(time='MS').mean(dim='time')
#data_xr_monthmean.time.encoding = data_xr.time.encoding
time_passed = (dt.datetime.now()-dtstart).total_seconds()
print(f'>>time passed: {time_passed:.2f} sec')

# for varname in data_frommap_merged.data_vars: #TODO: avoid adding time dimension to time-independent variables: https://github.com/pydata/xarray/issues/2145
#     if data_frommap_merged[varname].dims != data_xr_monthmean[varname].dims:
#         #print(f'{varname} dims changed from {data_xr[varname].dims} to {data_xr_monthmean[varname].dims}')
#         data_xr_monthmean[varname] = data_xr_monthmean[varname].isel(time=0) #selecting first time of variables that did not have time dimension before
#         if data_frommap_merged[varname].dims != data_xr_monthmean[varname].dims:
#             print(f'{varname} dims changed from {data_frommap_merged[varname].dims} to {data_xr_monthmean[varname].dims}')

if overwrite:
    print('writing file')
    dtstart = dt.datetime.now()
    data_xr_monthmean.to_netcdf(file_out)
    time_passed = (dt.datetime.now()-dtstart).total_seconds()
    print(f'>>time passed: {time_passed:.2f} sec')


timestep = 0
layer = 33
clim_bl = None
clim_wl = [-0.5,1]
clim_sal = [28,30.2]
clim_tem = [4,10]
crs = "EPSG:28992"

data_frommap = dfmt.open_partitioned_dataset(file_out) #TODO: file_out is not suitable for xugrid anymore, so how to read it? (probably corrupt) "AttributeError: Can only access grid topology via `.grid` if dataset contains exactly one grid. Dataset contains 0 grids. Use `.grids` instead."


print('plot grid and values from mapdata (salinity on layer, 3dim, on cell centers)')
fig, ax = plt.subplots()
pc = data_frommap.isel(time=timestep,layer=layer).ugrid.plot(linewidth=0.5, cmap="jet", edgecolor='face')
pc.set_clim(clim_sal)
# cbar = fig.colorbar(pc, ax=ax)
# cbar.set_label('%s [%s]'%(data_frommap.var_ncattrs['long_name'], data_frommap.var_ncattrs['units']))
# #ax.set_xlabel('x [%s]'%(gridunits))
# #ax.set_ylabel('y [%s]'%(gridunits))
ax.set_aspect('equal')
fig.tight_layout()
#fig.savefig(os.path.join(dir_output,'%s_mesh2d_sa1'%(os.path.basename(file_nc).replace('.',''))))

