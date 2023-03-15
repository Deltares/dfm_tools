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

#TODO: experiment with monthly/daily means or depth average of his/map fields (to show power of pandas/xarray)

overwrite = True

dir_testinput = r'c:\DATA\dfm_tools_testdata'
file_nc = os.path.join(dir_testinput,r'DFM_3D_z_Grevelingen\computations\run01\DFM_OUTPUT_Grevelingen-FM\Grevelingen-FM_0000_map.nc')

dir_output = '.'
postfix = '_'.join(os.path.basename(file_nc).split('_')[-2:])
file_out = os.path.join(dir_output,os.path.basename(file_nc).replace(postfix,f'monthlymean_{postfix}'))

data_frommap_merged = dfmt.open_partitioned_dataset(file_nc.replace('_0000_','_0*_'))

print('>> computing monthly means: ', end='')
dtstart = dt.datetime.now()
if 0: #on unique month numbers
    data_xr_monthmean = data_frommap_merged.groupby('time.month').mean() 
    data_xr_monthmean.rename({'month':'time'})
else: # on unique year+month combinations
    data_xr_monthmean = data_frommap_merged.resample(time='MS').mean(dim='time')
print(f'{(dt.datetime.now()-dtstart).total_seconds():.2f} sec')

for varname in data_frommap_merged.data_vars: #TODO: avoid adding time dimension to time-independent variables: https://github.com/pydata/xarray/issues/2145
    if data_frommap_merged[varname].dims != data_xr_monthmean[varname].dims:
        #print(f'{varname} dims changed from {data_xr[varname].dims} to {data_xr_monthmean[varname].dims}')
        data_xr_monthmean[varname] = data_xr_monthmean[varname].isel(time=0) #selecting first time of variables that did not have time dimension before
        if data_frommap_merged[varname].dims != data_xr_monthmean[varname].dims:
            print(f'{varname} dims changed from {data_frommap_merged[varname].dims} to {data_xr_monthmean[varname].dims}')

#reconnect data and grid #TODO: support resampling/groupby in xugrid?
import xugrid as xu
data_xr_monthmean = xu.UgridDataset(data_xr_monthmean, grids=[data_frommap_merged.grid])

if overwrite:
    print('>> writing file: ', end='')
    dtstart = dt.datetime.now()
    data_xr_monthmean.ugrid.to_netcdf(file_out)
    print(f'{(dt.datetime.now()-dtstart).total_seconds():.2f} sec')



data_frommap = dfmt.open_partitioned_dataset(file_out)
crs = "EPSG:28992"

print('plot grid and values from mapdata (salinity on layer, 3dim, on cell centers)')
fig, ax = plt.subplots()
pc = data_frommap['mesh2d_sa1'].isel(nmesh2d_layer=33,month=0,time=0,missing_dims='ignore').ugrid.plot(linewidth=0.5, cmap="jet", edgecolor='face')
pc.set_clim([28,30.2])
ax.set_aspect('equal')
fig.tight_layout()
#fig.savefig(os.path.join(dir_output,'%s_mesh2d_sa1'%(os.path.basename(file_nc).replace('.',''))))

