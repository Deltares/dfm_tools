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

#TODO: experiment with monthly/daily means or depth average of his/map fields (to show convenience of pandas/xarray)

overwrite = True

file_nc = dfmt.data.fm_grevelingen_map(return_filepath=True)
file_basename = os.path.basename(file_nc.replace('_0*_','_0000_'))

dir_output = '.'
postfix = '_'.join(file_basename.split('_')[-2:])
file_out = os.path.join(dir_output,file_basename.replace(postfix,f'monthlymean_{postfix}'))

data_frommap_merged = dfmt.open_partitioned_dataset(file_nc)
data_frommap_merged = dfmt.Dataset_varswithdim(data_frommap_merged, dimname='time') #drop all variables without time dimension

print('>> computing monthly means: ', end='')
dtstart = dt.datetime.now()
data_xr_monthmean = data_frommap_merged.resample(time='MS').mean(dim='time')
# on unique month numbers
# data_xr_monthmean = data_frommap_merged.groupby('time.month').mean() 
# data_xr_monthmean = data_xr_monthmean.rename({'month':'time'})
print(f'{(dt.datetime.now()-dtstart).total_seconds():.2f} sec')

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
pc = data_frommap['mesh2d_sa1'].isel(nmesh2d_layer=33,month=0,time=0,missing_dims='ignore').ugrid.plot(linewidth=0.5, cmap="jet")
pc.set_clim([28,30.2])
ax.set_aspect('equal')
fig.tight_layout()
#fig.savefig(os.path.join(dir_output,'%s_mesh2d_sa1'%(file_basename.replace('.',''))))

