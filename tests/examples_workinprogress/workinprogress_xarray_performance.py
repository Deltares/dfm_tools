# -*- coding: utf-8 -*-
"""
Created on Tue Apr  4 11:41:56 2023

@author: veenstra

module load netcdf/v4.9.0_v4.6.0_intel22.2.0
nccopy -k 'netCDF-4' 'DCSM-FM_0_5nm_0000_his.nc' 'DCSM-FM_0_5nm_0000_his_netcdf4.nc'
"""

import datetime as dt
import dfm_tools as dfmt
import xarray as xr
from dask.diagnostics import ProgressBar


file_nc = r'P:\archivedprojects\11208054-004-dcsm-fm\models\3D_DCSM-FM\2013-2017\B04_EOT20_RHO1_H1_H2\DFM_OUTPUT_DCSM-FM_0_5nm\DCSM-FM_0_5nm_0000_his.nc' #does not seem to be dask? >> no chunks visible
#file_nc = r'P:\archivedprojects\11208054-004-dcsm-fm\models\3D_DCSM-FM\2013-2017\B04_EOT20_RHO1_H1_H2\DFM_OUTPUT_DCSM-FM_0_5nm\DCSM-FM_0_5nm_0000_his_netcdf4.nc' #does not seem to be dask? >> no chunks visible


chunks = {}
chunks = {'time':1,'stations':10}
chunks = {'time':1000,'stations':10}
chunks = 'auto'

#performance measurements can be influenced by partly caching in memory
#netcdf3
#Open: 1581 sec (26 min) #Plot: 12646 sec (210 min)
#Open: 864  sec (14 min) #Plot: 2422 sec (40 min)
#Open: 705  sec (12 min) #Plot: 1323 sec (22 min)
#netcdf4
#Open: 0.3 sec (0 min) #Plot: ?? sec (?? min) #killed after 33 minutes of plotting (0% progress)
#Open: 20  sec (0 min) #Plot: 7505 sec (125 min)
#Open: 0.3 sec (0 min) #Plot: 3405 sec (56 min)


print('>> performance test opening: ',end='')
dtstart = dt.datetime.now()
ds = xr.open_dataset(file_nc,chunks=chunks)#,decode_cf=False,decode_times=False,decode_coords=False)
print(f'{(dt.datetime.now()-dtstart).total_seconds():.2f} sec')

ds = dfmt.preprocess_hisnc(ds)

print('>> performance test plotting: ',end='')
dtstart = dt.datetime.now()
ds_toplot = ds.salinity.isel(laydim=39).sel(stations='HOEKVHLD')
with ProgressBar(): #contains not all time needed
    ds_toplot.plot()
print(f'{(dt.datetime.now()-dtstart).total_seconds():.2f} sec')

