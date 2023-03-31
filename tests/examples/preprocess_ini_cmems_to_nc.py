# -*- coding: utf-8 -*-
"""
Created on Tue Sep 27 12:14:14 2022

@author: veenstra
"""

import os
import glob
import datetime as dt
import xarray as xr

#TODO: merge with other ini script and make generic for getting an inifield out of CMEMS/etc regulargrid Dataset or a 2D/3D FM map/rst Dataset

tSimStart = dt.datetime(1998,1,1)
dir_data  = r'p:\i1000668-tcoms\03_newModel\01_input\02_bnd\data_opendap' #folder containing CMEMS so and thetao netcdf files

dir_out = '.'

file_nc_list_so = glob.glob(f'{dir_data}\\cmems_so_199*.nc')
file_nc_list_thetao = glob.glob(f'{dir_data}\\cmems_thetao_199*.nc')
file_nc_list = file_nc_list_so + file_nc_list_thetao

print(f'opening {len(file_nc_list)} datasets')
data_xr = xr.open_mfdataset(file_nc_list)
data_xr_times_pd = data_xr.time.to_series()

interp = False
if interp: #this would be the proper way to do it, but FM needs two timesteps for some reason
    print('ds.interp()')
    data_xr_ontime = data_xr.interp(time=[tSimStart],kwargs=dict(bounds_error=True)) #bounds_error makes sure, outofbounds time results in "ValueError: A value in x_new is below the interpolation range."
else:
    print('ds.sel()')
    data_xr_ontime = data_xr.sel(time=slice(tSimStart-dt.timedelta(days=1),tSimStart+dt.timedelta(days=1)))

print('writing file')
outFile = os.path.join(dir_out,f'InitialField_{tSimStart.strftime("%Y-%m-%d_%H-%M-%S")}.nc')
data_xr_ontime.to_netcdf(outFile,format="NETCDF4_CLASSIC")

