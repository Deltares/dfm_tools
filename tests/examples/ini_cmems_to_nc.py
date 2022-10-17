# -*- coding: utf-8 -*-
"""
Created on Tue Sep 27 12:14:14 2022

@author: veenstra
"""

import os
import glob
import datetime as dt
import xarray as xr

tSimStart = dt.datetime(1998,1,1)
dir_data  = r'p:\i1000668-tcoms\03_newModel\01_input\02_bnd\CMEMS\download\data' #folder containing CMEMS so and thetao netcdf files

dir_out = '.'

file_nc_list_so = glob.glob(f'{dir_data}\\so_*.nc')
file_nc_list_thetao = glob.glob(f'{dir_data}\\thetao_*.nc')
file_nc_list = file_nc_list_so + file_nc_list_thetao

data_xr = xr.open_mfdataset(file_nc_list)
data_xr_ontime = data_xr.interp(time=tSimStart)

outFile = os.path.join(dir_out,f'InitialField_{tSimStart.strftime("%Y-%m-%d_%H-%M-%S")}.nc')

data_xr_ontime.to_netcdf(outFile,format="NETCDF4_CLASSIC")

