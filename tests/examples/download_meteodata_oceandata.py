# -*- coding: utf-8 -*-
"""
Created on Tue Oct 18 15:20:06 2022

@author: veenstra
"""

import os
import xarray as xr
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
plt.close('all')
from dfm_tools.download import download_ERA5, download_CMEMS

#download ERA5/CMEMS data for given domain, time extent and variable lists
#TODO: add CMCC, GFDL, HYCOM

overwrite = False # always set to True when changing the domain

# domain
longitude_min, longitude_max, latitude_min, latitude_max =    2,   4,  50, 52 #test domain
#longitude_min, longitude_max, latitude_min, latitude_max = -180, 180, -90, 90 #global

#dates as understood by pandas.date_range(). ERA5 has freq='M' (month) and CMEMS has freq='D' (day)
date_min = '2010-01-01'
date_max = '2010-01-02'

#variables per model will be written to separate netcdf files. Set to [] to skip model.
variables_era5 = ['v10n'] # variables, supply arbitrary string to get error with available variable names
varlist_cmems = ['bottomT','thetao']#,'so','zos','uo','vo'], #['o2','no3','po4','si','nppv','chl'],

#output directories per model
dir_output_era5 = './era5_temp'
dir_output_cmems = './cmems_temp'


#ERA5
dir_output = dir_output_era5
for varkey in variables_era5:
    if not os.path.isdir(dir_output):
        os.mkdir(dir_output)
    
    download_ERA5(varkey, 
                  date_min=date_min, date_max=date_max,
                  longitude_min=longitude_min, longitude_max=longitude_max, latitude_min=latitude_min, latitude_max=latitude_max,
                  dir_output=dir_output, overwrite=overwrite)

    #open mfdataset to check folder contents
    ds = xr.open_mfdataset(os.path.join(dir_output,f'era5_{varkey}_*.nc'))
    ds.close()
    fig,ax = plt.subplots()
    ds[varkey].isel(time=0).plot(ax=ax)


#CMEMS
dir_output = dir_output_cmems
date_min_cmems = pd.Timestamp(date_min)-pd.Timedelta(days=1) #CMEMS has daily noon values (not midnight), so subtract one day from date_min to cover desired time extent
for varkey in varlist_cmems:
    Path(dir_output).mkdir(parents=True, exist_ok=True)
    
    download_CMEMS(varkey,
                   date_min=date_min_cmems, date_max=date_max,
                   longitude_min=longitude_min, longitude_max=longitude_max, latitude_min=latitude_min, latitude_max=latitude_max,
                   dir_output=dir_output, overwrite=overwrite,
                   #source_combination='multiyear_physchem', #str from source_dict.keys(). Or provide motu_url/service/product arguments
                   motu_url='http://my.cmems-du.eu', service='GLOBAL_MULTIYEAR_PHY_001_030-TDS', product='cmems_mod_glo_phy_my_0.083_P1D-m', #or provide source_combination argument
                   credentials=None, #credentials=['username','password'], or create "%USER%/motucredentials.txt" with username on line 1 and password on line 2. Register at: https://resources.marine.copernicus.eu/registration-form'
                   timeout=200) #default is 30, but this might be not always be enough to avoid TimeoutExpired error.
    
    #open mfdataset to check folder contents and plot first field of each variable
    ds = xr.open_mfdataset(os.path.join(dir_output,f'cmems_{varkey}_*.nc'))
    ds.close()
    fig,ax = plt.subplots()
    if 'depth' in ds[varkey].dims:
        ds[varkey].isel(time=0,depth=0).plot(ax=ax)
    else:
        ds[varkey].isel(time=0).plot(ax=ax)


