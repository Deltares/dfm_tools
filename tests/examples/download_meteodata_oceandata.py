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
import contextily as ctx
import dfm_tools as dfmt

#download ERA5/CMEMS/HYCOM data for given domain, time extent and variables
#TODO: add CMCC, GFD
#TODO: add climatedata cmip6
#TODO: add GFS and other NOAA models (https://www.ncei.noaa.gov/products/weather-climate-models/global-forecast > NCEI > TDS)
#TODO: add click?

overwrite = True # always set to True when changing the domain

# domain
longitude_min, longitude_max, latitude_min, latitude_max =    2,   4,  50, 52 #test domain
#longitude_min, longitude_max, latitude_min, latitude_max = -180, 180, -90, 90 #global

#dates as understood by pandas.period_range(). ERA5 has freq='M' (month) and CMEMS has freq='D' (day)
date_min = '2010-01-01'
date_max = '2010-01-02'

#variables per model will be written to separate netcdf files. Set to [] to skip model.
variables_era5 = []#'v10n'] # supply arbitrary string to get error with available variable names
varlist_cmems = []#['bottomT','thetao','no3'] # avaliable variables differ per source_combination, check cmems loop for some options
varlist_hycom = ['surf_el']#'water_temp'] #['tau','water_u','water_v','water_temp','salinity','surf_el']

#output directories per model
dir_output_era5 = './era5_temp'
dir_output_cmems = './cmems_temp'
dir_output_hycom = './hycom_temp'


#ERA5
dir_output = dir_output_era5
for varkey in variables_era5:
    if not os.path.isdir(dir_output):
        os.mkdir(dir_output)
    
    dfmt.download_ERA5(varkey, 
                       longitude_min=longitude_min, longitude_max=longitude_max, latitude_min=latitude_min, latitude_max=latitude_max,
                       date_min=date_min, date_max=date_max,
                       dir_output=dir_output, overwrite=overwrite)

    #open mfdataset to check folder contents
    ds = xr.open_mfdataset(os.path.join(dir_output,f'era5_{varkey}_*.nc'))
    ds.close()
    fig,ax = plt.subplots()
    ds[varkey].isel(time=0).plot(ax=ax)
    ctx.add_basemap(ax=ax,crs="EPSG:4326",attribution=False)


#CMEMS
dir_output = dir_output_cmems
date_min_cmems = pd.Timestamp(date_min)-pd.Timedelta(days=1) #CMEMS has daily noon values (not midnight), so subtract one day from date_min to cover desired time extent
for varkey in varlist_cmems:
    Path(dir_output).mkdir(parents=True, exist_ok=True)
    
    if varkey in ['bottomT','mlotst','siconc','sithick','so','thetao','uo','usi','vo','vsi','zos']: #for physchem
        dataset_url = 'https://my.cmems-du.eu/thredds/dodsC/cmems_mod_glo_phy_my_0.083_P1D-m'
        #dataset_url = 'https://nrt.cmems-du.eu/thredds/dodsC/global-analysis-forecast-phy-001-024'
    else: # for bio
        dataset_url = 'https://my.cmems-du.eu/thredds/dodsC/cmems_mod_glo_bgc_my_0.25_P1D-m' #contains ['chl','no3','nppv','o2','po4','si']
        #dataset_url = 'https://nrt.cmems-du.eu/thredds/dodsC/global-analysis-forecast-bio-001-028-daily' #contains ['chl','fe','no3','nppv','o2','ph','phyc','po4','si','spco2']
    file_prefix = 'cmems_'
    
    dfmt.download_OPeNDAP(dataset_url=dataset_url,
                          credentials=None, #credentials=['username','password'], or create "%USER%/motucredentials.txt" with username on line 1 and password on line 2. Register at: https://resources.marine.copernicus.eu/registration-form'
                          varkey=varkey,
                          longitude_min=longitude_min, longitude_max=longitude_max, latitude_min=latitude_min, latitude_max=latitude_max,
                          date_min=date_min_cmems, date_max=date_max,
                          dir_output=dir_output, file_prefix=file_prefix, overwrite=overwrite)
    
    #open mfdataset to check folder contents and plot first field of each variable
    ds = xr.open_mfdataset(os.path.join(dir_output,f'{file_prefix}{varkey}_*.nc'))
    fig,ax = plt.subplots()
    if 'depth' in ds[varkey].dims:
        ds[varkey].isel(time=0,depth=0).plot(ax=ax)
    else:
        ds[varkey].isel(time=0).plot(ax=ax)
    ds.close()
    ctx.add_basemap(ax=ax,crs="EPSG:4326",attribution=False)
    

#HYCOM
dir_output = dir_output_hycom
for varkey in varlist_hycom:
    Path(dir_output).mkdir(parents=True, exist_ok=True)
    
    period_range_years = pd.period_range(date_min,date_max,freq='Y')
    dataset_url = [f'https://tds.hycom.org/thredds/dodsC/GLBu0.08/expt_19.1/{year}' for year in period_range_years] #list is possible with hycom, since it uses xr.open_mfdataset()
    file_prefix = 'hycom_'
    
    dfmt.download_OPeNDAP(dataset_url=dataset_url,
                          varkey=varkey,
                          longitude_min=longitude_min, longitude_max=longitude_max, latitude_min=latitude_min, latitude_max=latitude_max,
                          date_min=date_min, date_max=date_max,
                          dir_output=dir_output, file_prefix=file_prefix, overwrite=overwrite)
    
    #open mfdataset to check folder contents and plot first field of each variable
    ds = xr.open_mfdataset(os.path.join(dir_output,f'{file_prefix}{varkey}_*.nc'))
    fig,ax = plt.subplots()
    if 'depth' in ds[varkey].dims:
        ds[varkey].isel(time=0,depth=0).plot(ax=ax)
    else:
        ds[varkey].isel(time=0).plot(ax=ax)
    ds.close()
    ctx.add_basemap(ax=ax,crs="EPSG:4326",attribution=False)


"""
Download CMCC data (possible with opendap? urls seem really specific)
userguide: https://esgf.github.io/esgf-user-support/user_guide.html
server: https://esgf-data.dkrz.de/search/cmip6-dkrz/ (Sanne: exacte kopie op https://esgf-node.llnl.gov/search/cmip6/)

Example selection:
source-id: CMCC-ESM2
experiment-id: historical
table-id: Omon (=monthly)
Variable: no3/o2

> Search
"""

"""
CMIP6
https://esgf-node.llnl.gov/search/cmip6/
there is also a esgf python client: https://esgf-pyclient.readthedocs.io/en/latest/
notebook: https://esgf-pyclient.readthedocs.io/en/latest/notebooks/examples/download.html

example cdsapi request
request_dict = {
                'format': 'zip',
                'experiment': 'historical',
                'variable': 'storm_surge_residual',
                'model': 'CMCC-CM2-VHR4',
                'year': '1982',
                'month': '04',
                'temporal_aggregation': '10_min',
                }

file_out = 'download.zip'
c.retrieve(name='sis-water-level-change-timeseries-cmip6', request=request_dict, target=file_out)

"""
