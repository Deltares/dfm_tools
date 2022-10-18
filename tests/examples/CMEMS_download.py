# -*- coding: utf-8 -*-
"""
Created on Mon Oct 17 16:55:46 2022

@author: veenstra
"""

import os
from pathlib import Path
import xarray as xr
from dfm_tools.download import download_CMEMS

username = 'slaan' #TODO: remove credentials
password = '"D3lt4r3S!"'# keep double quotation marks for passwords with special characters

#make the /data/tmp directory if it does not exist
dir_output = './cmems_temp'
Path(dir_output).mkdir(parents=True, exist_ok=True)


download_CMEMS(username=username, password=password, #register at: https://resources.marine.copernicus.eu/registration-form'
               dir_output=dir_output,
               longitude_min=2, longitude_max=4, latitude_min=50, latitude_max=52,
               date_min='2010-01-01', date_max='2010-01-02', #'%Y-%m-%d'
               varlist=['bottomT'], #['thetao','so','zos','bottomT','uo','vo'], ['o2','no3','po4','si','nppv','chl'],
               reanalysis_forecast='reanalysis', quantitygroup='physchem', #or provide motu_url/service/product arguments
               #motu_url='http://my.cmems-du.eu/motu-web/Motu', service='GLOBAL_MULTIYEAR_PHY_001_030-TDS', product='cmems_mod_glo_phy_my_0.083_P1D-m',
               )


#open mfdataset to check folder contents
ds = xr.open_mfdataset(os.path.join(dir_output,f'cmems_*.nc'))#, combine='by_coords', decode_times=False)
ds.close()
