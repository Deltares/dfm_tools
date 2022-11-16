# -*- coding: utf-8 -*-
"""
Created on Mon Nov 14 14:34:12 2022

@author: veenstra
"""

import matplotlib.pyplot as plt
plt.close('all')
import xarray as xr
import dfm_tools as dfmt

file_nc_bc = r'p:\1204257-dcsmzuno\data\CMEMS\bnd\NorthSeaAndBaltic_1993-2019_20210510\temperature_extra_rand_dcsm.nc'

data_xr = xr.open_dataset(file_nc_bc)
data_xr = data_xr.set_coords(['station_id'])#,'station_names'])
data_xr = dfmt.preprocess_hisnc(data_xr)

fig,ax = plt.subplots()
data_xr.temperature.sel(node='extra_rand_dcsm_0001').T.plot()
fig.tight_layout()

