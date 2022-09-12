# -*- coding: utf-8 -*-
"""
Created on Wed Oct 27 14:47:34 2021

@author: veenstra

this test retrieves sobek observation data and plots it
"""

import os
import matplotlib.pyplot as plt
plt.close('all')
import xarray as xr

#from dfm_tools.get_nc import get_ncmodeldata
from dfm_tools.get_nc_helpers import get_stationid_fromstationlist#, get_hisstationlist

dir_testinput = r'c:\DATA\dfm_tools_testdata'
dir_output = '.'

file_nc = os.path.join(dir_testinput,'KenmerkendeWaarden','observations.nc')

#station_names = get_hisstationlist(file_nc=file_nc, varname='observation_id')
#station_names = get_hisstationlist(file_nc=file_nc, varname='water_level')
#data_fromsobek = get_ncmodeldata(file_nc=file_nc, varname='water_level', station=['Maasmond','HKVHLD','MO_1035.00'], timestep='all')

data_xr = xr.open_dataset(file_nc)
#stations_pd = data_xr.observation_id.astype(str).to_pandas()
station_list = ['Maasmond','HKVHLD','MO_1035.00']
idx_stations = get_stationid_fromstationlist(data_xr, stationlist=station_list, station_varname='observation_id')

fig, ax = plt.subplots()
for iS,stat_name in zip(idx_stations,station_list):
    data_fromsobek = data_xr.water_level.isel(id=iS)
    ax.plot(data_fromsobek.time,data_fromsobek,'-', label=stat_name)
ax.legend()
plt.savefig(os.path.join(dir_output,'%s_waterlevel'%(os.path.basename(file_nc).replace('.',''))))




    
    
