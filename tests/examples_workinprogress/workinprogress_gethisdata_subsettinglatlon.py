# -*- coding: utf-8 -*-
"""
Created on Mon Oct 17 15:40:49 2022

@author: veenstra
"""

from pathlib import Path
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
plt.close('all')
from dfm_tools.hydrolib_helpers import pointlike_to_DataFrame #install dfm_tools via https://github.com/openearth/dfm_tools
#from dfm_tools.get_nc_helpers import get_hisstationlist, get_stationid_fromstationlist
from hydrolib.core.io.polyfile.models import PolyFile #automatically installed with dfm_tools

file_pli = r'p:\1230882-emodnet_hrsm\GTSMv5.0\runs\reference_GTSMv4.1_wiCA\world.ldb'
polyfile_object = PolyFile(Path(file_pli))
polyobject_pd = pointlike_to_DataFrame(polyfile_object.objects[0])
polyobject_pd[polyobject_pd==999.999] = np.nan

file_nc = r'p:\1230882-emodnet_hrsm\GTSMv5.0\runs\reference_GTSMv4.1_wiCA\output\gtsm_model_0000_his.nc'
#file_nc = r'p:\1230882-emodnet_hrsm\GTSMv5.0\runs\GM42_2000m_eu0900m_ITfac5p5_wx\output\gtsm_model_0000_his.nc'
data_xr = xr.open_dataset(file_nc,chunks={'time':-1})

data_xr['station_name_str'] = data_xr['station_name'].to_pandas().str.decode('utf-8',errors='ignore').str.strip() #to_pandas is essential otherwise resulting bool might not be correct. .str.strip() to remove spaces left/right from station_name (necessary for sobek models)
data_xr = data_xr.set_coords('station_name_str')

data_xr_relvars = data_xr[['waterlevel','x_velocity','y_velocity']]

latlon_bool = ((data_xr['station_x_coordinate']>-10) & (data_xr['station_x_coordinate']<-8) & 
               (data_xr['station_y_coordinate']>40) & (data_xr['station_y_coordinate']<43))


data_xr_relvars_port = data_xr_relvars.sel(stations=latlon_bool).sel(time=slice('2014-01-01','2014-01-15')).load()

data_xr_relvars_port['velocity_direction'] = np.arctan2(data_xr_relvars_port['y_velocity'],data_xr_relvars_port['x_velocity'])
data_xr_relvars_port['velocity_direction'].attrs['standard_name'] = 'velocity_direction'
data_xr_relvars_port['velocity_direction'].attrs['long_name'] = 'velocity_direction'
data_xr_relvars_port['velocity_direction'].attrs['units'] = 'degrees'
data_xr_relvars_port['velocity_magnitude'] = np.sqrt(data_xr_relvars_port['y_velocity']**2+data_xr_relvars_port['x_velocity']**2)
data_xr_relvars_port['velocity_magnitude'].attrs['units'] = 'm/s'
#stations_pd = get_hisstationlist(file_nc)
#stat_port_list = data_xr_relvars_port.station_name_str.data

fig,(ax1,ax2) = plt.subplots(2,1,figsize=(10,7),sharex=True)
data_xr_relvars_port['velocity_magnitude'].plot.line(x='time',ax=ax1,linewidth=0.8)
data_xr_relvars_port['velocity_direction'].plot.line(x='time',ax=ax2,linewidth=0.8)
ax1.legend(data_xr_relvars_port['station_name_str'].to_series(),fontsize=7,loc=1)
ax2.legend(data_xr_relvars_port['station_name_str'].to_series(),fontsize=7,loc=1)
fig.tight_layout()

figmap,axmap = plt.subplots()
axmap.plot(polyobject_pd['x'],polyobject_pd['y'],'k-',linewidth=0.8)
axmap.plot(data_xr_relvars_port['station_x_coordinate'],data_xr_relvars_port['station_y_coordinate'],'xr')
axmap.set_xlim(-15,5)
axmap.set_ylim(34,47)
