# -*- coding: utf-8 -*-
"""
Created on Mon Oct 24 22:20:28 2022

@author: veenstra
"""

import matplotlib.pyplot as plt
plt.close('all')
import xarray as xr
import cartopy.crs as ccrs
import numpy as np
import cartopy.feature as cf
 
file_list_nc = ['p:\\11206304-futuremares\\data\\CMIP6_BC\\CMCC-ESM2\\uo_Omon_CMCC-ESM2_ssp126_r1i1p1f1_gn_201501-203412.nc',
                'p:\\11206304-futuremares\\data\\CMIP6_BC\\CMCC-ESM2\\vo_Omon_CMCC-ESM2_ssp126_r1i1p1f1_gn_201501-203412.nc',
                ]
file_list_nc = ['p:\\11206304-futuremares\\data\\CMIP6_BC\\CMCC-ESM2\\uo_Omon_CMCC-ESM2_historical_r1i1p1f1_gn_201001-201412.nc',
                'p:\\11206304-futuremares\\data\\CMIP6_BC\\CMCC-ESM2\\vo_Omon_CMCC-ESM2_historical_r1i1p1f1_gn_201001-201412.nc',
                ]

def change_axes(ax1,ax2):
    for ax in (ax1,ax2):
        ax.set_xlim(-10,30)
        ax.set_ylim(30,70)
        ax.coastlines() #TODO: add this to plot_background instead of more complex option
        ax.gridlines(crs=crs, draw_labels=True, linewidth=.6, color='gray', alpha=0.5, linestyle='-.')
        
crs = ccrs.PlateCarree()
projection = ccrs.PlateCarree()
vmin = -1
vmax = 1

#separate lat/lon
data_xr_uo = xr.open_mfdataset(file_list_nc[0],chunks={'time':1})
data_xr_vo = xr.open_mfdataset(file_list_nc[1],chunks={'time':1})
fig, (ax1,ax2) = plt.subplots(1,2,figsize=(18,6),sharex=True,sharey=True,subplot_kw=dict(projection=projection))
data_xr_uo["uo"].isel(time=0,lev=0).plot.pcolormesh(ax=ax1, transform=crs, x='longitude', y='latitude', vmin=vmin, vmax=vmax, cmap='jet')
data_xr_vo["vo"].isel(time=0,lev=0).plot.pcolormesh(ax=ax2, transform=crs, x='longitude', y='latitude', vmin=vmin, vmax=vmax, cmap='jet')
fig.tight_layout()
change_axes(ax1,ax2)

#using lat/lon of uo for vo
data_xr = xr.open_mfdataset(file_list_nc,chunks={'time':1},compat='override')
fig, (ax1,ax2) = plt.subplots(1,2,figsize=(18,6),sharex=True,sharey=True,subplot_kw=dict(projection=projection))
data_xr["uo"].isel(time=0,lev=0).plot.pcolormesh(ax=ax1, transform=crs, x='longitude', y='latitude', vmin=vmin, vmax=vmax, cmap='jet')
data_xr["vo"].isel(time=0,lev=0).plot.pcolormesh(ax=ax2, transform=crs, x='longitude', y='latitude', vmin=vmin, vmax=vmax, cmap='jet')
fig.tight_layout()
change_axes(ax1,ax2)


print('i', (data_xr_uo.i.to_numpy() == data_xr_vo.i.to_numpy()).all())
print('j', (data_xr_uo.j.to_numpy() == data_xr_vo.j.to_numpy()).all())
print('latitude', (data_xr_uo.latitude.to_numpy() == data_xr_vo.latitude.to_numpy()).all())
print('latitude_diff', np.max(np.abs(data_xr_uo.latitude.to_numpy() - data_xr_vo.latitude.to_numpy())))
print('longitude', (data_xr_uo.longitude.to_numpy() == data_xr_vo.longitude.to_numpy()).all())
print('longitude_maxdiff', np.max(np.abs(data_xr_uo.longitude.to_numpy() - data_xr_vo.longitude.to_numpy())))

