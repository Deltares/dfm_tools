# -*- coding: utf-8 -*-
"""
Created on Wed Oct 27 22:03:01 2021

@author: veenstra

"""

import os
import numpy as np
import matplotlib.pyplot as plt
plt.close('all')
from pathlib import Path
import xarray as xr

#from dfm_tools.get_nc import get_ncmodeldata
#from dfm_tools.get_nc_helpers import get_ncvarproperties#, get_hisstationlist#, get_varname_fromnc
from dfm_tools.hydrolib_helpers import pointlike_to_DataFrame
from hydrolib.core.io.polyfile.models import PolyFile

dir_testinput = r'c:\DATA\dfm_tools_testdata'
dir_output = '.'


#COSMO
file_nc = r'p:\archivedprojects\1220688-lake-kivu\2_data\COSMO\COSMOCLM_2012_out02_merged_4Wouter.nc'
#vars_pd = get_ncvarproperties(file_nc=file_nc)
data_xr = xr.open_dataset(file_nc)
data_xr = data_xr.drop(['height_10m','height_2m']) #gives cleaner figure title
data_U10M = data_xr['U_10M']
data_V10M = data_xr['V_10M']

file_ldb = Path(r'p:\archivedprojects\1220688-lake-kivu\3_modelling\1_FLOW\4_CH4_CO2_included\008\lake_kivu_geo.ldb')
polyfile_object = PolyFile(file_ldb)
data_ldb = pointlike_to_DataFrame(polyfile_object.objects[0])
data_ldb[data_ldb==999.999] = np.nan


fig, axs = plt.subplots(1,3, figsize=(16,6))
for iT, timestep in enumerate([0,1,10]):
    ax=axs[iT]
    data_U10M_tsel = data_U10M.isel(time=timestep)
    data_V10M_tsel = data_V10M.isel(time=timestep)
    data_magn = np.sqrt(data_U10M_tsel**2 + data_V10M_tsel**2)
    pc = data_magn.plot.pcolormesh(x='lon',y='lat', cmap='jet',ax=ax)
    pc.set_clim([0,5])
    ax.set_aspect('equal')
    ax.plot(data_ldb['x'], data_ldb['y'], 'k', linewidth=0.5)
    thinning = 2
    data_U10M_thin = data_U10M_tsel.loc[::thinning,::thinning]
    data_V10M_thin = data_V10M_tsel.loc[::thinning,::thinning]
    ax.quiver(data_U10M_thin.lon, data_U10M_thin.lat, data_U10M_thin, data_V10M_thin, 
              color='w',scale=50,width=0.008)#, edgecolor='face', cmap='jet')
fig.tight_layout()
plt.savefig(os.path.join(dir_output,'COSMO_magn_pcolorquiver'))


dist = 0.106
reg_x_vec = np.arange(data_xr.lon.min(),data_xr.lon.max()+dist,dist)
reg_y_vec = np.arange(data_xr.lat.min(),data_xr.lat.max()+dist,dist)
reg_grid_X,reg_grid_Y = np.meshgrid(reg_x_vec,reg_y_vec)

from scipy.interpolate import griddata
from dfm_tools.modplot import velovect

fig, axs = plt.subplots(1,3, figsize=(16,6),sharex=True,sharey=True)
for iT, timestep in enumerate([0,1,10]):
    ax=axs[iT]
    data_U10M_tsel = data_U10M.isel(time=timestep)
    data_V10M_tsel = data_V10M.isel(time=timestep)
    lonvals_flat = data_U10M_tsel.lon.to_numpy().flatten()
    latvals_flat = data_U10M_tsel.lat.to_numpy().flatten()
    U = griddata((lonvals_flat,latvals_flat),data_U10M_tsel.to_numpy().flatten(),(reg_grid_X,reg_grid_Y),method='nearest') #TODO: this is probably easier with xarray
    V = griddata((lonvals_flat,latvals_flat),data_V10M_tsel.to_numpy().flatten(),(reg_grid_X,reg_grid_Y),method='nearest')
    speed = np.sqrt(U*U + V*V)
    quiv_curved = velovect(ax,reg_grid_X,reg_grid_Y,U,V, arrowstyle='fancy', scale = 5, grains = 25, color=speed)#, cmap='jet')
    quiv_curved.lines.set_clim(0,5)
    #quiv_curved = ax.streamplot(reg_x_vec,reg_y_vec,U,V, arrowstyle='-|>', integration_direction='forward',broken_streamlines=False, color=speed, density=1.5)
    cbar = fig.colorbar(quiv_curved.lines, ax=ax)
    cbar.set_label('velocity magnitude (%s)'%(data_V10M.attrs['units']))
    ax.set_title(data_V10M_tsel.time.to_pandas())
    #ax.set_aspect('equal')
    ax.plot(data_ldb['x'], data_ldb['y'], 'k', linewidth=0.5)
fig.tight_layout()
plt.savefig(os.path.join(dir_output,'COSMO_magn_curvedquiver'))



#SFINCS
file_nc = r'p:\11202255-sfincs\Testbed\Original_tests\01_Implementation\08_restartfile\sfincs_map.nc'
#file_nc = r'p:\11202255-sfincs\Testbed\Original_tests\03_Application\22_Tsunami_Japan_Sendai\sfincs_map.nc'
#vars_pd = get_ncvarproperties(file_nc=file_nc)

data_xr = xr.open_dataset(file_nc)
data_xr = data_xr.set_coords(['x','y','edge_x','edge_y'])
data_fromnc_zs = data_xr['zs']

fig, axs = plt.subplots(3,1, figsize=(14,9))
for iT, timestep in enumerate([0,1,10]):
    ax=axs[iT]
    pc = data_fromnc_zs.isel(time=timestep).plot.pcolormesh(x='x',y='y',cmap='jet', ax=ax)
    pc.set_clim([0,0.15])
    ax.set_aspect('equal')
fig.tight_layout()
plt.savefig(os.path.join(dir_output,'SFINCS_zs_pcolor'))


data_u = data_xr['u']
data_v = data_xr['v']

fig, axs = plt.subplots(3,1, figsize=(14,9))
for iT, timestep in enumerate([0,1,10]):
    ax=axs[iT]
    data_u_tsel = data_u.isel(time=timestep)
    data_v_tsel = data_v.isel(time=timestep)
    vel_magn = np.sqrt(data_u_tsel**2 + data_v_tsel**2) 
    pc = vel_magn.plot.pcolormesh(x='edge_x',y='edge_y',cmap='jet',ax=ax)
    pc.set_clim([0,0.6])
    ax.set_aspect('equal')
    thinning = 5
    data_u_thin = data_u_tsel.loc[::thinning,::thinning]
    data_v_thin = data_v_tsel.loc[::thinning,::thinning]
    ax.quiver(data_u_thin.edge_x, data_u_thin.edge_y, data_u_thin, data_v_thin, 
              color='w')#,scale=3,width=0.005)#, edgecolor='face', cmap='jet')
fig.tight_layout()
plt.savefig(os.path.join(dir_output,'SFINCS_velocity_pcolorquiver'))


#SFINCS HIS
#file_nc = r'p:\11202255-sfincs\Testbed\Original_tests\01_Implementation\14_restartfile\sfincs_his.nc'
file_nc = r'p:\11202255-sfincs\Testbed\Original_tests\03_Application\04_Tsunami_Japan_Sendai\sfincs_his.nc'
#vars_pd = get_ncvarproperties(file_nc=file_nc)
data_xr = xr.open_dataset(file_nc)
#station_names = get_hisstationlist(file_nc=file_nc, varname='point_zs')
stations_pd = data_xr.station_name.astype(str).to_pandas()

fig, ax = plt.subplots()
for iS,stat_name in enumerate(stations_pd):
    data_sel = data_xr.point_zs.isel(stations=iS)
    ax.plot(data_sel.time, data_sel, label=stat_name)
ax.legend()
plt.savefig(os.path.join(dir_output,'SFINCS_hiszs'))

