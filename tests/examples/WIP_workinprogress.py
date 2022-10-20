# -*- coding: utf-8 -*-
"""
Created on Wed Oct 27 22:03:01 2021

@author: veenstra


WARNING: THIS TEST IS NOT YET FINISHED, WILL BE IMPROVED AND LINKED TO INTERNAL FUNCTIONS ASAP
"""

import os
import numpy as np
import matplotlib.pyplot as plt
plt.close('all')
from pathlib import Path
import xarray as xr

from dfm_tools.get_nc import get_ncmodeldata
from dfm_tools.get_nc_helpers import get_ncvarproperties#, get_hisstationlist#, get_varname_fromnc
from dfm_tools.hydrolib_helpers import pointlike_to_DataFrame

from hydrolib.core.io.polyfile.models import PolyFile

dir_testinput = r'c:\DATA\dfm_tools_testdata'
dir_output = '.'


#TODO: move to xarray and see if functions like center2corner are still used

#plt.close('all')
from dfm_tools.regulargrid import center2corner
#COSMO
file_nc = r'p:\archivedprojects\1220688-lake-kivu\2_data\COSMO\COSMOCLM_2012_out02_merged_4Wouter.nc'
vars_pd = get_ncvarproperties(file_nc=file_nc)

xcen = get_ncmodeldata(file_nc=file_nc, varname='lon')
ycen = get_ncmodeldata(file_nc=file_nc, varname='lat')
xcor = center2corner(xcen)
ycor = center2corner(ycen)
data_U10M = get_ncmodeldata(file_nc=file_nc, varname='U_10M', timestep=range(20))
data_V10M = get_ncmodeldata(file_nc=file_nc, varname='V_10M', timestep=range(20))
#xcen, ycen = np.meshgrid(data_lon, data_lat)
magn = np.sqrt(data_U10M**2 + data_V10M**2)


file_ldb = Path(r'p:\archivedprojects\1220688-lake-kivu\3_modelling\1_FLOW\4_CH4_CO2_included\008\lake_kivu_geo.ldb')
polyfile_object = PolyFile(file_ldb)
data_ldb = pointlike_to_DataFrame(polyfile_object.objects[0])
data_ldb[data_ldb==999.999] = np.nan

fig, axs = plt.subplots(1,3, figsize=(16,6))
for iT, timestep in enumerate([0,1,10]):
    ax=axs[iT]
    #pc = ax.pcolor(xcen, ycen, magn[timestep,:,:], cmap='jet')
    pc = ax.pcolor(xcor, ycor, magn[timestep,:,:], cmap='jet')
    pc.set_clim([0,5])
    cbar = fig.colorbar(pc, ax=ax)
    cbar.set_label('velocity magnitude (%s)'%(data_V10M.var_ncattrs['units']))
    ax.set_title('t=%d (%s)'%(timestep, data_V10M.var_times.iloc[timestep]))
    ax.set_aspect('equal')
    ax.plot(data_ldb['x'], data_ldb['y'], 'k', linewidth=0.5)
    thinning = 2
    ax.quiver(xcen[::thinning,::thinning], ycen[::thinning,::thinning], data_U10M[timestep,::thinning,::thinning], data_V10M[timestep,::thinning,::thinning], 
              color='w',scale=50,width=0.008)#, edgecolor='face', cmap='jet')
fig.tight_layout()
plt.savefig(os.path.join(dir_output,'COSMO_magn_pcolorquiver'))



dist = 0.1
reg_x_vec = np.linspace(np.min(xcen),np.max(xcen),int(np.ceil((np.max(xcen)-np.min(xcen))/dist)))
reg_y_vec = np.linspace(np.min(ycen),np.max(ycen),int(np.ceil((np.max(ycen)-np.min(ycen))/dist)))
reg_grid = np.meshgrid(reg_x_vec,reg_y_vec)
X = reg_grid[0]
Y = reg_grid[1]
from scipy.interpolate import griddata
from dfm_tools.modplot import velovect

fig, axs = plt.subplots(1,3, figsize=(16,6))
for iT, timestep in enumerate([0,1,10]):
    ax=axs[iT]
    #pc = ax.pcolor(xcen, ycen, magn[timestep,:,:], cmap='jet')
    #pc.set_clim([0,5])
    U = griddata((xcen.flatten(),ycen.flatten()),data_U10M[timestep,:,:].flatten(),tuple(reg_grid),method='nearest')
    V = griddata((xcen.flatten(),ycen.flatten()),data_V10M[timestep,:,:].flatten(),tuple(reg_grid),method='nearest')
    speed = np.sqrt(U*U + V*V)
    quiv_curved = velovect(ax,X,Y,U,V, arrowstyle='fancy', scale = 5, grains = 25, color=speed)#, cmap='jet')
    ax.set_aspect('equal')
    cbar = fig.colorbar(quiv_curved.lines, ax=ax)
    cbar.set_label('velocity magnitude (%s)'%(data_V10M.var_ncattrs['units']))
    ax.set_title('t=%d (%s)'%(timestep, data_V10M.var_times.iloc[timestep]))
    ax.set_aspect('equal')
    ax.plot(data_ldb['x'], data_ldb['y'], 'k', linewidth=0.5)
    thinning = 2
    ax.quiver(xcen[::thinning,::thinning], ycen[::thinning,::thinning], data_U10M[timestep,::thinning,::thinning], data_V10M[timestep,::thinning,::thinning], 
              color='w',scale=50,width=0.008)#, edgecolor='face', cmap='jet')
fig.tight_layout()
plt.savefig(os.path.join(dir_output,'COSMO_magn_curvedquiver'))





#SFINCS
file_nc = r'p:\11202255-sfincs\Testbed\Original_tests\01_Implementation\08_restartfile\sfincs_map.nc'
#file_nc = r'p:\11202255-sfincs\Testbed\Original_tests\03_Application\22_Tsunami_Japan_Sendai\sfincs_map.nc'
vars_pd = get_ncvarproperties(file_nc=file_nc)

data_fromnc_x = get_ncmodeldata(file_nc=file_nc, varname='x')
data_fromnc_y = get_ncmodeldata(file_nc=file_nc, varname='y')
data_fromnc_zs = get_ncmodeldata(file_nc=file_nc, varname='zs', timestep='all')

fig, axs = plt.subplots(3,1, figsize=(14,9))
for iT, timestep in enumerate([0,1,10]):
    ax=axs[iT]
    pc = ax.pcolor(data_fromnc_x, data_fromnc_y, data_fromnc_zs[timestep,:,:],cmap='jet')
    pc.set_clim([0,0.15])
    cbar = fig.colorbar(pc, ax=ax)
    cbar.set_label('%s (%s)'%(data_fromnc_zs.var_varname, data_fromnc_zs.var_ncattrs['units']))
    ax.set_title('t=%d (%s)'%(timestep, data_fromnc_zs.var_times.iloc[timestep]))
    ax.set_aspect('equal')
fig.tight_layout()
plt.savefig(os.path.join(dir_output,'SFINCS_zs_pcolor'))


data_fromnc_edgex = get_ncmodeldata(file_nc=file_nc, varname='edge_x')
data_fromnc_edgey = get_ncmodeldata(file_nc=file_nc, varname='edge_y')
data_fromnc_u = get_ncmodeldata(file_nc=file_nc, varname='u', timestep='all')
data_fromnc_v = get_ncmodeldata(file_nc=file_nc, varname='v', timestep='all')    
vel_magn = np.sqrt(data_fromnc_u**2 + data_fromnc_v**2) 

fig, axs = plt.subplots(3,1, figsize=(14,9))
for iT, timestep in enumerate([0,1,10]):
    ax=axs[iT]
    pc = ax.pcolor(data_fromnc_edgex, data_fromnc_edgey,vel_magn[timestep,:,:],cmap='jet')
    pc.set_clim([0,0.6])
    cbar = fig.colorbar(pc, ax=ax)
    cbar.set_label('velocity magnitude (%s)'%(data_fromnc_u.var_ncattrs['units']))
    ax.set_title('t=%d (%s)'%(timestep, data_fromnc_u.var_times.iloc[timestep]))
    ax.set_aspect('equal')
    thinning = 5
    ax.quiver(data_fromnc_edgex[::thinning,::thinning], data_fromnc_edgey[::thinning,::thinning], data_fromnc_u[timestep,::thinning,::thinning], data_fromnc_v[timestep,::thinning,::thinning], 
              color='w')#,scale=3,width=0.005)#, edgecolor='face', cmap='jet')
fig.tight_layout()
plt.savefig(os.path.join(dir_output,'SFINCS_velocity_pcolorquiver'))


#SFINCS HIS
#file_nc = r'p:\11202255-sfincs\Testbed\Original_tests\01_Implementation\14_restartfile\sfincs_his.nc'
file_nc = r'p:\11202255-sfincs\Testbed\Original_tests\03_Application\04_Tsunami_Japan_Sendai\sfincs_his.nc'
vars_pd = get_ncvarproperties(file_nc=file_nc)
data_xr = xr.open_dataset(file_nc)
#station_names = get_hisstationlist(file_nc=file_nc, varname='point_zs')
stations_pd = data_xr.station_name.astype(str).to_pandas()

fig, ax = plt.subplots()
for iS,stat_name in enumerate(stations_pd):
    data_sel = data_xr.point_zs.isel(stations=iS)
    ax.plot(data_sel.time, data_sel, label=stat_name)
ax.legend()
plt.savefig(os.path.join(dir_output,'SFINCS_hiszs'))

