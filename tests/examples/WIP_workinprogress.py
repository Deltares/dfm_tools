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

from dfm_tools.get_nc import get_netdata, get_ncmodeldata, plot_netmapdata
from dfm_tools.get_nc_helpers import get_ncvarproperties#, get_hisstationlist#, get_varname_fromnc
from dfm_tools.hydrolib_helpers import pointlike_to_DataFrame
from dfm_tools.xarray_helpers import preprocess_hirlam

from hydrolib.core.io.polyfile.models import PolyFile

dir_testinput = r'c:\DATA\dfm_tools_testdata'
dir_output = '.'

#print gridinfo of several files to compare
#file_nc = r'p:\1204257-dcsmzuno\2014\data\meteo\HIRLAM72_2018\h72_201803.nc' #TODO: xarray MissingDimensionsError
#print('\nfile = %s'%(file_nc))
#data_dummy = get_ncmodeldata(file_nc=file_nc, varname='northward_wind', timestep=0, get_linkedgridinfo=True)
file_nc = r'p:\archivedprojects\1220688-lake-kivu\2_data\COSMO\COSMOCLM_2012_out02_merged_4Wouter.nc'
print('\nfile = %s'%(file_nc))
data_dummy = get_ncmodeldata(file_nc=file_nc, varname='U_10M', timestep=0, get_linkedgridinfo=True)
file_nc = r'p:\11200665-c3s-codec\2_Hydro\ECWMF_meteo\meteo\ERA-5\2000\ERA5_metOcean_atm_19991201_19991231.nc'
print('\nfile = %s'%(file_nc))
data_dummy = get_ncmodeldata(file_nc=file_nc, varname='msl', timestep=0, get_linkedgridinfo=True)
file_nc = r'p:\11202255-sfincs\Testbed\Original_tests\01_Implementation\08_restartfile\sfincs_map.nc'
print('\nfile = %s'%(file_nc))
data_dummy = get_ncmodeldata(file_nc=file_nc, varname='zs', timestep=0, get_linkedgridinfo=True)
print('\nfile = %s'%(file_nc))
data_dummy = get_ncmodeldata(file_nc=file_nc, varname='u', timestep=0, get_linkedgridinfo=True)
file_nc = r'p:\archivedprojects\1220688-lake-kivu\3_modelling\1_FLOW\7_heatfluxinhis\063_netcdf\trim-thiery_002_coarse.nc'
print('\nfile = %s'%(file_nc))
data_dummy = get_ncmodeldata(file_nc=file_nc, varname='S1', timestep=0, get_linkedgridinfo=True)
print('\nfile = %s'%(file_nc))
data_dummy = get_ncmodeldata(file_nc=file_nc, varname='U1', timestep=0, layer=0, get_linkedgridinfo=True)
print('\nfile = %s'%(file_nc))
data_dummy = get_ncmodeldata(file_nc=file_nc, varname='V1', timestep=0, layer=0, get_linkedgridinfo=True)
file_nc = r'p:\1204257-dcsmzuno\2019\DCSMv6\A01\SDS-A01_map.nc'
print('\nfile = %s'%(file_nc))
data_dummy = get_ncmodeldata(file_nc=file_nc, varname='SEP', timestep=0, get_linkedgridinfo=True)
print('\nfile = %s'%(file_nc))
data_dummy = get_ncmodeldata(file_nc=file_nc, varname='VELU', timestep=0, layer=0, get_linkedgridinfo=True)
file_nc = r'p:\11203869-morwaqeco3d\05-Tidal_inlet\02_FM_201910\FM_MF10_Max_30s\wave\wavm-inlet.nc'
print('\nfile = %s'%(file_nc))
data_dummy = get_ncmodeldata(file_nc=file_nc, varname='veloc-x', timestep=0, get_linkedgridinfo=True)
file_nc = os.path.join(dir_testinput,r'DFM_3D_z_Grevelingen\computations\run01\DFM_OUTPUT_Grevelingen-FM\Grevelingen-FM_0000_map.nc')
print('\nfile = %s'%(file_nc))
data_dummy = get_ncmodeldata(file_nc=file_nc, varname='mesh2d_s1', timestep=0, multipart=False, get_linkedgridinfo=True)
data_dummy = get_ncmodeldata(file_nc=file_nc, varname='mesh2d_u1', timestep=0, layer=0, multipart=False, get_linkedgridinfo=True)
data_dummy = get_ncmodeldata(file_nc=file_nc, varname='mesh2d_flowelem_bl', multipart=False, get_linkedgridinfo=True)


# test Grevelingen (integrated example, where all below should move towards)
file_nc = os.path.join(dir_testinput,r'DFM_3D_z_Grevelingen\computations\run01\DFM_OUTPUT_Grevelingen-FM\Grevelingen-FM_0000_map.nc')
vars_pd = get_ncvarproperties(file_nc=file_nc)
ugrid = get_netdata(file_nc=file_nc)
fig, ax = plt.subplots()
plot_netmapdata(ugrid.verts, values=None, ax=None, linewidth=0.5, color="crimson", facecolor="None")
ax.set_aspect('equal')



#HIRLAM
file_nc = r'p:\1204257-dcsmzuno\2014\data\meteo\HIRLAM72_2018\h72_201803.nc'
data_xr = xr.open_mfdataset(file_nc,drop_variables=['x','y'],preprocess=preprocess_hirlam)
timestep = 0
coarsefac = 2 #coarsen dataset for more performance, but not necessary

data_u = data_xr['eastward_wind'].isel(time=timestep)
data_v = data_xr['northward_wind'].isel(time=timestep)
magn = np.sqrt(data_u**2 + data_v**2)

fig, ax = plt.subplots()
ax.plot(magn.longitude,magn.latitude,'-b',linewidth=0.2)
ax.plot(magn.longitude.T,magn.latitude.T,'-b',linewidth=0.2)
plt.savefig(os.path.join(dir_output,'hirlam_mesh'))

fig, ax = plt.subplots()
ax.pcolormesh(magn.longitude,magn.latitude,magn)
plt.savefig(os.path.join(dir_output,'hirlam_magn_pcolor'))



#plt.close('all')
from dfm_tools.regulargrid import center2corner
#COSMO
file_nc = r'p:\archivedprojects\1220688-lake-kivu\2_data\COSMO\COSMOCLM_2012_out02_merged_4Wouter.nc'
vars_pd = get_ncvarproperties(file_nc=file_nc)

xcen = get_ncmodeldata(file_nc=file_nc, varname='lon')
ycen = get_ncmodeldata(file_nc=file_nc, varname='lat')
xcor = center2corner(xcen)
ycor = center2corner(ycen)
data_U10M = get_ncmodeldata(file_nc=file_nc, varname='U_10M', timestep=range(20), get_linkedgridinfo=True)
data_V10M = get_ncmodeldata(file_nc=file_nc, varname='V_10M', timestep=range(20), get_linkedgridinfo=True)
#xcen, ycen = np.meshgrid(data_lon, data_lat)
magn = np.sqrt(data_U10M**2 + data_V10M**2)

fig, ax = plt.subplots()
ax.plot(xcen, ycen, '-b', linewidth=0.2)
ax.plot(xcen.T, ycen.T, '-b', linewidth=0.2)
plt.savefig(os.path.join(dir_output,'COSMO_mesh'))

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
    ax.set_title('t=%d (%s)'%(timestep, data_V10M.var_times.loc[timestep]))
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
    ax.set_title('t=%d (%s)'%(timestep, data_V10M.var_times.loc[timestep]))
    ax.set_aspect('equal')
    ax.plot(data_ldb['x'], data_ldb['y'], 'k', linewidth=0.5)
    thinning = 2
    ax.quiver(xcen[::thinning,::thinning], ycen[::thinning,::thinning], data_U10M[timestep,::thinning,::thinning], data_V10M[timestep,::thinning,::thinning], 
              color='w',scale=50,width=0.008)#, edgecolor='face', cmap='jet')
fig.tight_layout()
plt.savefig(os.path.join(dir_output,'COSMO_magn_curvedquiver'))



#ERA5
file_nc = r'p:\11200665-c3s-codec\2_Hydro\ECWMF_meteo\meteo\ERA-5\2000\ERA5_metOcean_atm_19991201_19991231.nc'
vars_pd = get_ncvarproperties(file_nc=file_nc)
data_lon = get_ncmodeldata(file_nc=file_nc, varname='longitude')
data_lat = get_ncmodeldata(file_nc=file_nc, varname='latitude')
data_psl = get_ncmodeldata(file_nc=file_nc, varname='msl',timestep=10, get_linkedgridinfo=True)

lons,lats = np.meshgrid(data_lon,data_lat)
fig, ax = plt.subplots()
ax.plot(lons, lats,'-b',linewidth=0.2)
ax.plot(lons.T, lats.T,'-b',linewidth=0.2)
plt.savefig(os.path.join(dir_output,'ERA5_mesh'))

fig, ax = plt.subplots()
#ax.pcolor(lons, lats, data_psl[0,:,:])
ax.pcolor(data_lon, data_lat, data_psl[0,:,:])
#plt.pcolor(mesh2d_node_x,mesh2d_node_y,airp,linewidth=0.5)
plt.savefig(os.path.join(dir_output,'ERA5_msl_pcolor'))






#SFINCS
file_nc = r'p:\11202255-sfincs\Testbed\Original_tests\01_Implementation\08_restartfile\sfincs_map.nc'
#file_nc = r'p:\11202255-sfincs\Testbed\Original_tests\03_Application\22_Tsunami_Japan_Sendai\sfincs_map.nc'
vars_pd = get_ncvarproperties(file_nc=file_nc)

data_fromnc_x = get_ncmodeldata(file_nc=file_nc, varname='x')
data_fromnc_y = get_ncmodeldata(file_nc=file_nc, varname='y')
data_fromnc_zs = get_ncmodeldata(file_nc=file_nc, varname='zs', timestep='all', get_linkedgridinfo=True)

fig, ax = plt.subplots()
ax.plot(data_fromnc_x, data_fromnc_y,'-b',linewidth=0.2)
ax.plot(data_fromnc_x.T, data_fromnc_y.T,'-b',linewidth=0.2)
plt.savefig(os.path.join(dir_output,'SFINCS_mesh'))    

fig, axs = plt.subplots(3,1, figsize=(14,9))
for iT, timestep in enumerate([0,1,10]):
    ax=axs[iT]
    pc = ax.pcolor(data_fromnc_x, data_fromnc_y, data_fromnc_zs[timestep,:,:],cmap='jet')
    pc.set_clim([0,0.15])
    cbar = fig.colorbar(pc, ax=ax)
    cbar.set_label('%s (%s)'%(data_fromnc_zs.var_varname, data_fromnc_zs.var_ncattrs['units']))
    ax.set_title('t=%d (%s)'%(timestep, data_fromnc_zs.var_times.loc[timestep]))
    ax.set_aspect('equal')
fig.tight_layout()
plt.savefig(os.path.join(dir_output,'SFINCS_zs_pcolor'))


data_fromnc_edgex = get_ncmodeldata(file_nc=file_nc, varname='edge_x')
data_fromnc_edgey = get_ncmodeldata(file_nc=file_nc, varname='edge_y')
data_fromnc_u = get_ncmodeldata(file_nc=file_nc, varname='u', timestep='all')
data_fromnc_v = get_ncmodeldata(file_nc=file_nc, varname='v', timestep='all')    
vel_magn = np.sqrt(data_fromnc_u**2 + data_fromnc_v**2)

fig, ax = plt.subplots()
ax.plot(data_fromnc_edgex, data_fromnc_edgey,'-b',linewidth=0.2)
ax.plot(data_fromnc_edgex.T, data_fromnc_edgey.T,'-b',linewidth=0.2)
plt.savefig(os.path.join(dir_output,'SFINCS_meshedge'))    

fig, axs = plt.subplots(3,1, figsize=(14,9))
for iT, timestep in enumerate([0,1,10]):
    ax=axs[iT]
    pc = ax.pcolor(data_fromnc_edgex, data_fromnc_edgey,vel_magn[timestep,:,:],cmap='jet')
    pc.set_clim([0,0.6])
    cbar = fig.colorbar(pc, ax=ax)
    cbar.set_label('velocity magnitude (%s)'%(data_fromnc_u.var_ncattrs['units']))
    ax.set_title('t=%d (%s)'%(timestep, data_fromnc_u.var_times.loc[timestep]))
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
#data_fromnc_his = get_ncmodeldata(file_nc=file_nc, varname='point_zs', station='all', timestep='all')

fig, ax = plt.subplots()
for iS,stat_name in enumerate(stations_pd):
    data_sel = data_xr.point_zs.isel(stations=iS)
    ax.plot(data_sel.time, data_sel, label=stat_name)
ax.legend()
plt.savefig(os.path.join(dir_output,'SFINCS_hiszs'))

