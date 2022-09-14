# -*- coding: utf-8 -*-
"""
Created on Wed Oct 27 22:01:16 2021

@author: veenstra
"""

import os
import numpy as np
import matplotlib.pyplot as plt
plt.close('all')
import xarray as xr

from dfm_tools.get_nc import get_netdata, get_ncmodeldata, plot_netmapdata#, get_xzcoords_onintersection
from dfm_tools.get_nc_helpers import get_ncvardimlist, get_hisstationlist
from dfm_tools.regulargrid import scatter_to_regulargrid#, meshgridxy2verts, center2corner

dir_testinput = r'c:\DATA\dfm_tools_testdata'
dir_output = '.'


#MAPFILE
file_nc = r'p:\11203869-morwaqeco3d\05-Tidal_inlet\02_FM_201910\FM_MF10_Max_30s\fm\DFM_OUTPUT_inlet\inlet_map.nc'
vars_pd, dims_pd = get_ncvardimlist(file_nc=file_nc)
vars_pd.to_csv(os.path.join(dir_output,'vars_pd.csv'))
vars_pd_sel = vars_pd[vars_pd['long_name'].str.contains('transport')]
#vars_pd_sel = vars_pd[vars_pd['dimensions'].str.contains('mesh2d_nFaces') & vars_pd['long_name'].str.contains('wave')]

ugrid = get_netdata(file_nc=file_nc)

varname = 'mesh2d_mor_bl'
var_clims = [-50,0]
var_longname = vars_pd['long_name'][vars_pd['nc_varkeys']==varname].iloc[0]
fig, axs = plt.subplots(3,1, figsize=(6,9))
fig.suptitle('%s (%s)'%(varname, var_longname))

ax = axs[0]
data_frommap_0 = get_ncmodeldata(file_nc=file_nc, varname=varname, timestep=0, get_linkedgridinfo=True)
pc = plot_netmapdata(ugrid.verts, values=data_frommap_0.flatten(), ax=ax, linewidth=0.5, cmap='jet', clim=var_clims)
cbar = fig.colorbar(pc, ax=ax)
cbar.set_label('%s (%s)'%(data_frommap_0.var_varname, data_frommap_0.var_ncattrs['units']))
ax.set_title('t=0 (%s)'%(data_frommap_0.var_times.iloc[0]))

ax = axs[1]
data_frommap_end = get_ncmodeldata(file_nc=file_nc, varname=varname, timestep=-1)
pc = plot_netmapdata(ugrid.verts, values=data_frommap_end.flatten(), ax=ax, linewidth=0.5, cmap='jet', clim=var_clims)
cbar = fig.colorbar(pc, ax=ax)
cbar.set_label('%s (%s)'%(data_frommap_end.var_varname, data_frommap_end.var_ncattrs['units']))
ax.set_title('t=end (%s)'%(data_frommap_end.var_times.iloc[0]))

ax = axs[2]
pc = plot_netmapdata(ugrid.verts, values=(data_frommap_end-data_frommap_0).flatten(), ax=ax, linewidth=0.5, cmap='jet', clim=[-3,3])
cbar = fig.colorbar(pc, ax=ax)
cbar.set_label('%s (%s)'%(data_frommap_0.var_varname, data_frommap_0.var_ncattrs['units']))
ax.set_title('t=end-0 (difference)')

for ax in axs:
    ax.set_aspect('equal')
    #ax.set_ylim(val_ylim)
fig.tight_layout()
plt.savefig(os.path.join(dir_output,'%s_%s'%(os.path.basename(file_nc).replace('.',''), varname)))



varname = 'mesh2d_hwav'
var_longname = vars_pd['long_name'][vars_pd['nc_varkeys']==varname].iloc[0]
fig, ax = plt.subplots(1,1)
fig.suptitle('%s (%s)'%(varname, var_longname))

data_frommap = get_ncmodeldata(file_nc=file_nc, varname=varname, timestep=-1)
pc = plot_netmapdata(ugrid.verts, values=data_frommap.flatten(), ax=ax, linewidth=0.5, cmap='jet')
cbar = fig.colorbar(pc, ax=ax)
cbar.set_label('%s (%s)'%(data_frommap.var_varname, data_frommap.var_ncattrs['units']))
ax.set_title('t=end (%s)'%(data_frommap.var_times.iloc[0]))
ax.set_aspect('equal')

fig.tight_layout()
plt.savefig(os.path.join(dir_output,'%s_%s'%(os.path.basename(file_nc).replace('.',''), varname)))




#file_nc = r'p:\11203869-morwaqeco3d\05-Tidal_inlet\02_FM_201910\FM_MF10_Max_30s\fm\DFM_OUTPUT_inlet\inlet_com.nc'
"""
#COMFILE
vars_pd, dims_pd = get_ncvardimlist(file_nc=file_nc)
vars_pd_sel = vars_pd[vars_pd['long_name'].str.contains('wave')]
#vars_pd_sel = vars_pd[vars_pd['dimensions'].str.contains('mesh2d_nFaces') & vars_pd['long_name'].str.contains('wave')]

ugrid = get_netdata(file_nc=file_nc)

#construct different ugrid (with bnds?)
data_fromnc_FlowElemContour_x = get_ncmodeldata(file_nc=file_nc, varname='FlowElemContour_x')
data_fromnc_FlowElemContour_y = get_ncmodeldata(file_nc=file_nc, varname='FlowElemContour_y')
data_fromnc_FlowElemContour_xy = np.stack([data_fromnc_FlowElemContour_x,data_fromnc_FlowElemContour_y],axis=2)

varname_list = ['hrms', 'tp', 'dir']#, 'distot', 'wlen']
for varname in varname_list:
    var_longname = vars_pd['long_name'][vars_pd['nc_varkeys']==varname].iloc[0]
    fig, ax = plt.subplots()#fig, axs = plt.subplots(2,1, figsize=(6,8))
    fig.suptitle('%s (%s)'%(varname, var_longname))
    
    timestep = 0
    data_frommap = get_ncmodeldata(file_nc=file_nc, varname=varname, timestep=timestep)
    pc = plot_netmapdata(data_fromnc_FlowElemContour_xy, values=data_frommap.flatten(), ax=ax, linewidth=0.5, cmap='jet')
    cbar = fig.colorbar(pc, ax=ax)
    cbar.set_label('%s (%s)'%(data_frommap.var_varname, data_frommap.var_ncattrs['units']))
    ax.set_title('t=%d (%s)'%(timestep, data_frommap.var_times.iloc[0]))
    ax.set_aspect('equal')
    
    fig.tight_layout()
    plt.savefig(os.path.join(dir_output,'%s_%s'%(os.path.basename(file_nc).replace('.',''), varname)))
"""

#WAVM FILE
file_nc = r'p:\11203869-morwaqeco3d\05-Tidal_inlet\02_FM_201910\FM_MF10_Max_30s\wave\wavm-inlet.nc'
vars_pd, dims_pd = get_ncvardimlist(file_nc=file_nc)
vars_pd_sel = vars_pd[vars_pd['long_name'].str.contains('dissi')]
#vars_pd_sel = vars_pd[vars_pd['dimensions'].str.contains('mesh2d_nFaces') & vars_pd['long_name'].str.contains('wave')]


#get cell center coordinates from regular grid, convert to grid_verts on corners
data_fromnc_x = get_ncmodeldata(file_nc=file_nc, varname='x')
data_fromnc_y = get_ncmodeldata(file_nc=file_nc, varname='y')
#x_cen_withbnd = center2corner(data_fromnc_x)
#y_cen_withbnd = center2corner(data_fromnc_y)
#grid_verts = meshgridxy2verts(x_cen_withbnd, y_cen_withbnd)

#plt.close('all')
varname_list = ['hsign', 'dir', 'period', 'dspr', 'dissip']
var_clim = [[0,2], [0,360], [0,7.5], [0,35], [0,20]]
for iV, varname in enumerate(varname_list):
    var_longname = vars_pd['long_name'][vars_pd['nc_varkeys']==varname].iloc[0]
    
    fig, axs = plt.subplots(2,1, figsize=(12,7))
    fig.suptitle('%s (%s)'%(varname, var_longname))

    timestep = 10
    data_frommap = get_ncmodeldata(file_nc=file_nc, varname=varname, timestep=timestep, get_linkedgridinfo=True)
    ax = axs[0]
    pc = ax.pcolor(data_fromnc_x, data_fromnc_y, data_frommap[0,:,:], cmap='jet')
    pc.set_clim(var_clim[iV])
    cbar = fig.colorbar(pc, ax=ax)
    cbar.set_label('%s (%s)'%(data_frommap.var_varname, data_frommap.var_ncattrs['units']))
    ax.set_title('t=%d (%s)'%(timestep, data_frommap.var_times.iloc[0]))
    ax.set_aspect('equal')
    
    timestep = -1
    data_frommap = get_ncmodeldata(file_nc=file_nc, varname=varname, timestep=timestep)        
    ax = axs[1]
    pc = ax.pcolor(data_fromnc_x, data_fromnc_y, data_frommap[0,:,:], cmap='jet')
    pc.set_clim(var_clim[iV])
    cbar = fig.colorbar(pc, ax=ax)
    cbar.set_label('%s (%s)'%(data_frommap.var_varname, data_frommap.var_ncattrs['units']))
    ax.set_title('t=%d (%s)'%(timestep, data_frommap.var_times.iloc[0]))
    ax.set_aspect('equal')
    
    fig.tight_layout()
    plt.savefig(os.path.join(dir_output,'%s_%s'%(os.path.basename(file_nc).replace('.',''), varname)))
    
    if varname == 'dir':
        #also plot with vectors
        ax = axs[0]
        ax.clear()
        pc = ax.quiver(data_fromnc_x, data_fromnc_y, 1,1,data_frommap[0,:,:],
                       angles=90-data_frommap[0,:,:], cmap='jet', scale=100)
        for ax in axs:
            ax.set_title('t=%d (%s)'%(timestep, data_frommap.var_times.iloc[0]))
            ax.set_aspect('equal')
        plt.savefig(os.path.join(dir_output,'%s_%s_vec'%(os.path.basename(file_nc).replace('.',''), varname)))
        for ax in axs:
            ax.set_xlim([25000,65000])
            ax.set_ylim([2500,15000])
        plt.savefig(os.path.join(dir_output,'%s_%s_veczoom'%(os.path.basename(file_nc).replace('.',''), varname)))



#HISFILE
file_nc = r'p:\11203869-morwaqeco3d\05-Tidal_inlet\02_FM_201910\FM_MF10_Max_30s\fm\DFM_OUTPUT_inlet\inlet_his.nc'
vars_pd, dims_pd = get_ncvardimlist(file_nc=file_nc)
vars_pd_sel = vars_pd[vars_pd['long_name'].str.contains('level')]
stat_list = get_hisstationlist(file_nc,varname='station_name')
crs_list = get_hisstationlist(file_nc,varname='cross_section_name')
data_xr = xr.open_dataset(file_nc)
stations_pd = data_xr.station_name.astype(str).to_pandas()

var_names = ['waterlevel','bedlevel']#,'mesh2d_ssn']
for iV, varname in enumerate(var_names):
    #data_fromhis = get_ncmodeldata(file_nc=file_nc, varname=varname, timestep='all', station='all')

    fig, ax = plt.subplots(1,1, figsize=(10,5))
    for iS, stat in enumerate(stations_pd):
        data_fromhis = data_xr.get(varname).isel(stations=iS,time=slice(0,3000))
        var_longname = data_fromhis.attrs['long_name'] #vars_pd['long_name'][vars_pd['nc_varkeys']==varname].iloc[0]
        ax.plot(data_fromhis.time, data_fromhis, linewidth=1, label=stat)
    ax.legend()
    ax.set_ylabel('%s (%s)'%(data_fromhis.name,data_fromhis.attrs['units']))
    ax.set_xlim(data_fromhis.time[[0,-1]].to_numpy())
    fig.tight_layout()
    plt.savefig(os.path.join(dir_output,'%s_%s'%(os.path.basename(file_nc).replace('.',''), varname)))






#MAPFILE TRANSPORT
file_nc = r'p:\11203869-morwaqeco3d\05-Tidal_inlet\02_FM_201910\FM_MF10_Max_30s\fm\DFM_OUTPUT_inlet\inlet_map.nc'
#file_nc = r'p:\11203869-morwaqeco3d\04-Breakwater\02_FM_201910\01_FM_MF25_Max_30s_User_1200s\fm\DFM_OUTPUT_straight_coast\straight_coast_map.nc'
vars_pd, dims_pd = get_ncvardimlist(file_nc=file_nc)
data_xr = xr.open_dataset(file_nc)
#vars_pd_sel = vars_pd[vars_pd['long_name'].str.contains('transport')]
#vars_pd_sel = vars_pd[vars_pd['dimensions'].str.contains('mesh2d_nFaces') & vars_pd['long_name'].str.contains('wave')]

ugrid = get_netdata(file_nc=file_nc)
timestep = 10
data_frommap_facex = data_xr.mesh2d_face_x # get_ncmodeldata(file_nc=file_nc, varname='mesh2d_face_x')
data_frommap_facey = data_xr.mesh2d_face_y # get_ncmodeldata(file_nc=file_nc, varname='mesh2d_face_y')
data_frommap_transx = data_xr.mesh2d_sxtot.isel(time=[-1]) #get_ncmodeldata(file_nc=file_nc, varname='mesh2d_sxtot', timestep=timestep, station='all')
data_frommap_transy = data_xr.mesh2d_sytot.isel(time=[-1]) #get_ncmodeldata(file_nc=file_nc, varname='mesh2d_sytot', timestep=timestep, station='all')
magnitude = (data_frommap_transx ** 2 + data_frommap_transy ** 2) ** 0.5

#plt.close('all')
fig, ax = plt.subplots(1,1, figsize=(14,8))
quiv = ax.quiver(data_frommap_facex, data_frommap_facey, data_frommap_transx[0,0,:], data_frommap_transy[0,0,:],
                 magnitude[0,0,:])#, scale=0.015)
cbar = fig.colorbar(quiv, ax=ax)
cbar.set_label('%s and %s (%s)'%(data_frommap_transx.name, data_frommap_transy.name, data_frommap_transy.attrs['units']))
ax.set_title('t=%d (%s)'%(timestep, data_frommap_transx.time.to_pandas().iloc[0]))
ax.set_aspect('equal')
fig.tight_layout()
plt.savefig(os.path.join(dir_output,'%s_%s_%s_t%d'%(os.path.basename(file_nc).replace('.',''), data_frommap_transx.name, data_frommap_transy.name,timestep)))
xlim_get = ax.get_xlim()
ylim_get = ax.get_ylim()

#interpolate to regular grid
X,Y,U = scatter_to_regulargrid(xcoords=data_frommap_facex, ycoords=data_frommap_facey, ncellx=29, ncelly=20, values=data_frommap_transx[0,0,:])
X,Y,V = scatter_to_regulargrid(xcoords=data_frommap_facex, ycoords=data_frommap_facey, ncellx=29, ncelly=20, values=data_frommap_transy[0,0,:])
speed = np.sqrt(U*U + V*V)

fig, ax = plt.subplots(1,1, figsize=(14,8))
quiv = ax.quiver(X, Y, U, V, speed)
cbar = fig.colorbar(quiv, ax=ax)
cbar.set_label('%s and %s (%s)'%(data_frommap_transx.name, data_frommap_transy.name, data_frommap_transy.attrs['units']))
ax.set_title('t=%d (%s)'%(timestep, data_frommap_transx.time.to_pandas().iloc[0]))
ax.set_xlim(xlim_get)
ax.set_ylim(ylim_get)
ax.set_aspect('equal')
fig.tight_layout()
plt.savefig(os.path.join(dir_output,'%s_%s_%s_t%d_regquiver'%(os.path.basename(file_nc).replace('.',''), data_frommap_transx.name, data_frommap_transy.name,timestep)))

#xs = X.flatten()
#ys = Y.flatten()
#seed_points = np.array([list(xs), list(ys)])
fig, ax = plt.subplots(1,1, figsize=(14,8))
strm = ax.streamplot(X, Y, U, V, color=speed, density=2, linewidth=1+2*speed/np.max(speed))#, cmap='winter', 
#                      minlength=0.01, maxlength = 2, arrowstyle='fancy')#,
#                      integration_direction='forward')#, start_points = seed_points.T)
#strm = ax.streamplot(X, Y, U, V, color=speed, linewidth=1+2*speed/np.max(speed), density=10,# cmap='winter',
#                     minlength=0.0001, maxlength = 0.07, arrowstyle='fancy',
#                     integration_direction='forward', start_points = seed_points.T)
cbar = fig.colorbar(strm.lines)
cbar.set_label('%s and %s (%s)'%(data_frommap_transx.name, data_frommap_transy.name, data_frommap_transy.attrs['units']))
ax.set_title('t=%d (%s)'%(timestep, data_frommap_transx.time.to_pandas().iloc[0]))
ax.set_xlim(xlim_get)
ax.set_ylim(ylim_get)
ax.set_aspect('equal')
fig.tight_layout()
plt.savefig(os.path.join(dir_output,'%s_%s_%s_t%d_regstreamplot'%(os.path.basename(file_nc).replace('.',''), data_frommap_transx.name, data_frommap_transy.name,timestep)))

from dfm_tools.modplot import velovect
fig, ax = plt.subplots(1,1, figsize=(14,8))
quiv_curved = velovect(ax,X,Y,U,V, arrowstyle='fancy', scale = 5, grains = 25, color=speed)
cbar = fig.colorbar(quiv_curved.lines)
cbar.set_label('%s and %s (%s)'%(data_frommap_transx.name, data_frommap_transy.name, data_frommap_transy.attrs['units']))
ax.set_title('t=%d (%s)'%(timestep, data_frommap_transx.time.to_pandas().iloc[0]))
ax.set_xlim(xlim_get)
ax.set_ylim(ylim_get)
ax.set_aspect('equal')
fig.tight_layout()
plt.savefig(os.path.join(dir_output,'%s_%s_%s_t%d_curvedquiver'%(os.path.basename(file_nc).replace('.',''), data_frommap_transx.name, data_frommap_transy.name,timestep)))


