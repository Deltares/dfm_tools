# -*- coding: utf-8 -*-
"""
Created on Wed Oct 27 21:54:57 2021

@author: veenstra


to get delft3D to write netCDF output instead of .dat files, add these lines to your mdf:
    FlNcdf= #maphis#
    ncFormat=4

"""

import os
import numpy as np
import matplotlib.pyplot as plt
plt.close('all')
import xarray as xr
import dfm_tools as dfmt

#TODO: are corner2center and center2corner functions still needed?

dir_testinput = r'c:\DATA\dfm_tools_testdata'
dir_output = '.'

file_nc = r'p:\archivedprojects\1220688-lake-kivu\3_modelling\1_FLOW\7_heatfluxinhis\062_netcdf\trim-thiery_002_coarse.nc'
ds_thiery = xr.open_dataset(file_nc)
vars_pd = dfmt.get_ncvarproperties(ds_thiery)


def d3d_findnanval(data_nc_XZ,data_nc_YZ):
    values, counts = np.unique(data_nc_XZ, return_counts=True)
    X_nanval = values[np.argmax(counts)]
    values, counts = np.unique(data_nc_YZ, return_counts=True)
    Y_nanval = values[np.argmax(counts)]
    if X_nanval!=Y_nanval:
        XY_nanval = None
    else:
        XY_nanval = X_nanval
    return XY_nanval

def preprocess_delft3d4(ds):
    ds = ds.isel(MC=slice(None,-1),NC=slice(None,-1),M=slice(None,-1),N=slice(None,-1)) #last values are all nan/fill values
    
    #find and set nans in XZ/YZ arrays
    data_nc_XZ = ds.XZ
    data_nc_YZ = ds.YZ
    XY_nanval = d3d_findnanval(data_nc_XZ,data_nc_YZ)
    if XY_nanval is not None:
        mask_XY = (data_nc_XZ==XY_nanval) & (data_nc_YZ==XY_nanval)
        ds['XZ'] = data_nc_XZ.where(~mask_XY)
        ds['YZ'] = data_nc_YZ.where(~mask_XY)

    #find and set nans in XCOR/YCOR arrays
    data_nc_XCOR = ds.XCOR
    data_nc_YCOR = ds.YCOR
    XY_nanval = d3d_findnanval(data_nc_XCOR,data_nc_YCOR)
    if XY_nanval is not None:
        mask_XYCOR = (data_nc_XCOR==XY_nanval) & (data_nc_YCOR==XY_nanval)
        ds['XCOR'] = data_nc_XCOR.where(~mask_XYCOR)
        ds['YCOR'] = data_nc_YCOR.where(~mask_XYCOR)
    
    return ds

ds_thiery = preprocess_delft3d4(ds_thiery)
data_nc_XZ = ds_thiery.XZ
data_nc_YZ = ds_thiery.YZ
data_nc_XCOR = ds_thiery.XCOR
data_nc_YCOR = ds_thiery.YCOR

layno=-2


fig, ax = plt.subplots()
ax.plot(data_nc_XCOR,data_nc_YCOR,'-b',linewidth=0.2)
ax.plot(data_nc_XCOR.T,data_nc_YCOR.T,'-b',linewidth=0.2)
ax.set_aspect('equal')
dfmt.plot_coastlines(ax=ax)
plt.savefig(os.path.join(dir_output,'kivu_mesh'))

fig, axs = plt.subplots(1,3, figsize=(16,7))
for iT, timestep in enumerate([1,10,15]):
    ax=axs[iT]
    ds_thiery_tsel = ds_thiery.isel(KMAXOUT_RESTR=layno,time=timestep)
    U1 = ds_thiery_tsel.U1.to_numpy()
    V1 = ds_thiery_tsel.V1.to_numpy()
    ALFAS = ds_thiery_tsel.ALFAS.to_numpy() #contains rotation of all cells wrt real world
    timestr = ds_thiery_tsel.time.to_numpy()
    vel_x, vel_y, vel_magn, direction_naut_deg = dfmt.uva2xymagdeg(U1=U1,V1=V1,ALFAS=ALFAS)
    pc = ax.pcolor(data_nc_XCOR,data_nc_YCOR,vel_magn[1:,1:],cmap='jet')
    #pc = ax.pcolormesh(data_nc_XCOR,data_nc_YCOR,vel_magn[1:,1:],cmap='jet')
    pc.set_clim([0,0.15])
    cbar = fig.colorbar(pc, ax=ax)
    cbar.set_label('velocity magnitude (%s)'%(ds_thiery_tsel.U1.attrs['units']))
    ax.set_title(f't={timestep} ({timestr})')
    ax.set_aspect('equal')
    ax.quiver(data_nc_XZ[::2,::2], data_nc_YZ[::2,::2], vel_x[::2,::2], vel_y[::2,::2], 
              scale=3,color='w',width=0.005)
    ax.plot(data_nc_XCOR,data_nc_YCOR,'-b',linewidth=0.2)
    ax.plot(data_nc_XCOR.T,data_nc_YCOR.T,'-b',linewidth=0.2)
    dfmt.plot_coastlines(ax=ax)
fig.tight_layout()
plt.savefig(os.path.join(dir_output,'kivu_velocity_pcolor'))


fig, axs = plt.subplots(1,3, figsize=(16,7))
for iT, timestep in enumerate([1,10,15]):
    ax.plot(data_nc_XCOR,data_nc_YCOR,'-b',linewidth=0.2)
    ax.plot(data_nc_XCOR.T,data_nc_YCOR.T,'-b',linewidth=0.2)
    ax=axs[iT]
    ds_thiery_tsel = ds_thiery.isel(KMAXOUT_RESTR=layno,time=timestep)
    U1 = ds_thiery_tsel.U1.to_numpy()
    V1 = ds_thiery_tsel.V1.to_numpy()
    ALFAS = ds_thiery_tsel.ALFAS.to_numpy()
    timestr = ds_thiery_tsel.time.to_numpy()
    vel_x, vel_y, vel_magn, direction_naut_deg = dfmt.uva2xymagdeg(U1=U1,V1=V1,ALFAS=ALFAS) #TODO: clean up this function and make more convenient
    ax.set_title(f't={timestep} ({timestr})')
    ax.set_aspect('equal')
    pc = ax.quiver(data_nc_XZ[::2,::2], data_nc_YZ[::2,::2], vel_x[::2,::2], vel_y[::2,::2], vel_magn[::2,::2],
              scale=3,color='w',width=0.005, cmap='jet')
    pc.set_clim([0,0.15])
    cbar = fig.colorbar(pc, ax=ax)
    cbar.set_label('velocity magnitude (%s)'%(ds_thiery_tsel.U1.attrs['units']))
    dfmt.plot_coastlines(ax=ax)
fig.tight_layout()
plt.savefig(os.path.join(dir_output,'kivu_velocity'))

#QNET
fig, axs = plt.subplots(1,3, figsize=(16,7))
for iT, timestep in enumerate([1,10,15]):
    ax=axs[iT]
    ds_thiery_tsel = ds_thiery.isel(KMAXOUT_RESTR=layno,time=timestep)
    data_nc_QNET = ds_thiery_tsel.QNET.to_numpy()[1:,1:] #TODO: part of data is removed, not desireable. also using cor instead of cen.
    timestr = ds_thiery_tsel.time.to_numpy()
    pc = ax.pcolor(data_nc_XCOR,data_nc_YCOR,data_nc_QNET,cmap='jet')
    pc.set_clim([-60,60])
    cbar = fig.colorbar(pc, ax=ax)
    cbar.set_label('%s (%s)'%(ds_thiery_tsel.QNET.name, ds_thiery_tsel.QNET.attrs['units']))
    ax.set_title(f't={timestep} ({timestr})')
    ax.set_aspect('equal')
    ax.plot(data_nc_XCOR,data_nc_YCOR,'-b',linewidth=0.2)
    ax.plot(data_nc_XCOR.T,data_nc_YCOR.T,'-b',linewidth=0.2)
    dfmt.plot_coastlines(ax=ax)
fig.tight_layout()
plt.savefig(os.path.join(dir_output,'kivu_Qnet'))

#BED
fig, ax = plt.subplots(figsize=(6,8))

data_nc_DPV0 = ds_thiery.DPV0.to_numpy()[1:,1:] #TODO: part of data is removed, not desireable. also using cor instead of cen.
pc = ax.pcolor(data_nc_XCOR,data_nc_YCOR,data_nc_DPV0,cmap='jet')
cbar = fig.colorbar(pc, ax=ax)
cbar.set_label('%s (%s)'%(ds_thiery.DPV0.name, ds_thiery.DPV0.attrs['units']))
ax.set_aspect('equal')
ax.plot(data_nc_XCOR,data_nc_YCOR,'-b',linewidth=0.2)
ax.plot(data_nc_XCOR.T,data_nc_YCOR.T,'-b',linewidth=0.2)
dfmt.plot_coastlines(ax=ax)
fig.tight_layout()
plt.savefig(os.path.join(dir_output,'kivu_bedlevel'))



#from MAP DATA CURVEDBEND
file_nc = os.path.join(dir_testinput,'D3D_3D_sigma_curved_bend_nc\\trim-cb2-sal-added-3d.nc')
ds_cb2 = xr.open_dataset(file_nc)
vars_pd = dfmt.get_ncvarproperties(ds_cb2)

ds_cb2 = preprocess_delft3d4(ds_cb2)

data_nc_XZ = ds_cb2.XZ
data_nc_YZ = ds_cb2.YZ
data_nc_XCOR = ds_cb2.XCOR
data_nc_YCOR = ds_cb2.YCOR


fig, ax = plt.subplots()
ax.plot(data_nc_XCOR,data_nc_YCOR,'-b',linewidth=0.2)
ax.plot(data_nc_XCOR.T,data_nc_YCOR.T,'-b',linewidth=0.2)
#ax.plot(data_nc_XZ,data_nc_YZ,'-r',linewidth=0.2)
#ax.plot(data_nc_XZ.T,data_nc_YZ.T,'-r',linewidth=0.2)
ax.set_aspect('equal')
lim_x = [0,4100]
lim_y = [0,4100]
ax.set_xlim(lim_x)
ax.set_ylim(lim_y)
plt.savefig(os.path.join(dir_output,'curvedbend_mesh'))


txt_abcd = 'abcdefgh'
var_clim = [0,1.2]
ncols=2
fig, axs = plt.subplots(2,ncols, figsize=(8.6,8))
fig.suptitle('velocity magnitude on four times')
for iT, timestep in enumerate([0,1,2,4]):
    id0 = int(np.floor(iT/ncols))
    id1 = iT%ncols
    ax=axs[id0,id1]
    ds_cb2_tsel = ds_cb2.isel(KMAXOUT_RESTR=9,time=timestep)
    U1 = ds_cb2_tsel.U1.to_numpy()
    V1 = ds_cb2_tsel.V1.to_numpy()
    ALFAS = ds_cb2_tsel.ALFAS.to_numpy() #contains rotation of all cells wrt real world
    timestr = ds_cb2_tsel.time.to_pandas()
    vel_x, vel_y, vel_magn, direction_naut_deg = dfmt.uva2xymagdeg(U1=U1,V1=V1,ALFAS=ALFAS)
    pc = ax.pcolormesh(data_nc_XCOR,data_nc_YCOR,vel_magn[1:,1:],cmap='jet') #TODO: part of data is remove, not desireable
    pc.set_clim(var_clim)
    ax.set_aspect('equal')
    ax.quiver(data_nc_XZ, data_nc_YZ, vel_x, vel_y,
              scale=25,color='w',width=0.003)
    #add grid
    ax.plot(data_nc_XCOR,data_nc_YCOR,'-b',linewidth=0.2)
    ax.plot(data_nc_XCOR.T,data_nc_YCOR.T,'-b',linewidth=0.2)
    #additional figure formatting to tweak the details
    ax.grid(alpha=0.4)
    if id1 != 0:
        ax.get_yaxis().set_ticklabels([])
    else:
        ax.set_ylabel('y dist', labelpad=-0.5)
    if id0 == 0:
        ax.get_xaxis().set_ticklabels([])
    else:
        ax.set_xlabel('x dist')
    ax.tick_params(axis='x', labelsize=9)
    ax.tick_params(axis='y', labelsize=9)
    lim_x = [0,4100]
    lim_y = [0,4100]
    ax.set_xlim(lim_x)
    ax.set_ylim(lim_y)
    ax.text(lim_x[0]+0.02*np.diff(lim_x), lim_y[0]+0.95*np.diff(lim_y), '%s) t=%d (%s)'%(txt_abcd[iT],timestep, timestr),fontweight='bold',fontsize=12)
fig.tight_layout()
#additional figure formatting to tweak the details
plt.subplots_adjust(left=0.07, right=0.90, bottom=0.065, top=0.95, wspace=0.03, hspace=0.04)
cbar_ax = fig.add_axes([0.91, 0.065, 0.02, 0.885])
cbar = fig.colorbar(pc, cax=cbar_ax, ticks=np.linspace(var_clim[0],var_clim[1],7))
cbar_ax.set_ylabel('velocity magnitude [%s]'%(ds_cb2.U1.attrs['units']))
plt.savefig(os.path.join(dir_output,'curvedbend_velocity_pcolor'))


