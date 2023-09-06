# -*- coding: utf-8 -*-
"""
Created on Thu Mar 30 22:06:29 2023

@author: veenstra
"""

import os
import datetime as dt
import numpy as np
import matplotlib.pyplot as plt
plt.close('all')
import xarray as xr
import dfm_tools as dfmt

dir_output = '.'

mode = 'gtsm'
if mode=='testdata':
    x = np.linspace(-4,4,120)
    y = np.linspace(-3,3,100)
    X,Y = np.meshgrid(x,y)
    U = -1 - X**2 + Y
    V = 1 + X - Y**2
elif mode=='eraint_uvz':
    ds = xr.tutorial.load_dataset("eraint_uvz").isel(month=0,level=0)
    ds = ds.sortby('latitude') #.sortby() to avoid "ValueError: 'y' must be strictly increasing"
    x = ds.longitude.to_numpy()
    y = ds.latitude.to_numpy()
    X,Y = np.meshgrid(x,y)
    U = ds.u.to_numpy()
    V = ds.v.to_numpy()
elif mode=='gtsm':
    file_nc = r'p:\1230882-emodnet_hrsm\GTSMv5.0\runs\reference_GTSMv4.1_wiCA\output\gtsm_model_0000_map.nc'
    drop_vars = ['waterdepth','TidalPotential_without_SAL','SALPotential','internal_tides_dissipation','czs','czu','FlowElemContour_x','FlowElemContour_y']
    uds = dfmt.open_partitioned_dataset(file_nc,drop_variables=drop_vars)
    uds_sel = uds.sel(time='2014-01-31 10:00:00').ugrid.sel(x=slice(-3,10),y=slice(48,58))
    ds = dfmt.rasterize_ugrid(uds_sel,resolution=0.1)
    x = ds.x.to_numpy()
    y = ds.y.to_numpy()
    X,Y = np.meshgrid(x,y)
    U = ds.ucx.to_numpy()
    V = ds.ucy.to_numpy()
    #plot
    fig,ax = plt.subplots()
    uds_sel_speed = np.sqrt(uds_sel.ucx*uds_sel.ucx + uds_sel.ucy*uds_sel.ucy)
    uds_sel_speed.ugrid.plot(cmap='jet',vmax=1)
    #ax.quiver(X, Y, U, V, color='w')
    speed = np.sqrt(U*U + V*V)
    strm = dfmt.velovect(ax, X, Y, U, V, color='w', cmap='winter', arrowstyle='fancy', linewidth=speed*2, integration_direction='forward',
                         density=5, grains=30, arrowsize=0.7)
speed = np.sqrt(U*U + V*V)

grains = 15

def get_start_points(x,y,grains:(tuple,int) = 15): #TODO: is comparable to dfmt.modplot._gen_starting_points()
    if isinstance(grains,tuple):
        nx, ny = grains
    elif isinstance(grains,int):
        nx = ny = grains
    x_coarse = np.linspace(x.min(), x.max(), nx)
    y_coarse = np.linspace(y.min(), y.max(), ny)
    X_coarse,Y_coarse = np.meshgrid(x_coarse,y_coarse)
    xs = X_coarse.ravel()
    ys = Y_coarse.ravel()
    start_points = np.c_[xs, ys]
    return start_points
start_points = get_start_points(x=x,y=y,grains=grains)

# Varying color along a streamline
fig,(ax1,ax2,ax3) = plt.subplots(1,3,figsize=(14,5),sharex=True,sharey=True)

strm = ax1.streamplot(X, Y, U, V, color=speed, cmap='winter', arrowstyle='fancy', linewidth=speed/5, integration_direction='both')
#strm = ds.plot.streamplot(ax=ax1, x='longitude',y='latitude',u='u',v='v',color=speed) #TODO: eventually do it with xarray directly
fig.colorbar(strm.lines)

strm = ax2.streamplot(X, Y, U, V, color=speed, cmap='winter', arrowstyle='fancy', linewidth=speed/5, integration_direction='both',
                      density=5, #density=1 default, higher gives equally spaced view
                      minlength=0.01, maxlength = 0.07,
                      start_points=start_points)
fig.colorbar(strm.lines)


print('>> dfmt.velovect(): ',end='')
dtstart = dt.datetime.now()
strm = dfmt.velovect(ax3, X, Y, U, V, color=speed, cmap='winter', arrowstyle='fancy', linewidth=speed/5, integration_direction='forward',
                     density=5, grains=grains)
print(f'{(dt.datetime.now()-dtstart).total_seconds():.2f} sec')
fig.colorbar(strm.lines)

fig.tight_layout()
plt.savefig(os.path.join(dir_output,'curved_quiver'))
