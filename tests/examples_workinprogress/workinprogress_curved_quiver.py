# -*- coding: utf-8 -*-
"""
Created on Thu Mar 30 22:06:29 2023

@author: veenstra
"""

import datetime as dt
import numpy as np
import matplotlib.pyplot as plt
#plt.close('all')
import xarray as xr
import dfm_tools as dfmt

testdata = True
if testdata:
    x = np.linspace(-4,4,120)
    y = np.linspace(-3,3,100)
    X,Y = np.meshgrid(x,y)
    U = -1 - X**2 + Y
    V = 1 + X - Y**2
else:
    ds = xr.tutorial.load_dataset("eraint_uvz").isel(month=0,level=0)
    ds = ds.sortby('latitude') #.sortby() to avoid "ValueError: 'y' must be strictly increasing"
    x = ds.longitude.to_numpy()
    y = ds.latitude.to_numpy()
    X,Y = np.meshgrid(x,y)
    U = ds.u.to_numpy()
    V = ds.v.to_numpy()
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
                     scale=5, grains=grains)
print(f'{(dt.datetime.now()-dtstart).total_seconds():.2f} sec')
fig.colorbar(strm.lines)

plt.tight_layout()
plt.show()