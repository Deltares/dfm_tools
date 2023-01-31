# -*- coding: utf-8 -*-
"""
Created on Tue Jan 31 13:57:55 2023

@author: veenstra
"""

#https://github.com/pydata/xarray/issues/7014
#https://github.com/pydata/xarray/issues/4061


import matplotlib.pyplot as plt
plt.close('all')
import matplotlib as mpl
import numpy as np
import xarray as xr

cmap = mpl.cm.get_cmap("viridis")
clim_air = [227,302]
start,stop = clim_air
boundaries = np.array([i for i in np.arange(start, stop, step=4)])
norm = None #this works for both plots
norm = mpl.colors.BoundaryNorm(boundaries, cmap.N) #this fails for da.plot(), but works for ax.pcolormesh()
norm = mpl.colors.BoundaryNorm(boundaries, ncolors=len(boundaries)-1) #this works for da.plot(), but fails for ax.pcolormesh()
norm = mpl.colors.BoundaryNorm(boundaries, ncolors=100) #ncolors value in between, both methods fail.

ds = xr.tutorial.load_dataset("air_temperature")
da = ds.air.isel(time=0)

#get bedlevel and create plot with ugrid and cross section line
fig, (ax1,ax2) = plt.subplots(2,1,figsize=(6,6))

pc1 = da.plot(cmap=cmap, ax=ax1, norm=norm, add_colorbar=False)#, vmin=norm.vmin,vmax=norm.vmax)
fig.colorbar(pc1,ax=ax1, extend = 'max')
#pc1.set_clim(clim_air)
ax1.set_title('da.plot()')

pc2 = ax2.pcolormesh(ds.lon,ds.lat,da.values, cmap=cmap, norm=norm)
fig.colorbar(pc2,ax=ax2, extend = 'max')
#pc2.set_clim(clim_air)
ax2.set_title('ax.pcolormesh(x,y,da)')

fig.tight_layout()
