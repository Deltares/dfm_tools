# -*- coding: utf-8 -*-
"""
Created on Tue Jan 24 10:42:34 2023

@author: veenstra
"""


import os
import matplotlib.pyplot as plt
plt.close('all')
import dfm_tools as dfmt
import matplotlib as mpl
import numpy as np

file_nc = os.path.join(r'c:\DATA\dfm_tools_testdata','DFM_3D_z_Grevelingen','computations','run01','DFM_OUTPUT_Grevelingen-FM','Grevelingen-FM_0*_map.nc')

cmap = mpl.cm.get_cmap("viridis")
clim_bl = [-40,10]
start,stop = clim_bl
step = 3
bounds = [i for i in np.arange(start, stop +1E-6, step)]
boundaries = np.array(bounds)
norm = mpl.colors.BoundaryNorm(boundaries, cmap.N) #setting this to None results in equal plots, when not None the ugrid plot is wrong

data_frommap_merged = dfmt.open_partitioned_dataset(file_nc)
ugrid_all_verts = dfmt.get_ugrid_verts(data_frommap_merged)

#get bedlevel and create plot with ugrid and cross section line
fig, (ax1,ax2) = plt.subplots(2,1,figsize=(7,8))

pc1 = mpl.collections.PolyCollection(ugrid_all_verts, linewidth=0.5, edgecolors='face', cmap=cmap, norm=norm)
pc1.set_array(data_frommap_merged['mesh2d_flowelem_bl'].to_numpy())
ax1.add_collection(pc1)
ax1.autoscale() #necessary to call manually
fig.colorbar(pc1,ax=ax1, extend = 'max')
pc1.set_clim(clim_bl)
ax1.set_aspect('equal')

pc2 = data_frommap_merged['mesh2d_flowelem_bl'].ugrid.plot(linewidth=0.5, edgecolor='face',cmap=cmap, ax=ax2, norm=norm, add_colorbar=False)
fig.colorbar(pc2,ax=ax2, extend = 'max')
pc2.set_clim(clim_bl)
ax2.set_aspect('equal')
