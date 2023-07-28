# -*- coding: utf-8 -*-
"""
Created on Fri Jul 28 16:05:10 2023

@author: veenstra
"""

import os
import matplotlib.pyplot as plt
plt.close('all')
import xugrid as xu
import dfm_tools as dfmt


dir_issue = r'p:\dflowfm\maintenance\JIRA\06000-06999\06548'

file_nc = os.path.join(dir_issue,'set_bathy_2.17.06','baltic_nobathy_net.nc') #1D/minimal network

#uds = xu.open_dataset(file_nc) #with new grid with 2D topology: "ValueError: edge_node_connectivity is invalid, the following edges are not associated with any face"
uds = dfmt.open_partitioned_dataset(file_nc)

#check if indeed 1D grid object
assert isinstance(uds.grid, xu.ugrid.ugrid1d.Ugrid1d)

mk = uds.grid.meshkernel #TODO: is_geographic=False is currently hardcoded: https://github.com/Deltares/xugrid/issues/128

# uds.ugrid.sel(x=slice_x,y=slice_x) 
fig,ax = plt.subplots()
print(uds.grid.to_dataset().Mesh2D) # 1D topology has to be converted to 2D
uds.grid.plot(ax=ax,linewidth=0.5,color='crimson')
ax.set_title('1D network')

# compute 2D topology
import meshkernel
mk_mesh1d = mk.mesh1d_get()
input_mesh2d = meshkernel.Mesh2d(mk_mesh1d.node_x, mk_mesh1d.node_y, mk_mesh1d.edge_nodes)

# simple approach, but results in cartesian grid
# mk.mesh2d_set(input_mesh2d)
# uds_ugrid2d = xu.Ugrid2d.from_meshkernel(mk.mesh2d_get())

#now with is_geographic=True and EPSG
mk2 = meshkernel.MeshKernel(is_geographic=True)
mk2.mesh2d_set(input_mesh2d)
xu_grid_uds = dfmt.meshkernel_to_UgridDataset(mk=mk2, crs='EPSG:4326')
# xu_grid_uds.ugrid.to_netcdf('out_net.nc')