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


uds_withcellinfo = dfmt.add_network_cellinfo(uds)

assert isinstance(uds.grid, xu.Ugrid1d)
assert isinstance(uds_withcellinfo.grid, xu.Ugrid2d)
assert hasattr(uds_withcellinfo.grid, "face_node_connectivity")

fig,ax = plt.subplots()
print(uds.grid.to_dataset().Mesh2D) # 1D topology has to be converted to 2D
uds.grid.plot(ax=ax,linewidth=0.5,color='crimson')
ax.set_title('1D network')

fig,ax = plt.subplots()
uds_withcellinfo.grid.plot(ax=ax,linewidth=0.5,color='crimson')
ax.set_title('2D network')

