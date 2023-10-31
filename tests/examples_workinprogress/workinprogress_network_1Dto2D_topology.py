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
file_nc = os.path.join(dir_issue,'set_bathy_2.17.06','baltic_nobathy_net.nc') #1D/minimal network, with hanging edges
# file_nc = r"p:\dflowfm\maintenance\JIRA\06000-06999\06548\set_bathy_2.25.07\DFM_interpreted_network_dummy_model_net.nc"

uds = xu.open_dataset(file_nc)

uds_withcellinfo = dfmt.add_network_cellinfo(uds)

# fig,ax = plt.subplots()
# uds.grid.plot(ax=ax,linewidth=0.02)
# ax.set_title('1D network')

# fig,ax = plt.subplots()
# uds_withcellinfo.grid.plot(ax=ax,linewidth=0.02)
# ax.set_title('2D network')

fig,ax = plt.subplots()
uds.NetNode_z.ugrid.plot(ax=ax)
ax.set_title('1D network')

fig,ax = plt.subplots()
uds_withcellinfo.NetNode_z.ugrid.plot(ax=ax)
ax.set_title('2D network')

