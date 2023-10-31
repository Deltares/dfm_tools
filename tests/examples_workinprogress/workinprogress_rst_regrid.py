# -*- coding: utf-8 -*-
"""
Created on Tue Oct 31 10:22:33 2023

@author: veenstra
"""


import dfm_tools as dfmt
import matplotlib.pyplot as plt
plt.close('all')
import xugrid as xu

file_nc = r"p:\11203618-uae-wqps\2_models\3_WAQ\03_Simulations\05_OptimizationOctoberCoarse\05_RefTestSmartLoads\DFM_OUTPUT_AGM3D\AGM3D_0*_20210102_000000_rst.nc"

uds_rst = dfmt.open_partitioned_dataset(file_nc, preprocess=dfmt.enrich_rst_with_map)

# fig,ax = plt.subplots()
# uds_rst.s1.isel(time=0).ugrid.plot()

#open fine network and convert from Ugrid1d to Ugrid2d
file_net_fine = r"p:\11203618-uae-wqps\2_models\4_models_to_FEWS\preVisit2\ArabianGulf\1_FM\AGM3DWQ\agm_v20231015a_net.nc"
uds_fine = xu.open_dataset(file_net_fine)

uds_fine_2d = dfmt.add_network_cellinfo(uds_fine)

# #plot grids
# fig,ax = plt.subplots()
# uds_fine_2d.grid.plot(ax=ax)
# uds.grid.plot(ax=ax, color='r')

uds_regridder = xu.CentroidLocatorRegridder(source=uds_rst.grid, target=uds_fine_2d.grid)
uda_regrid = uds_regridder.regrid(uds_rst.sa1)

# plot regridded result
fig,ax = plt.subplots()
uda_regrid.isel(time=-1, laydim=-5).ugrid.plot(ax=ax, robust=True, vmax=50)
