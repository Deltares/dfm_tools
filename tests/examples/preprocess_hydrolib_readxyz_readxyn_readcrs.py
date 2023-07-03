# -*- coding: utf-8 -*-
"""
Created on Fri Oct  7 09:36:01 2022

@author: veenstra
"""

import os
import matplotlib.pyplot as plt
plt.close('all')
import dfm_tools as dfmt
import hydrolib.core.dflowfm as hcdfm

dir_output = '.'

file_xyz = r'p:\archivedprojects\11205258-006-kpp2020_rmm-g6\C_Work\08_RMM_FMmodel\geometry_j19_NL_6-v2\rmm_vzm_v1p1_initial_water_level.xyz'
#file_xyz = r'p:\archivedprojects\11206813-006-kpp2021_rmm-2d\C_Work\31_RMM_FMmodel\general\diffusivity_rivzee_v6.xyz'
data_xyz = hcdfm.XYZModel(file_xyz)
fig,ax = plt.subplots()
xyz_gpd = dfmt.pointlike_to_geodataframe_points(data_xyz, crs=None) #TODO: make crs=None default?
xyz_gpd_sel = xyz_gpd.cx[55000:83000,432000:450000]
xyz_gpd_sel.plot(ax=ax,column='z',markersize=0.5,cmap='jet',legend=True)
fig.savefig(os.path.join(dir_output,os.path.basename(file_xyz).replace('.','')))


file_xyn = r'p:\archivedprojects\11206813-006-kpp2021_rmm-2d\C_Work\31_RMM_FMmodel\geometry_j19_6-v2\output_locations\rmm_vzm-j19_6-v2b_3_measurement_obs.xyn'
data_xyn1 = hcdfm.ObservationPointModel(file_xyn) #TODO: this should raise an error, but it returns an empty list: https://github.com/Deltares/HYDROLIB-core/issues/502
data_xyn2 = hcdfm.XYNModel(file_xyn)
data_xyn2_gpd = dfmt.pointlike_to_geodataframe_points(data_xyn2,crs=None)
data_xyn2_gpd.plot()


file_crs = r'p:\archivedprojects\11206813-006-kpp2021_rmm-2d\C_Work\31_RMM_FMmodel\geometry_j19_6-v2\cross_sections\rmm_vzm-j19_6-v2b_3_measurement_crs.pli'
data_crs1 = hcdfm.CrossLocModel(file_crs) #TODO: this should raise an error, but it returns an empty list: https://github.com/Deltares/HYDROLIB-core/issues/502
data_crs2 = hcdfm.CrossDefModel(file_crs) #TODO: same
data_crs3 = hcdfm.PolyFile(file_crs) #works with polyfile
polyobject_gpd = dfmt.PolyFile_to_geodataframe_linestrings(data_crs3,crs=None)
polyobject_gpd.plot(color='r')
