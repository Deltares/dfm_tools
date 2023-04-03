# -*- coding: utf-8 -*-
"""
Created on Fri Oct  7 09:36:01 2022

@author: veenstra
"""

import os
from pathlib import Path
import matplotlib.pyplot as plt
plt.close('all')
import dfm_tools as dfmt
import hydrolib.core.dflowfm as hcdfm

dir_output = '.'

file_xyz = r'p:\archivedprojects\11205258-006-kpp2020_rmm-g6\C_Work\08_RMM_FMmodel\geometry_j19_NL_6-v2\rmm_vzm_v1p1_initial_water_level.xyz'
#file_xyz = r'p:\archivedprojects\11206813-006-kpp2021_rmm-2d\C_Work\31_RMM_FMmodel\general\diffusivity_rivzee_v6.xyz'

data_xyz = hcdfm.XYZModel(Path(file_xyz))
xyz_pd = dfmt.pointlike_to_DataFrame(data_xyz)

fig,ax = plt.subplots()
xyz_pd.plot.scatter(x='x',y='y',c='z',s=0.5,ax=ax)#,vmin=-1,vmax=1)
fig.tight_layout()
fig.savefig(os.path.join(dir_output,os.path.basename(file_xyz).replace('.','')))

file_xyn = r'p:\archivedprojects\11206813-006-kpp2021_rmm-2d\C_Work\31_RMM_FMmodel\geometry_j19_6-v2\output_locations\rmm_vzm-j19_6-v2b_3_measurement_obs.xyn'
file_xyn = r'c:\Users\veenstra\Downloads\test_obs.xyn' #contains: 1.23  3.45  'some text with spaces' And everything beyond does not matter \n 1.23 3.45 Only_the_first_word and the test is ignored
data_xyn1 = hcdfm.ObservationPointModel(Path(file_xyn)) #TODO: this should raise an error, but it returns an empy list: https://github.com/Deltares/HYDROLIB-core/issues/502
data_xyn2 = hcdfm.XYNModel(Path(file_xyn)) #TODO: test_obs.xyn should raise an error: https://github.com/Deltares/HYDROLIB-core/issues/459
data_xyn2_pd = dfmt.pointlike_to_DataFrame(data_xyn2) #works as expected #TODO: add DataFrame_to_xyn() conversion?
data_xyn2.save(file_xyn.replace('_obs.xyn','_out_obs.xyn')) #TODO: "TypeError: serialize() takes 3 positional arguments but 4 were given",  https://github.com/Deltares/HYDROLIB-core/issues/459

file_crs = r'p:\archivedprojects\11206813-006-kpp2021_rmm-2d\C_Work\31_RMM_FMmodel\geometry_j19_6-v2\cross_sections\rmm_vzm-j19_6-v2b_3_measurement_crs.pli'
data_crs1 = hcdfm.CrossLocModel(Path(file_crs)) #TODO: this should raise an error, but it returns an empty list: https://github.com/Deltares/HYDROLIB-core/issues/502
data_crs2 = hcdfm.CrossDefModel(Path(file_crs)) #TODO: https://github.com/Deltares/HYDROLIB-core/issues/364
data_crs3 = hcdfm.PolyFile(Path(file_crs)) #works with polyfile
polyobject_pd = [dfmt.pointlike_to_DataFrame(x) for x in data_crs3.objects] #TODO: how to get a useful+plottable object?
