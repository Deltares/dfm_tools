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
try: #0.3.1 release
    from hydrolib.core.io.xyz.models import XYZModel
    from hydrolib.core.io.obs.models import ObservationPointModel
    from hydrolib.core.io.crosssection.models import CrossLocModel, CrossDefModel
except: #main branch and next release #TODO: move to easy imports after https://github.com/Deltares/HYDROLIB-core/issues/410
    from hydrolib.core.io.dflowfm.xyz.models import XYZModel
    from hydrolib.core.io.dflowfm.obs.models import ObservationPointModel
    from hydrolib.core.io.dflowfm.crosssection.models import CrossLocModel, CrossDefModel

dir_output = '.'

file_xyz = r'p:\11205258-006-kpp2020_rmm-g6\C_Work\08_RMM_FMmodel\geometry_j19_NL_6-v2\rmm_vzm_v1p1_initial_water_level.xyz'
#file_xyz = r'p:\11206813-006-kpp2021_rmm-2d\C_Work\31_RMM_FMmodel\general\diffusivity_rivzee_v6.xyz' #TODO: xyz file with initial whitespaces is not read properly: https://github.com/Deltares/HYDROLIB-core/issues/415

data_xyz = XYZModel(Path(file_xyz))
xyz_pd = dfmt.pointlike_to_DataFrame(data_xyz)

fig,ax = plt.subplots()
xyz_pd.plot.scatter(x='x',y='y',c='z',s=0.5,ax=ax,vmin=-1,vmax=1)
fig.tight_layout()
fig.savefig(os.path.join(dir_output,os.path.basename(file_xyz).replace('.','')))


#TODO: below does not work, 'old' xyn/crs is not supported: https://github.com/Deltares/HYDROLIB-core/issues/364

file_xyn = r'p:\11206813-006-kpp2021_rmm-2d\C_Work\31_RMM_FMmodel\geometry_j19_6-v2\output_locations\rmm_vzm-j19_6-v2b_3_measurement_obs.xyn'
data_xyn = ObservationPointModel(Path(file_xyn))

file_crs = r'p:\11206813-006-kpp2021_rmm-2d\C_Work\31_RMM_FMmodel\geometry_j19_6-v2\cross_sections\rmm_vzm-j19_6-v2b_3_measurement_crs.pli'
data_crs1 = CrossLocModel(Path(file_crs))
data_crs2 = CrossDefModel(Path(file_crs))
