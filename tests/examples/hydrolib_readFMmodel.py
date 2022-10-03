# -*- coding: utf-8 -*-
"""
Created on Mon Oct  3 12:07:18 2022

@author: veenstra
"""

from pathlib import Path
from hydrolib.core.io.mdu.models import FMModel, NetworkModel, ExtModel

file_mdu = Path(r'p:\11206813-006-kpp2021_rmm-2d\C_Work\31_RMM_FMmodel\computations\model_setup\run_206_HYDROLIB\RMM_dflowfm.mdu') #model with all but one structure and all but one lateral commented, reduces validation errors from >200 to 5. TODO: resolve validation errors
#file_mdu = Path(r'c:\DATA\dfm_tools_testdata\DFM_3D_z_Grevelingen\computations\run01\Grevelingen-FM.mdu')

print('reading fmmodel')
fm = FMModel(file_mdu)
print('done')


#ext = ExtModel(fm.external_forcing.extforcefilenew)

#dimr = DIMR("dimr_config.xml")
