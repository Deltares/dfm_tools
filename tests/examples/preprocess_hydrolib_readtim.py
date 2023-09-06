# -*- coding: utf-8 -*-
"""
Created on Tue Apr 11 17:45:54 2023

@author: veenstra
"""

import dfm_tools as dfmt
import hydrolib.core.dflowfm as hcdfm

file_tim = r'p:\archivedprojects\11205258-006-kpp2020_rmm-g6\C_Work\08_RMM_FMmodel\structures_toRTC\Algerakering_GateLowerEdgeLevel.tim'
file_tim = r'p:\archivedprojects\11206811-002-d-hydro-grevelingen\simulaties\Jaarsom2017_dfm_006_zlayer\boundary_conditions\hist\jaarsom_2017\sources_sinks\FlakkeeseSpuisluis.tim'

data_tim = hcdfm.TimModel(file_tim)

refdate = '2007-01-01'
tim_pd = dfmt.TimModel_to_DataFrame(data_tim, parse_column_labels=True, refdate=refdate)

data_tim_out = dfmt.DataFrame_to_TimModel(tim_pd, refdate=refdate)
#data_tim_out.save(file_tim.replace('.tim','_out.tim'))
