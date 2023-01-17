# -*- coding: utf-8 -*-
"""
Created on Mon Jan 16 12:09:43 2023

@author: veenstra
"""

import xarray as xr

file_nc = r'p:\1204257-dcsmzuno\2006-2012\3D-DCSM-FM\A18b_ntsu1\DFM_OUTPUT_DCSM-FM_0_5nm\DCSM-FM_0_5nm_0000_map.nc'
#file_nc = r'p:\11206813-006-kpp2021_rmm-2d\C_Work\31_RMM_FMmodel\computations\model_setup\run_207\results\RMM_dflowfm_0000_map.nc'
#file_nc = r'p:\archivedprojects\11203379-005-mwra-updated-bem\03_model\02_final\A72_ntsu0_kzlb2\DFM_OUTPUT_MB_02\MB_02_0000_map.nc'
#file_nc = r'c:\DATA\MB_02_0000_map.nc'

ds = xr.open_dataset(file_nc)
