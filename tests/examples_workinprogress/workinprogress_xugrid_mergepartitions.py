# -*- coding: utf-8 -*-
"""
Created on Tue Jan 17 09:27:12 2023

@author: veenstra
"""

import os
#import xugrid as xu
#import xarray as xr
#xr.set_options(file_cache_maxsize=1)
import dfm_tools as dfmt

file_nc_list = ['p:\\1204257-dcsmzuno\\2006-2012\\3D-DCSM-FM\\A18b_ntsu1\\DFM_OUTPUT_DCSM-FM_0_5nm\\DCSM-FM_0_5nm_0*_map.nc', #3D DCSM
                'p:\\11206813-006-kpp2021_rmm-2d\\C_Work\\31_RMM_FMmodel\\computations\\model_setup\\run_207\\results\\RMM_dflowfm_0*_map.nc', #RMM 2D
                'p:\\1230882-emodnet_hrsm\\GTSMv5.0\\runs\\reference_GTSMv4.1_wiCA_2.20.06_mapformat4\\output\\gtsm_model_0*_map.nc', #GTSM 2D
                'p:\\11208053-005-kpp2022-rmm3d\\C_Work\\01_saltiMarlein\\RMM_2019_computations_02\\computations\\theo_03\\DFM_OUTPUT_RMM_dflowfm_2019\\RMM_dflowfm_2019_0*_map.nc', #RMM 3D
                'p:\\archivedprojects\\11203379-005-mwra-updated-bem\\03_model\\02_final\\A72_ntsu0_kzlb2\\DFM_OUTPUT_MB_02\\MB_02_0*_map.nc',
                ]

for file_nc in file_nc_list:
    print('processing %s'%(os.path.basename(file_nc)))
    
    ds_merged_xu = dfmt.open_partitioned_dataset(file_nc)#.replace('_0*_','_0000_'))
    
    
    #ds_merged_xu.load()
    #ds_merged_xu.mesh2d_sa1.mean(dim='time').compute() #to flood memory
