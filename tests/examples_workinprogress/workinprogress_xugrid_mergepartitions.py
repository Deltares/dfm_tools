# -*- coding: utf-8 -*-
"""
Created on Tue Jan 17 09:27:12 2023

@author: veenstra
"""

import os
import datetime as dt
import glob
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
    dtstart_all = dt.datetime.now()
    if isinstance(file_nc,list):
        file_nc_list = file_nc
    else:
        file_nc_list = glob.glob(file_nc)
    if len(file_nc_list)==0:
        raise Exception('file(s) not found, empty file_nc_list')
    
    ds_merged_xu = dfmt.open_partitioned_dataset(file_nc_list, decode_times_perfile=False)
    
    # print(f'>> xr.open_dataset() with {len(file_nc_list)} partition(s): ',end='')
    # dtstart = dt.datetime.now()
    # partitions = []
    # for iF, file_nc_one in enumerate(file_nc_list[:4]):
    #     print(iF+1,end=' ')
    #     uds = xu.open_dataset(file_nc_one, chunks={'time':1},decode_times=False)#, cache=False)
    #     partitions.append(uds)
    # print(': ',end='')
    # print(f'{(dt.datetime.now()-dtstart).total_seconds():.2f} sec')
    
    # print(f'>> xu.merge_partitions() with {len(file_nc_list)} partition(s): ',end='')
    # dtstart = dt.datetime.now()
    # ds_merged_xu = xu.merge_partitions(partitions)
    # print(f'{(dt.datetime.now()-dtstart).total_seconds():.2f} sec')
    
    #ds_merged_xu.load()
    #ds_merged_xu.mesh2d_sa1.mean(dim='time').compute() #to flood memory
