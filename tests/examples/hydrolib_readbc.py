# -*- coding: utf-8 -*-
"""
Created on Wed Aug 17 11:19:51 2022

@author: veenstra
"""


#import os
from pathlib import Path
import datetime as dt
import pandas as pd
#from netCDF4 import num2date
#import cftime
import matplotlib.pyplot as plt
plt.close('all')
import numpy as np
from hydrolib.core.io.bc.models import ForcingModel
from dfm_tools.hydrolib_helpers import forcingobject_to_dataframe

#NOTE: for examples with writing bc files, check dfm_tools.interpolate_grid2bnd.*

if 0: #read in bc file with Timeseries objects (waterlevels, discharges)
    #file_bc = Path(r'p:\11208053-004-kpp2022-rmm1d2d\C_Work\09_Validatie2018_2020\dflowfm2d-rmm_vzm-j19_6-v2d\boundary_conditions\2020\flow\rmm_zeerand_v3_2020.bc') #>100 timeseries
    file_bc = Path(r'p:\11208053-004-kpp2022-rmm1d2d\C_Work\09_Validatie2018_2020\dflowfm2d-rmm_vzm-j19_6-v2d\boundary_conditions\rmm_rivdis_meas_20171220_20210102_MET.bc') #TODO: why can it not be str? #three timeseries
    #file_bc = Path(r'p:\11208053-004-kpp2022-rmm1d2d\C_Work\09_Validatie2018_2020\dflowfm2d-rmm_vzm-j19_6-v2d\boundary_conditions\2018\flow\rmm_discharge_laterals_20171220_20190101_MET.bc')
    
    #Load .bc-file using HydroLib object ForcingModel.
    m = ForcingModel(file_bc)
    
    """
    #print(m.forcing) #TODO: spyder crasht bij printen van m.forcing, 25% CPU en veel werkgeheugen, maar vooral een freeze. https://github.com/Deltares/HYDROLIB-core/issues/315
    type(m) #hydrolib.core.io.bc.models.ForcingModel
    type(m.forcing) #list of timeseries from bc file
    type(m.forcing[0]) # hydrolib.core.io.bc.models.TimeSeries
    type(m.forcing[0].__dict__) #dict
    m.forcing[0].__dict__.keys() # dict_keys(['comments', 'datablock', 'name', 'function', 'quantityunitpair', 'timeinterpolation', 'offset', 'factor'])
    type(m.forcing[0].__dict__['datablock']) #list #TODO: can this be retrieved more efficient or easier conversion to np.array?
    m.forcing[0].datablock #this is also possible without __dict__ method, but might be confusing to mix this and above method
    
    #TODO: reading improvement possible? (more intuitive coding) >> suggestion is below
    df = pd.DataFrame([f.__dict__ for f in m.forcing]) #this is a very heavy dataframe to interact with, but I guess it makes sense
    print(df[['datablock','quantityunitpair','name']])
    """
    
    df_data_list = [forcingobject_to_dataframe(forcingobj, convert_time=True) for forcingobj in m.forcing]
    
    #plot
    fig, ax = plt.subplots(figsize=(12, 6))
    for df_data in df_data_list[:5]:
        ax.plot(df_data, label=df_data.index.name,linewidth=0.7)
    ax.legend(loc=1)

    #calculate average waterlevel per point on boundary:
    df_data_all = pd.concat(df_data_list,axis=1)
    #list_names = [x.index.name for x in df_data_list]
    #list_colnamestr = [f"{x[0]} [{x[1]}]" for x in df_data_all.columns]
    mean_bndvals = df_data_all.mean(axis=0)
    
    #Plot mean waterlevel in BC-file over the complete period for each point on the boundary
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.plot([i for i in range(len(mean_bndvals))], mean_bndvals.values)
    ax.set_xlabel('Point on boundary')
    ax.set_ylabel(f"mean {mean_bndvals.index[0][0]} [{mean_bndvals.index[0][1]}]")


if 0: #read bc file with T3D blocks
    file_bc_3D = Path(r'n:\My Documents\werkmap\hydrolib_test\haixia\salinity_bc_simple_newnames.bc')
    #file_bc_3D = Path(r'n:\My Documents\werkmap\hydrolib_test\haixia\uxuy_bc_South_v2_firstpoint.bc')
    #file_bc_3D = Path(r'n:\My Documents\werkmap\hydrolib_test\haixia\salinity_bc_South_v2.bc')
    #file_bc_3D = Path(r'n:\My Documents\werkmap\hydrolib_test\haixia\uxuy_bc_South_v2.bc')
    #TODO: files from p:\11208154-002-haixia\02-hydrodynamics\02_Model_set_up\02_Make_forcing\CMEMS\bc\V2\Southv2\, bc files for other boundaries are 4/5 times bigger
    
    m = ForcingModel(file_bc_3D) #TODO: crashes on validation error since [verticalPositions,verticalInterpolation,verticalPositionType] are not accepted as missing (while dflowfm has default values for them) >> verticalPositions are in file as "Vertical position". Issue created: https://github.com/Deltares/HYDROLIB-core/issues/306


if 0: #read bc file with harmonic blocks and compare
    
    file_bc_orig = Path(r'n:\My Documents\werkmap\hydrolib_test\DCSM\tide_OB_all_20181108.bc')
    file_bc_new = Path(r'n:\My Documents\werkmap\hydrolib_test\DCSM\tide_DCSM-FM_OB_all_20181108_nocomments.bc')
    
    iPoint = 1
    
    list_pd = []
    for iF,file_bc in enumerate([file_bc_orig,file_bc_new]):
        #Load .bc-file using HydroLib object ForcingModel.
        m = ForcingModel(file_bc) #TODO: this takes a while, can it be spead up? (does all fit in memory? maybe only read when requested like netcdf variable or inconvenient?)
    
        df = forcingobject_to_dataframe(m.forcing[iPoint], convert_time=False)
        name = df.index.name
        df = df.droplevel(1,axis=1)
        df = df.set_index('astronomic component')
        df.index.name = name
        list_pd.append(df)
    
    print(list_pd[0]-list_pd[1])
    
