# -*- coding: utf-8 -*-
"""
Created on Wed Aug 17 11:19:51 2022

@author: veenstra
"""

#import os
from pathlib import Path
#import datetime as dt
import matplotlib.pyplot as plt
plt.close('all')
from hydrolib.core.io.bc.models import ForcingModel
from dfm_tools.hydrolib_helpers import forcinglike_to_Dataset, forcinglike_to_DataFrame, Dataset_to_TimeSeries, Dataset_to_T3D, Dataset_to_T3Dvector#, Dataset_to_Astronomic

#TODO: merge this into hydrolib_readFMmodel.py after issues are resolved?
#NOTE: for examples with writing bc files, check dfm_tools.interpolate_grid2bnd.* and dfm_tools.hydrolib_helpers.

nPoints = 5 #None for all points

file_bc_list = [Path(r'n:\My Documents\werkmap\hydrolib_test\DCSM\tide_OB_all_20181108.bc'), #TODO: make better DataArray format for astronomic components
                #Path(r'p:\11208053-004-kpp2022-rmm1d2d\C_Work\09_Validatie2018_2020\dflowfm2d-rmm_vzm-j19_6-v2d\boundary_conditions\2020\flow\rmm_zeerand_v3_2020.bc'), #>100 timeseries
                Path(r'n:\My Documents\werkmap\hydrolib_test\rmm_zeerand_v3_2020_short3.bc'),
                Path(r'p:\11208053-004-kpp2022-rmm1d2d\C_Work\09_Validatie2018_2020\dflowfm2d-rmm_vzm-j19_6-v2d\boundary_conditions\rmm_rivdis_meas_20171101_20210102_MET.bc'), #TODO: why can it not be str? #three timeseries
                #Path(r'p:\11208053-004-kpp2022-rmm1d2d\C_Work\09_Validatie2018_2020\dflowfm2d-rmm_vzm-j19_6-v2d\boundary_conditions\2018\flow\rmm_discharge_laterals_20171201_20190101_MET.bc'),
                Path(r'n:\My Documents\werkmap\hydrolib_test\haixia\salinity_bc_South_v2_firstpoint.bc'),
                Path(r'n:\My Documents\werkmap\hydrolib_test\haixia\uxuy_bc_South_v2_firstpoint.bc'), #TODO SOLVED: not possible to read uxuy yet: https://github.com/Deltares/HYDROLIB-core/issues/316
                #Path(r'n:\My Documents\werkmap\hydrolib_test\haixia\salinity_bc_South_v2.bc'), #large file, takes time
                #Path(r'n:\My Documents\werkmap\hydrolib_test\haixia\uxuy_bc_South_v2.bc'),
                ]
file_bc_list = [#Path(r'n:\My Documents\werkmap\hydrolib_test\DCSM\tide_OB_all_20181108.bc'), #TODO: make better DataArray format for astronomic components
                Path(r'p:\11208053-004-kpp2022-rmm1d2d\C_Work\09_Validatie2018_2020\dflowfm2d-rmm_vzm-j19_6-v2d\boundary_conditions\rmm_rivdis_meas_20171101_20210102_MET.bc'), #TODO: why can it not be str? #three timeseries
                #Path(r'n:\My Documents\werkmap\hydrolib_test\haixia\salinity_bc_South_v2_firstpoint.bc'),
                #Path(r'n:\My Documents\werkmap\hydrolib_test\haixia\uxuy_bc_South_v2_firstpoint_short.bc'), #TODO SOLVED: not possible to read uxuy yet: https://github.com/Deltares/HYDROLIB-core/issues/316
                ]


for file_bc in file_bc_list:
    #Load .bc-file using HydroLib object ForcingModel.
    m = ForcingModel(file_bc)
    ForcingModel_object_out = ForcingModel()
    
    # m.general.comments = {'a':'aa'} #TODO: adding comments to top of file is not possible, only if using filetype or fileversion: https://github.com/Deltares/HYDROLIB-core/issues/130. Top file comment newfeature: https://github.com/Deltares/HYDROLIB-core/issues/362
    # m.save('test.bc')
    
    """
    type(m) #hydrolib.core.io.bc.models.ForcingModel
    type(m.forcing) #list of timeseries from bc file
    type(m.forcing[0]) # hydrolib.core.io.bc.models.TimeSeries
    type(m.forcing[0].__dict__) #dict
    m.forcing[0].__dict__.keys() # dict_keys(['comments', 'datablock', 'name', 'function', 'quantityunitpair', 'timeinterpolation', 'offset', 'factor'])
    type(m.forcing[0].datablock) #list
    """
    
    #plot
    fig, ax = plt.subplots(figsize=(12, 6))
    for iFO, forcingobj in enumerate(m.forcing[:nPoints]):
        forcing_xr = forcinglike_to_Dataset(forcingobj)
        if forcingobj.function=='t3d':
            if hasattr(forcingobj.quantityunitpair[1],'elementname'): #uxuy vector
                forcing_ts = Dataset_to_T3Dvector(forcing_xr)
            else:
                forcing_ts = Dataset_to_T3D(forcing_xr)
            for iV,varkey in enumerate(forcing_xr.data_vars): #can be multiple data_vars
                if iV>=1: #create new figure in case of multiple variables
                    fig, ax = plt.subplots(figsize=(12, 6))
                forcing_xr[varkey].T.plot(ax=ax)
            #forcing_xr_masked = forcing_xr.where(forcing_xr != forcing_xr.isel(depth=0)) #TODO: better way to mask data?
            #fig, ax = plt.subplots(figsize=(12, 6))
            #forcing_xr_masked.T.plot()
        elif forcingobj.function=='timeseries':
            forcing_ts = Dataset_to_TimeSeries(forcing_xr)
            for iV,varkey in enumerate(forcing_xr.data_vars): #always one data_var
                forcing_xr[varkey].plot(ax=ax, label=forcing_xr[varkey].attrs['name'], linewidth=0.8)
            ax.legend(loc=1)
        elif forcingobj.function=='astronomic':
            #TODO: also implement Dataset_to_Astronomic() definition
            if iFO==0:
                ax2 = ax.twinx()
                continue #skip first (invalid) point
            data_vars = list(forcing_xr.data_vars)
            forcing_xr[data_vars[0]].plot(ax=ax, linewidth=0.7)
            forcing_xr[data_vars[1]].plot(ax=ax2, linewidth=0.7)
            #ax.legend(loc=1)
            # forcing_pd = forcinglike_to_DataFrame(forcingobj) #TODO: only relevant for Astronomic, for TimeSeries/T3D it is equal to xr.to_pandas() (phase out this definition)
            # forcing_pd[forcing_pd.columns[0]].plot(ax=ax, linewidth=0.7)
            # forcing_pd[forcing_pd.columns[1]].plot(ax=ax, linewidth=0.7, secondary_y=True)
            # ax.legend(forcing_pd.columns,loc=1)
        else:
            forcing_xr.plot(ax=ax, label=forcing_xr.attrs['name'], linewidth=0.7)
            ax.legend(loc=1)
            raise Exception(f'non-defined function: {forcingobj.function}')
        #ForcingModel_object_out.forcing.append(forcing_ts)
        #ForcingModel_object_out.save(Path(str(file_bc).replace('.bc','_reproduced.bc')))
       


