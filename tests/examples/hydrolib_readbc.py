# -*- coding: utf-8 -*-
"""
Created on Wed Aug 17 11:19:51 2022

@author: veenstra
"""
#TODO: merge this into hydrolib_readFMmodel.py after issues are resolved

#import os
from pathlib import Path
import matplotlib.pyplot as plt
plt.close('all')
from hydrolib.core.io.bc.models import ForcingModel
from dfm_tools.hydrolib_helpers import forcinglike_to_DataArray, forcinglike_to_DataFrame


#NOTE: for examples with writing bc files, check dfm_tools.interpolate_grid2bnd.* and dfm_tools.hydrolib_helpers.

nPoints = 5 #None for all points

file_bc_list = [Path(r'n:\My Documents\werkmap\hydrolib_test\DCSM\tide_OB_all_20181108.bc'), #TODO: make better DataArray format for astronomic components
                #Path(r'p:\11208053-004-kpp2022-rmm1d2d\C_Work\09_Validatie2018_2020\dflowfm2d-rmm_vzm-j19_6-v2d\boundary_conditions\2020\flow\rmm_zeerand_v3_2020.bc'), #>100 timeseries
                #Path(r'p:\11208053-004-kpp2022-rmm1d2d\C_Work\09_Validatie2018_2020\dflowfm2d-rmm_vzm-j19_6-v2d\boundary_conditions\rmm_rivdis_meas_20171101_20210102_MET.bc'), #TODO: why can it not be str? #three timeseries
                Path(r'p:\11208053-004-kpp2022-rmm1d2d\C_Work\09_Validatie2018_2020\dflowfm2d-rmm_vzm-j19_6-v2d\boundary_conditions\2018\flow\rmm_discharge_laterals_20171201_20190101_MET.bc'),
                Path(r'n:\My Documents\werkmap\hydrolib_test\haixia\salinity_bc_South_v2_firstpoint.bc'), #TODO SOLVED: old keywords not supported yet
                #Path(r'n:\My Documents\werkmap\hydrolib_test\haixia\uxuy_bc_South_v2_firstpoint.bc') ,#TODO: uxuy still crashes
                #Path(r'n:\My Documents\werkmap\hydrolib_test\haixia\salinity_bc_South_v2.bc'), #large file, takes time
                #Path(r'n:\My Documents\werkmap\hydrolib_test\haixia\uxuy_bc_South_v2.bc'),
                ]


for file_bc in file_bc_list:
    #Load .bc-file using HydroLib object ForcingModel.
    m = ForcingModel(file_bc) #TODO SOLVED: in case of 3D this crashes on validation error since [verticalPositions,verticalInterpolation,verticalPositionType] are not accepted as missing (while dflowfm has default values for them) >> verticalPositions are in file as "Vertical position". Issue created: https://github.com/Deltares/HYDROLIB-core/issues/306
    
    #m.general.comments = {'a':'aa'} #TODO: adding comments to top of file is not possible, only if using filetype or fileversion: https://github.com/Deltares/HYDROLIB-core/issues/130. Top file comment newfeature: https://github.com/Deltares/HYDROLIB-core/issues/362
    #m.save('test.bc')
    
    #print(m.forcing) #TODO: spyder crasht bij printen van m.forcing, 25% CPU en veel werkgeheugen, maar vooral een freeze. pretty print suggestion: https://github.com/Deltares/HYDROLIB-core/issues/315
    """
    type(m) #hydrolib.core.io.bc.models.ForcingModel
    type(m.forcing) #list of timeseries from bc file
    type(m.forcing[0]) # hydrolib.core.io.bc.models.TimeSeries
    type(m.forcing[0].__dict__) #dict
    m.forcing[0].__dict__.keys() # dict_keys(['comments', 'datablock', 'name', 'function', 'quantityunitpair', 'timeinterpolation', 'offset', 'factor'])
    type(m.forcing[0].datablock) #list
    """
    
    pli_quan, pli_unit = m.forcing[0].quantityunitpair[1].quantity, m.forcing[0].quantityunitpair[1].unit

    #plot
    fig, ax = plt.subplots(figsize=(12, 6))
    for iFO, forcingobj in enumerate(m.forcing[:nPoints]):
        forcing_xr = forcinglike_to_DataArray(forcingobj)
        if forcingobj.function=='t3d':
            forcing_xr.T.plot() #TODO: this overwites previous plot, so does not make sense
            #ax.set_ylim(-500,5)
        elif forcingobj.function=='astronomic':
            if iFO==0:
                ax2 = ax.twinx()
                continue #skip first (invalid) point
            forcing_xr.isel(quantity=0).plot(ax=ax, label='amplitude', linewidth=0.7)
            forcing_xr.isel(quantity=1).plot(ax=ax2, label='phase', linewidth=0.7)
            ax.legend(loc=1)
            # forcing_pd = forcinglike_to_DataFrame(forcingobj) #only relevant for Astronomic, for TimeSeries/T3D it is equal to xr.to_pandas()
            # forcing_pd[forcing_pd.columns[0]].plot(ax=ax, linewidth=0.7)
            # forcing_pd[forcing_pd.columns[1]].plot(ax=ax, linewidth=0.7, secondary_y=True)
            # ax.legend(forcing_pd.columns,loc=1)
        else:
            forcing_xr.plot(ax=ax, label=forcing_xr.attrs['name'], linewidth=0.7)
            ax.legend(loc=1)
        
    """
    #Plot mean value in BC-file over the complete period for each point on the boundary
    mean_list = [forcinglike_to_DataArray(forcingobj).mean().mean() for forcingobj in m.forcing[:nPoints]]
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.plot(mean_list,'-o')
    ax.set_xlabel('Point on boundary')
    ax.set_ylabel(f"mean {pli_quan} [{pli_unit}]")
    """



