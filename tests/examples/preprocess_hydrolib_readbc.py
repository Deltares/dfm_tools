# -*- coding: utf-8 -*-
"""
Created on Wed Aug 17 11:19:51 2022

@author: veenstra
"""

import os
import matplotlib.pyplot as plt
plt.close('all')
import dfm_tools as dfmt
import hydrolib.core.dflowfm as hcdfm


#TODO: merge this into preprocess_hydrolib_readFMmodel.py after issues are resolved?
#NOTE: for examples with writing bc files, check dfm_tools.interpolate_grid2bnd.* and dfm_tools.hydrolib_helpers.

file_bc_list = [r'p:\11208154-002-haixia\02-hydrodynamics\02_Model_set_up\02_Make_forcing\FES2014\New_bnd_lines_2022\bc_South_v2.bc',
                r'p:\11208053-004-kpp2022-rmm1d2d\C_Work\09_Validatie2018_2020\dflowfm2d-rmm_vzm-j19_6-v2d\boundary_conditions\rmm_rivdis_meas_20171101_20210102_MET.bc', #three timeseries
                r'p:\11208154-002-haixia\02-hydrodynamics\02_Model_set_up\02_Make_forcing\CMEMS\bc_2022\South_v2\nonan\salinitybnd_bc_South_v2_CMEMS.bc',
                r'p:\11208154-002-haixia\02-hydrodynamics\02_Model_set_up\02_Make_forcing\CMEMS\bc_2022\South_v2\nonan\uxuy_bc_South_v2_CMEMS.bc',
                ]

dir_output = '.'

for file_bc in file_bc_list:
    #Load .bc-file using HydroLib object ForcingModel.
    m = hcdfm.ForcingModel(file_bc)
    ForcingModel_object_out = hcdfm.ForcingModel()
    
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
    forcingobj = m.forcing[1]
    forcing_xr = dfmt.forcinglike_to_Dataset(forcingobj, convertnan=True)
    data_vars = list(forcing_xr.data_vars)
    if forcingobj.function=='t3d':
        forcing_ts = dfmt.Dataset_to_T3D(forcing_xr)
        if hasattr(forcingobj.quantityunitpair[1],'elementname'): #uxuy vector
            plt.close()
            fig, axes = plt.subplots(2,1,figsize=(12, 8),sharex=True,sharey=True)
            forcing_xr[data_vars[0]].T.plot(ax=axes[0])
            forcing_xr[data_vars[1]].T.plot(ax=axes[1])
        else: #salinitybnd/temperaturebnd
            forcing_xr[data_vars[0]].T.plot(ax=ax)
    elif forcingobj.function=='timeseries': #waterlevelbnd
        forcing_ts = dfmt.Dataset_to_TimeSeries(forcing_xr)
        forcing_xr[data_vars[0]].plot(ax=ax, label=forcing_xr[data_vars[0]].attrs['name'], linewidth=0.8)
        ax.legend(loc=1)
    elif forcingobj.function=='astronomic': #eg tidal components
        #forcing_ts = Dataset_to_Astronomic(forcing_xr) #TODO: implement Dataset_to_Astronomic() function
        ax2 = ax.twinx()
        forcing_xr[data_vars[0]].plot(ax=ax, linewidth=0.8)
        forcing_xr[data_vars[1]].plot(ax=ax2, linewidth=0.8)
    else:
        forcing_xr.plot(ax=ax, label=forcing_xr.attrs['name'], linewidth=0.8)
        ax.legend(loc=1)
        raise Exception(f'non-defined function: {forcingobj.function}')
    #ForcingModel_object_out.forcing.append(forcing_ts)
    #ForcingModel_object_out.save(Path(file_bc.replace('.bc','_reproduced.bc')))
    fig.tight_layout()
    fig.savefig(os.path.join(dir_output,os.path.basename(file_bc.replace('.',''))))


