# -*- coding: utf-8 -*-
"""
Created on Wed Oct 27 14:45:30 2021

@author: veenstra
"""

import os
import xarray as xr
import matplotlib.pyplot as plt
plt.close('all')
import dfm_tools as dfmt
import hydrolib.core.dflowfm as hcdfm


plotting = True

# Discharge locations 
crs_q = ['WL_934.8_QL_Zaltbommel', 'MA_230.61_QL_Heesbeen']
bc_q = ['Waal_0001','Maas_0001']
map_q = dict(zip(crs_q, bc_q))
# Water level locations
station_wl = ['BE_976.00', 'HD_983.41_R_LMW-Cl_Moerdijkbrug'] 
bc_wl = ['Beneden-Merwede_0001','Nieuwe-Merwede_0001']
map_wl = dict(zip(station_wl, bc_wl))
dir_input = r'p:\archivedprojects\11206813-006-kpp2021_rmm-2d\C_Work\07_Baseline\baseline-rmm_vzm-beno19_6-v1\models\dflowfm\dflowfm2d-rmm_vzm-beno19_6-v1b\computations\test'
dir_output = '.'

sims = ['tba']#,'tbb','tbc','tbd','tbe','tbf']

for sim in sims: 
    file_bc_output = os.path.join(dir_output,f'{sim}_bnd.bc')
    file_nc = os.path.join(dir_input,sim, 'results','RMM_VZM_0000_his.nc')
    data_xr = xr.open_mfdataset(file_nc,preprocess=dfmt.preprocess_hisnc)
    
    if plotting:
        print('plot waterlevel from his')
        fig, ax = plt.subplots()
        data_fromhis_xr = data_xr.waterlevel.sel(stations=station_wl)
        data_fromhis_xr.plot.line('-',ax=ax,x='time')
        ax.grid(True)
        #ax.tick_params('x',rotation=20)
        fig.tight_layout()
        fig.savefig(os.path.join(dir_output,f'{sim}_waterlevel.png'))
        
        print('plot discharge from his')
        fig, ax = plt.subplots()
        data_fromhis_xr = data_xr.cross_section_discharge.sel(cross_section=crs_q)
        data_fromhis_xr.plot.line('-',ax=ax,x='time')
        ax.grid(True)
        #ax.tick_params('x',rotation=20)
        fig.tight_layout()
        fig.savefig(os.path.join(dir_output,f'{sim}_discharge.png'))
        
    ForcingModel_object = hcdfm.ForcingModel()
    for stat_wl in station_wl: 
        data_fromhis_xr = data_xr.waterlevel.sel(stations=stat_wl)
        data_fromhis_xr.attrs['locationname'] = map_wl[stat_wl]
        data_fromhis_xr.name = 'waterlevelbnd'
        t = dfmt.Dataset_to_TimeSeries(data_fromhis_xr)
        forcingobject_one_xr = dfmt.forcinglike_to_Dataset(t,convertnan=True)
        ForcingModel_object.forcing.append(t)

    for crs_q_one in crs_q: 
        data_fromhis_xr = data_xr.cross_section_discharge.sel(cross_section=crs_q_one)
        data_fromhis_xr.attrs['locationname'] = map_q[crs_q_one]
        data_fromhis_xr.attrs['units'] = 'm3/s'
        data_fromhis_xr.name = 'dischargebnd'
        t = dfmt.Dataset_to_TimeSeries(data_fromhis_xr)
        forcingobject_one_xr = dfmt.forcinglike_to_Dataset(t,convertnan=True)
        ForcingModel_object.forcing.append(t)

    ForcingModel_object.save(file_bc_output)
    
# ext_new = ExtModel()
# #generate boundary object for the ext file (quantity, pli-filename, bc-filename)
# boundary_object = Boundary(quantity=quantity,
#                             locationfile=Path(dir_output,file_pli.name),
#                             forcingfile=ForcingModel_object,
#                             )
# ext_new.boundary.append(boundary_object)