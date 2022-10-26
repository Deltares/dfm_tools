# -*- coding: utf-8 -*-
"""
Created on Wed Oct 27 14:45:30 2021

@author: veenstra
"""

from cProfile import label
import os
import xarray as xr
import matplotlib.pyplot as plt

from dfm_tools.get_nc_helpers import get_stationid_fromstationlist, get_hisstationlist
import dfm_tools.hydrolib_helpers
from dfm_tools.hydrolib_helpers import forcinglike_to_Dataset, Dataset_to_TimeSeries, Dataset_to_T3D#, Dataset_to_Astronomic
from hydrolib.core.io.ext.models import Boundary, ExtModel, ForcingModel

plotting = True

# Discharge locations 
crs_q = ['WL_934.8_QL_Zaltbommel', 'MA_230.61_QL_Heesbeen']
bc_q = ['Waal_0001','Maas_0001']
map_q = dict(zip(crs_q, bc_q))
# Water level locations 
station_wl = ['BE_976.00', 'HD_983.41_R_LMW-Cl_Moerdijkbrug'] 
bc_wl = ['Beneden-Merwede_0001','Nieuwe-Merwede_0001']
map_wl = dict(zip(station_wl, bc_wl))
#dir_input = r'p:\11206813-006-kpp2021_rmm-2d\C_Work\07_Baseline\baseline-rmm_vzm-beno19_6-v1\models\dflowfm\dflowfm2d-rmm_vzm-beno19_6-v1b\computations\test'
dir_input = '.'
dir_output = '.'

sims = ['tba','tbb','tbc','tbd','tbe','tbf']

for sim in sims: 
    file_bc_output = [os.path.join(dir_output,f'{sim}_bnd.bc')]
    file_nc = os.path.join(dir_input,sim, 'results','RMM_VZM_0000_his.nc')
    data_xr = xr.open_dataset(file_nc) #TODO: maybe adding chunking argument like chunks={'time':-1,'station':200}) (https://github.com/pydata/xarray/discussions/6458)
    data_xr['station_name_str'] = data_xr['station_name'].str.decode('utf-8',errors='ignore').str.strip() #TODO: this is currently necessary but might not. Often, .astype(str) is enough, but .decode() and .strip() is necesary for this file: p:\\11208067-003-kpp-internationaal\\final_results_Snellius\\current_2006_25per\\DFM_OUTPUT_DCSM-FM_0_5nm_waq\\DCSM-FM_0_5nm_waq_0000_his.nc
    data_xr = data_xr.set_coords('station_name_str')
    
    #all_stations_wl = get_hisstationlist(file_nc)
    station_idx = get_stationid_fromstationlist(data_xr, stationlist=station_wl)
    
    #all_crs_q = get_hisstationlist(file_nc, varname='cross_section_discharge')    
    crs_idx = get_stationid_fromstationlist(data_xr, stationlist=crs_q, station_varname='cross_section_name')

    if plotting:
        print('plot waterlevel from his')
        fig, ax = plt.subplots()
        for counter in station_idx: 
            data_fromhis_xr = data_xr.waterlevel.isel(stations=counter) #TODO: also possible to index directly with station strings?
            #
            #plt.plot(data_fromhis_xr.time,data_fromhis_xr,label=station)
            data_fromhis_xr.plot.line('-',ax=ax,x='time',label=map_wl[str(data_fromhis_xr.station_name_str.values)])
        plt.legend()
        plt.title('')
        ax.grid(True)
        ax.tick_params('x',rotation=20)
        fig.savefig(os.path.join(dir_output,'%s_waterlevel'%(os.path.basename(file_nc).replace('.',''))))
        

        print('plot discharge from his')
        #data_fromhis = get_ncmodeldata(file_nc=file_nc, varname='bedlevel', station=station)#, multipart=False)
        data_fromhis_xr = data_xr.cross_section_discharge.isel(cross_section=crs_idx) #TODO: also possible to index directly with station strings?
        fig, ax = plt.subplots()
        for counter in crs_idx: 
            data_fromhis_xr = data_xr.cross_section_discharge.isel(cross_section=counter) #TODO: also possible to index directly with station strings?
            #ax.plot(data_fromhis.var_stations.iloc[:,0],data_fromhis,'-')
            #ax.plot(data_fromhis_xr.station_name,data_fromhis_xr,'-')
            data_fromhis_xr.plot.line('-',ax=ax,x='time',label=map_q[str(data_fromhis_xr.cross_section_name.values.astype('str')).strip()])
        plt.legend()
        plt.title('')
        ax.grid(True)
        ax.tick_params('x',rotation=20)
        fig.savefig(os.path.join(dir_output,'%s_discharge'%(os.path.basename(file_nc).replace('.',''))))


    ForcingModel_object = ForcingModel()
    for counter in station_idx: 
        data_fromhis_xr = data_xr.waterlevel.isel(stations=counter) #TODO: also possible to index directly with station strings?
        data_fromhis_xr.attrs['locationname']=map_wl[str(data_fromhis_xr.station_name_str.values)]
        data_fromhis_xr.name = 'waterlevelbnd'
        t = Dataset_to_TimeSeries(data_fromhis_xr)
        forcingobject_one_xr = forcinglike_to_Dataset(t,convertnan=True)
        ForcingModel_object.forcing.append(t)

    for counter in crs_idx: 
        data_fromhis_xr = data_xr.cross_section_discharge.isel(cross_section=counter) #TODO: also possible to index directly with station strings?
        data_fromhis_xr.attrs['locationname']=map_q[str(data_fromhis_xr.cross_section_name.values.astype('str')).strip()]
        data_fromhis_xr.attrs['units'] = 'm3/s'
        data_fromhis_xr.name = 'dischargebnd'
        t = Dataset_to_TimeSeries(data_fromhis_xr)
        forcingobject_one_xr = forcinglike_to_Dataset(t,convertnan=True)
        ForcingModel_object.forcing.append(t)

    ForcingModel_object.save('test.bc')

# ext_bnd = ExtModel()
# #generate boundary object for the ext file (quantity, pli-filename, bc-filename)
# boundary_object = Boundary(quantity=quantity,
#                             locationfile=Path(dir_output,file_pli.name),
#                             forcingfile=ForcingModel_object,
#                             )
# ext_bnd.boundary.append(boundary_object)
pass