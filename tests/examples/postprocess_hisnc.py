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

file_nc_list = [dfmt.data.fm_grevelingen_his(return_filepath=True),
                # r'p:\archivedprojects\11202512-h2020_impaqt\07_Mediterranean_model\MedSea_impaqt_model\computations_final\r013_waq\DFM_OUTPUT_MedSea_impaqt_FM\MedSea_impaqt_FM_0000_his.nc',
                # r'p:\archivedprojects\11206811-002-kpp-veerse-meer\grove-model\vm_grof-j19_6-v1a\scenarios\S0\S0_run4\DFM_OUTPUT_VM_WQ_3D_grof\VM_WQ_3D_grof_0000_his.nc', #veersemeer, hisfile with proper z-coordinates
                ]


for file_nc in file_nc_list:
    # plt.close('all')
    
    print('processing %s'%(os.path.basename(file_nc)))
    basename = os.path.basename(file_nc).replace('.','')
    
    data_xr = xr.open_mfdataset(file_nc, preprocess=dfmt.preprocess_hisnc) #TODO: maybe adding chunking argument like chunks={'time':-1,'station':200}) (https://github.com/pydata/xarray/discussions/6458)
    vars_pd = dfmt.get_ncvarproperties(data_xr)
    
    if 'Grevelingen-FM_0000' in file_nc:
        stations_requested = ['GTSO-01','GTSO-02','GTSO-03','GTSO-04','GTSO-05','GTSO-06','GTSO-07',
                              'GTSO-08','GTSO-09','GTSO-10','GTSO-11','GTSO-12','GTSO-13','GTSO-14',
                              'GTSO-15','GTSO-16','GTSO-17','GTSO-18','GTSO-19','GTSO-20',
                              'Bommenede','Grevelingen hevel West','Brouwerssluis binnen','Brouwerssluis binnen-hand']
    elif 'impaqt' in file_nc:
        stations_requested = ['MO_TS_MO_ATHOS','MO_TS_MO_LESVO','MO_TS_MO_SKYRO','IOC_thes','farm_impaqt']
    elif 'VM_WQ_3D_grof' in file_nc:
        stations_requested = ['TSO_VM-3','TSO_VM-10','TSO_VM-16']
    
    print('plot waterlevel from his')
    data_fromhis_xr = data_xr.waterlevel.sel(stations=stations_requested)
    fig, ax = plt.subplots(figsize=(10,6))
    data_fromhis_xr.plot.line('-',ax=ax,x='time')
    ax.legend(data_fromhis_xr.stations.to_series(),fontsize=9) #optional, to reduce legend font size
    data_fromhis_xr_dailymean = data_fromhis_xr.resample(time='D').mean(dim='time') #add daily mean values in the back #TODO: raises "TypeError: __init__() got an unexpected keyword argument 'base'" since py39 environment
    data_fromhis_xr_dailymean.plot.line('-',ax=ax,x='time',add_legend=False,zorder=0,linewidth=.8,color='grey')
    
   
    print('plot salinity from his')
    data_fromhis_xr = data_xr.salinity.sel(stations=stations_requested).isel(laydim=20)
    fig, ax = plt.subplots(figsize=(10,6))
    data_fromhis_xr.plot.line('-',ax=ax,x='time')
    ax.legend(data_fromhis_xr.stations.to_series(),fontsize=9) #optional, to reduce legend font size
    
    
    print('zt temperature plot and wl')
    data_xr_selzt = data_xr.isel(stations=2).isel(time=slice(40,100))
    fig, ax = plt.subplots(1,1,figsize=(12,7))
    data_xr_selzt.waterlevel.plot.line(ax=ax,color='r') #waterlevel line
    pc = dfmt.plot_ztdata(data_xr_sel=data_xr_selzt, varname='temperature', ax=ax, cmap='jet') #temperature pcolormesh
    fig.colorbar(pc,ax=ax)
    CS = dfmt.plot_ztdata(data_xr_sel=data_xr_selzt, varname='temperature', ax=ax, only_contour=True, levels=6, colors='k', linewidths=0.8, linestyles='solid') #temperature contour
    ax.clabel(CS, fontsize=10)
    
    print('zt temperature plot sliced at depth(s)')
    depths = [-1,-4,0,-6]
    data_fromhis_atdepths = dfmt.get_Dataset_atdepths(data_xr=data_xr, depths=depths, reference='z0') #depth w.r.t. z0/waterlevel/bedlevel
    data_xr_selzt = data_fromhis_atdepths.isel(stations=2).isel(time=slice(40,100))
    fig, ax = plt.subplots(1,1,figsize=(12,7))
    data_xr_selzt.temperature.plot(ax=ax, cmap='jet', x='time')

