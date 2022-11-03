# -*- coding: utf-8 -*-
"""
Created on Wed Oct 27 14:45:30 2021

@author: veenstra
"""

import os
import xarray as xr
import matplotlib.pyplot as plt
plt.close('all')

from dfm_tools.get_nc import plot_ztdata
from dfm_tools.get_nc_helpers import get_hisstationlist
from dfm_tools.xarray_helpers import preprocess_hisnc, Dataset_varswithdim

dir_testinput = r'c:\DATA\dfm_tools_testdata'
dir_output = '.'

file_nc_list = [os.path.join(dir_testinput,'vanNithin','tttz_0000_his.nc'),
                os.path.join(dir_testinput,'DFM_3D_z_Grevelingen\\computations\\run01\\DFM_OUTPUT_Grevelingen-FM\\Grevelingen-FM_0000_his.nc'),
                r'p:\11202512-h2020_impaqt\07_Mediterranean_model\MedSea_impaqt_model\computations_final\r013_waq\DFM_OUTPUT_MedSea_impaqt_FM\MedSea_impaqt_FM_0000_his.nc',
                r'p:\11206813-006-kpp2021_rmm-2d\C_Work\31_RMM_FMmodel\computations\model_setup\run_206\results\RMM_dflowfm_0000_his.nc', #added since there are duplicate stations which are dropped
                ]

for file_nc in file_nc_list:
    if 'Grevelingen-FM_0000' in file_nc:
        #file_nc = os.path.join(dir_testinput,r'DFM_3D_z_Grevelingen\computations\run01\DFM_OUTPUT_Grevelingen-FM\Grevelingen-FM_0000_his.nc')
        stations_requested = ['GTSO-01','GTSO-02','GTSO-03','GTSO-04','GTSO-05','GTSO-06','GTSO-07',
                              'GTSO-08','GTSO-09','GTSO-10','GTSO-11','GTSO-12','GTSO-13','GTSO-14',
                              'GTSO-15','GTSO-16','GTSO-17','GTSO-18','GTSO-19','GTSO-20',
                              'Bommenede','Grevelingen hevel West','Brouwerssluis binnen','Brouwerssluis binnen-hand']
        stations_requested_zt = ['GTSO-02']
    elif 'tttz' in file_nc: #NITHIN
        #file_nc = os.path.join(dir_testinput,'vanNithin','tttz_0000_his.nc')
        stations_requested = ['Peiraias', 'Ovrios_2','Ovrios','Ovrios','Ortholithi']
        stations_requested_zt = ['Ortholithi']
    elif 'impaqt' in file_nc:
        stations_requested = ['MO_TS_MO_ATHOS','MO_TS_MO_LESVO','MO_TS_MO_SKYRO','IOC_thes','farm_impaqt']
        stations_requested_zt = ['MO_TS_MO_ATHOS']
    elif 'RMM_dflowfm' in file_nc:
        stations_requested = ['WAQ_Vuren','NW_1030.19_R_LMW-H_Hoek-van-Holland','WAQ_TielWaal_waq']
    
    data_xr = xr.open_mfdataset(file_nc, preprocess=preprocess_hisnc) #TODO: maybe adding chunking argument like chunks={'time':-1,'station':200}) (https://github.com/pydata/xarray/discussions/6458)
    #data_xr_indexlist = list(data_xr.indexes.keys()) #TODO: also add waterbalance as index?
    #data_xr_perdim = {dimname: Dataset_varswithdim(data_xr,dimname=dimname) for dimname in data_xr.dims}
    statlist_pd = data_xr['stations'].to_dataframe() #alternatively use .to_series() for labels only, also possible for other indexed dimensions like cross_section and general_structures (list(data_xr.indexes.keys()))

    print('plot waterlevel from his')
    data_fromhis_xr = data_xr.waterlevel.sel(stations=stations_requested)
    fig, ax = plt.subplots(figsize=(10,6))
    data_fromhis_xr.plot.line('-',ax=ax,x='time')
    ax.legend(data_fromhis_xr.stations.to_series(),fontsize=9) #optional, to reduce legend font size
    data_fromhis_xr_dailymean = data_fromhis_xr.resample(time='D').mean(dim='time') #add daily mean values in the back
    data_fromhis_xr_dailymean.plot.line('-',ax=ax,x='time',add_legend=False,zorder=0,linewidth=.8,color='grey')
    fig.tight_layout()
    fig.savefig(os.path.join(dir_output,'%s_waterlevel'%(os.path.basename(file_nc).replace('.',''))))
    if 'RMM_dflowfm' in file_nc:
        continue

    print('plot bedlevel from his')
    data_fromhis_xr = data_xr.bedlevel.sel(stations=stations_requested)
    fig, ax = plt.subplots(figsize=(10,6))
    data_fromhis_xr.plot.line('-',ax=ax)
    ax.set_xticklabels(data_fromhis_xr.stations.to_series(),rotation=45,ha='right') #optional, to rotate x-labels
    fig.tight_layout()
    fig.savefig(os.path.join(dir_output,'%s_bedlevel'%(os.path.basename(file_nc).replace('.',''))))
    
    print('plot salinity from his')
    data_fromhis_xr = data_xr.salinity.sel(stations=stations_requested).isel(laydim=20)
    fig, ax = plt.subplots(figsize=(10,6))
    data_fromhis_xr.plot.line('-',ax=ax,x='time')
    ax.legend(data_fromhis_xr.stations.to_series(),fontsize=9) #optional, to reduce legend font size
    fig.tight_layout()
    fig.savefig(os.path.join(dir_output,'%s_salinity'%(os.path.basename(file_nc).replace('.',''))))
    
    print('plot salinity over depth')
    #depth retrieval is probably wrong
    data_fromhis_xr = data_xr.salinity.sel(stations=stations_requested).isel(time=4)
    fig, ax = plt.subplots(figsize=(10,6))
    data_fromhis_xr.T.plot.line('-',ax=ax,y='zcoordinate_c')
    ax.legend(data_fromhis_xr.stations.to_series(),fontsize=9) #optional, to reduce legend font size
    fig.tight_layout()
    fig.savefig(os.path.join(dir_output,'%s_salinityoverdepth'%(os.path.basename(file_nc).replace('.',''))))
    
    print('zt temperature plot and wl')
    station_zt = stations_requested_zt[0]
    data_xr_selzt = data_xr.sel(stations=station_zt).isel(time=slice(40,100))
    data_fromhis_wl_xr = data_xr_selzt['waterlevel']
    fig, (axwl,ax1) = plt.subplots(2,1,figsize=(12,7),gridspec_kw={'height_ratios':[1,2]},sharex=True,sharey=True)
    axwl.plot(data_xr_selzt.time[[0,-1]],[0,0],'k-',linewidth=0.5)
    ax1.plot(data_xr_selzt.time[[0,-1]],[0,0],'k-',linewidth=0.5)
    axwl.plot(data_xr_selzt.time,data_fromhis_wl_xr,'-',label=f'wl {station_zt}')
    c = plot_ztdata(data_xr_sel=data_xr_selzt, varname='temperature', ax=ax1, cmap='jet')
    fig.colorbar(c,ax=axwl)
    fig.colorbar(c,ax=ax1)
    #contour
    CS = plot_ztdata(data_xr_sel=data_xr_selzt, varname='temperature', ax=ax1, only_contour=True, levels=6, colors='k', linewidths=0.8, linestyles='solid')
    ax1.clabel(CS, fontsize=10)
    fig.tight_layout()
    fig.savefig(os.path.join(dir_output,'%s_zt_temp'%(os.path.basename(file_nc).replace('.',''))))
    axwl.set_ylim(-2,0.5)
    fig.savefig(os.path.join(dir_output,'%s_zt_temp_zoomwl'%(os.path.basename(file_nc).replace('.',''))))

