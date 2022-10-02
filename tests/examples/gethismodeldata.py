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
from dfm_tools.get_nc_helpers import get_stationid_fromstationlist#, get_hisstationlist

dir_testinput = r'c:\DATA\dfm_tools_testdata'
dir_output = '.'

file_nc_list = [os.path.join(dir_testinput,'vanNithin','tttz_0000_his.nc'),
                os.path.join(dir_testinput,'DFM_3D_z_Grevelingen\\computations\\run01\\DFM_OUTPUT_Grevelingen-FM\\Grevelingen-FM_0000_his.nc'),
                r'p:\11202512-h2020_impaqt\07_Mediterranean_model\MedSea_impaqt_model\computations_final\r013_waq\DFM_OUTPUT_MedSea_impaqt_FM\MedSea_impaqt_FM_0000_his.nc',
                ]

for file_nc in file_nc_list:
    if 'Grevelingen-FM_0000' in file_nc:
        #file_nc = os.path.join(dir_testinput,r'DFM_3D_z_Grevelingen\computations\run01\DFM_OUTPUT_Grevelingen-FM\Grevelingen-FM_0000_his.nc')
        station = ['GTSO-01','GTSO-02','GTSO-03','GTSO-04','GTSO-05','GTSO-06','GTSO-07',
                   'GTSO-08','GTSO-09','GTSO-10','GTSO-11','GTSO-12','GTSO-13','GTSO-14',
                   'GTSO-15','GTSO-16','GTSO-17','GTSO-18','GTSO-19','GTSO-20',
                   'Bommenede','Grevelingen hevel West','Brouwerssluis binnen','Brouwerssluis binnen-hand']
        station_zt = ['GTSO-02']
    elif 'tttz' in file_nc: #NITHIN
        #file_nc = os.path.join(dir_testinput,'vanNithin','tttz_0000_his.nc')
        station = ['Peiraias', 'Ovrios_2','Ovrios','Ovrios','Ortholithi']
        station_zt = ['Ortholithi']
    elif 'impaqt' in file_nc:
        station = ['MO_TS_MO_ATHOS','MO_TS_MO_LESVO','MO_TS_MO_SKYRO','IOC_thes','farm_impaqt']
        station_zt = ['MO_TS_MO_ATHOS']
    
    data_xr = xr.open_dataset(file_nc) #TODO: maybe adding chunking argument like chunks={'time':-1,'station':200}) (https://github.com/pydata/xarray/discussions/6458)
    #stations_pd = get_hisstationlist(file_nc)
    idx_stations = get_stationid_fromstationlist(data_xr, stationlist=station)
    idx_stations_zt = get_stationid_fromstationlist(data_xr, stationlist=station_zt)[0] #if provide single station (string, no list), the shape of the resulting xarray is correct
    
    print('plot bedlevel from his')
    #data_fromhis = get_ncmodeldata(file_nc=file_nc, varname='bedlevel', station=station)#, multipart=False)
    data_fromhis_xr = data_xr.bedlevel.isel(stations=idx_stations) #TODO: also possible to index directly with station strings?
    fig, ax = plt.subplots()
    #ax.plot(data_fromhis.var_stations.iloc[:,0],data_fromhis,'-')
    ax.plot(data_fromhis_xr.station_name,data_fromhis_xr,'-')
    ax.tick_params('x',rotation=90)
    fig.savefig(os.path.join(dir_output,'%s_bedlevel'%(os.path.basename(file_nc).replace('.',''))))
    
    print('plot waterlevel from his')
    #data_fromhis = get_ncmodeldata(file_nc=file_nc, varname='waterlevel', timestep='all', station=station)#, multipart=False)
    data_fromhis_xr = data_xr.waterlevel.isel(stations=idx_stations)
    fig, ax = plt.subplots()
    #ax.plot(data_fromhis.var_times,data_fromhis,'-')
    ax.plot(data_fromhis_xr.time,data_fromhis_xr.to_numpy(),'-')
    fig.savefig(os.path.join(dir_output,'%s_waterlevel'%(os.path.basename(file_nc).replace('.',''))))
    
    print('plot salinity from his')
    #data_fromhis = get_ncmodeldata(file_nc=file_nc, varname='salinity', timestep='all', layer=5, station=station)#, multipart=False)
    #data_fromhis_flat = data_fromhis[:,:,0]
    data_fromhis_xr = data_xr.salinity.isel(stations=idx_stations,laydim=5)
    fig, ax = plt.subplots()
    #ax.plot(data_fromhis.var_times,data_fromhis_flat,'-')
    ax.plot(data_fromhis_xr.time,data_fromhis_xr,'-')
    fig.savefig(os.path.join(dir_output,'%s_salinity'%(os.path.basename(file_nc).replace('.',''))))
    
    print('plot salinity over depth')
    #depth retrieval is probably wrong
    #data_fromhis_depth = get_ncmodeldata(file_nc=file_nc, varname='zcoordinate_c', timestep=4, layer='all', station=station)
    #data_fromhis = get_ncmodeldata(file_nc=file_nc, varname='salinity', timestep=4, layer='all', station=station)
    data_fromhis_depth_xr = data_xr.zcoordinate_c.isel(stations=idx_stations,time=4)
    data_fromhis_xr = data_xr.salinity.isel(stations=idx_stations,time=4)
    fig, ax = plt.subplots()
    #ax.plot(data_fromhis[0,:,:].T, data_fromhis_depth[0,:,:].T,'-')
    ax.plot(data_fromhis_xr.T, data_fromhis_depth_xr.T,'-')
    ax.legend(data_fromhis_xr.station_name.astype(str).to_numpy()) #TODO: maybe less complex via stations_pd? maybe decode upon open_dataset
    fig.savefig(os.path.join(dir_output,'%s_salinityoverdepth'%(os.path.basename(file_nc).replace('.',''))))
    
    print('zt temperature plot and wl')
    data_xr_selzt = data_xr.isel(stations=idx_stations_zt,time=slice(40,100))
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

