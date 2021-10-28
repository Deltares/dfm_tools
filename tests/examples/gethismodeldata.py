# -*- coding: utf-8 -*-
"""
Created on Wed Oct 27 14:45:30 2021

@author: veenstra
"""

import os
import matplotlib.pyplot as plt
plt.close('all')

from dfm_tools.get_nc import get_ncmodeldata, plot_ztdata

dir_testinput = r'c:\DATA\dfm_tools_testdata'
dir_output = '.'

file_nc_list = [os.path.join(dir_testinput,'vanNithin','tttz_0000_his.nc'),
                os.path.join(dir_testinput,'DFM_3D_z_Grevelingen\\computations\\run01\\DFM_OUTPUT_Grevelingen-FM\\Grevelingen-FM_0000_his.nc'),
                #'p:\\11202512-h2020_impaqt\\Mediterranean_model\\MedSea_impaqt_model\\computations\\r003_test\\DFM_OUTPUT_MedSea_impaqt_FM\\MedSea_impaqt_FM_0000_his.nc',
                ]

for file_nc in file_nc_list:
    def cen2cor(time_cen):
        #convert time centers2corners (more accurate representation in zt-plot, but can also be skipped)
        import pandas as pd
        time_int = (time_cen.iloc[2]-time_cen.iloc[1])
        time_cor = data_fromhis_temp.var_times-time_int/2
        time_cor = time_cor.append(pd.Series([time_cor.iloc[-1]+time_int]))
        return time_cor

    if 'Grevelingen-FM_0000' in file_nc:
        #file_nc = os.path.join(dir_testinput,r'DFM_3D_z_Grevelingen\computations\run01\DFM_OUTPUT_Grevelingen-FM\Grevelingen-FM_0000_his.nc')
        station = 'all'
        station_zt = 'GTSO-02'
    elif 'tttz' in file_nc: #NITHIN
        #file_nc = os.path.join(dir_testinput,'vanNithin','tttz_0000_his.nc')
        station = ['Peiraias', 'Ovrios_2','Ovrios','Ovrios','Ortholithi']
        station_zt = 'Ortholithi'
    elif 'impaqt' in file_nc:
        station_zt = 'IOC_thes' #['IOC_thes','farm_impaqt']

    print('plot bedlevel from his')
    data_fromhis = get_ncmodeldata(file_nc=file_nc, varname='bedlevel', station=station)#, multipart=False)
    fig, ax = plt.subplots()
    ax.plot(data_fromhis.var_stations.iloc[:,0],data_fromhis,'-')
    ax.tick_params('x',rotation=90)
    fig.savefig(os.path.join(dir_output,'%s_bedlevel'%(os.path.basename(file_nc).replace('.',''))))

    print('plot waterlevel from his')
    data_fromhis = get_ncmodeldata(file_nc=file_nc, varname='waterlevel', timestep='all', station=station)#, multipart=False)
    fig, ax = plt.subplots()
    ax.plot(data_fromhis.var_times,data_fromhis,'-')
    fig.savefig(os.path.join(dir_output,'%s_waterlevel'%(os.path.basename(file_nc).replace('.',''))))
    
    print('plot salinity from his')
    data_fromhis = get_ncmodeldata(file_nc=file_nc, varname='salinity', timestep='all', layer=5, station=station)#, multipart=False)
    data_fromhis_flat = data_fromhis[:,:,0]
    fig, ax = plt.subplots()
    ax.plot(data_fromhis.var_times,data_fromhis_flat,'-')
    fig.savefig(os.path.join(dir_output,'%s_salinity'%(os.path.basename(file_nc).replace('.',''))))

    print('plot salinity over depth')
    #depth retrieval is probably wrong
    data_fromhis_depth = get_ncmodeldata(file_nc=file_nc, varname='zcoordinate_c', timestep=4, layer='all', station=station)
    data_fromhis = get_ncmodeldata(file_nc=file_nc, varname='salinity', timestep=4, layer='all', station=station)
    fig, ax = plt.subplots()
    ax.plot(data_fromhis[0,:,:].T, data_fromhis_depth[0,:,:].T,'-')
    ax.legend(data_fromhis.var_stations.iloc[:,0])
    fig.savefig(os.path.join(dir_output,'%s_salinityoverdepth'%(os.path.basename(file_nc).replace('.',''))))

    print('zt temperature plot and wl')
    data_fromhis_temp = get_ncmodeldata(file_nc=file_nc, varname='temperature', timestep=range(40,100), layer= 'all', station=station_zt)
    data_fromhis_wl = get_ncmodeldata(file_nc=file_nc, varname='waterlevel', timestep=range(40,100), station=station_zt)#, multipart=False)
    fig, (axwl,ax1) = plt.subplots(2,1,figsize=(12,7),gridspec_kw={'height_ratios':[1,2]},sharex=True,sharey=True)
    statid_subset = data_fromhis_wl.var_stations['station_id'].tolist().index(station_zt)
    axwl.plot(data_fromhis_wl.var_times.iloc[[0,-1]],[0,0],'k-',linewidth=0.5)
    ax1.plot(data_fromhis_wl.var_times.iloc[[0,-1]],[0,0],'k-',linewidth=0.5)
    axwl.plot(data_fromhis_wl.var_times,data_fromhis_wl[:,statid_subset],'-',label='wl %s'%(station_zt))
    c = plot_ztdata(file_nc=file_nc, dfmtools_hisvar=data_fromhis_temp, ax=ax1, statid_subset=statid_subset, cmap='jet')
    fig.colorbar(c,ax=axwl)
    fig.colorbar(c,ax=ax1)
    #contour
    CS = plot_ztdata(file_nc=file_nc, dfmtools_hisvar=data_fromhis_temp, ax=ax1, statid_subset=statid_subset, only_contour=True, levels=6, colors='k', linewidths=0.8, linestyles='solid')
    ax1.clabel(CS, fontsize=10)
    fig.tight_layout()
    fig.savefig(os.path.join(dir_output,'%s_zt_temp'%(os.path.basename(file_nc).replace('.',''))))
    axwl.set_ylim(-2,0.5)
    fig.savefig(os.path.join(dir_output,'%s_zt_temp_zoomwl'%(os.path.basename(file_nc).replace('.',''))))

    
    

