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

dir_testinput = r'c:\DATA\dfm_tools_testdata'
dir_output = '.'

file_nc_list = [os.path.join(dir_testinput,'vanNithin','tttz_0000_his.nc'),
                os.path.join(dir_testinput,'DFM_3D_z_Grevelingen\\computations\\run01\\DFM_OUTPUT_Grevelingen-FM\\Grevelingen-FM_0000_his.nc'),
                r'p:\archivedprojects\11202512-h2020_impaqt\07_Mediterranean_model\MedSea_impaqt_model\computations_final\r013_waq\DFM_OUTPUT_MedSea_impaqt_FM\MedSea_impaqt_FM_0000_his.nc',
                r'p:\11206813-006-kpp2021_rmm-2d\C_Work\31_RMM_FMmodel\computations\model_setup\run_206\results\RMM_dflowfm_0000_his.nc', #contains duplicate station_names which are dropped
                r'p:\archivedprojects\11203869-morwaqeco3d\05-Tidal_inlet\02_FM_201910\FM_MF10_Max_30s\fm\DFM_OUTPUT_inlet\inlet_his.nc', #morphology
                r'p:\11202255-sfincs\course_material\DSD_INT_2022\02_hands-on\Charleston_subgrid_example_allforcing\pre-run_output\sfincs_his.nc', #SFINCS
                os.path.join(dir_testinput,'D3D_3D_sigma_curved_bend_nc\\trih-cb2-sal-added-3d.nc'), #DELFT3D4 netcdf
                r'p:\archivedprojects\1220688-lake-kivu\3_modelling\1_FLOW\7_heatfluxinhis\063_netcdf\trih-thiery_002_coarse.nc', #DELFT3D4 netcdf
                ]


for file_nc in file_nc_list:
    basename = os.path.basename(file_nc).replace('.','')
    
    data_xr = xr.open_mfdataset(file_nc, preprocess=dfmt.preprocess_hisnc) #TODO: maybe adding chunking argument like chunks={'time':-1,'station':200}) (https://github.com/pydata/xarray/discussions/6458)
    vars_pd = dfmt.get_ncvarproperties(data_xr)
    
    if 'Grevelingen-FM_0000' in file_nc:
        #file_nc = os.path.join(dir_testinput,r'DFM_3D_z_Grevelingen\computations\run01\DFM_OUTPUT_Grevelingen-FM\Grevelingen-FM_0000_his.nc')
        stations_requested = ['GTSO-01','GTSO-02','GTSO-03','GTSO-04','GTSO-05','GTSO-06','GTSO-07',
                              'GTSO-08','GTSO-09','GTSO-10','GTSO-11','GTSO-12','GTSO-13','GTSO-14',
                              'GTSO-15','GTSO-16','GTSO-17','GTSO-18','GTSO-19','GTSO-20',
                              'Bommenede','Grevelingen hevel West','Brouwerssluis binnen','Brouwerssluis binnen-hand']
    elif 'tttz' in file_nc: #NITHIN
        #file_nc = os.path.join(dir_testinput,'vanNithin','tttz_0000_his.nc')
        stations_requested = ['Peiraias', 'Ovrios_2','Ovrios','Ovrios','Ortholithi']
    elif 'impaqt' in file_nc:
        stations_requested = ['MO_TS_MO_ATHOS','MO_TS_MO_LESVO','MO_TS_MO_SKYRO','IOC_thes','farm_impaqt']
    elif 'RMM_dflowfm' in file_nc:
        stations_requested = ['WAQ_Vuren','NW_1030.19_R_LMW-H_Hoek-van-Holland','WAQ_TielWaal_waq']
    elif 'moergestels_broek' in file_nc:
        stations_requested = ['ObsPt1', 'ObsPt2', 'ObsPt2D1', 'ObsPt2D2']
    elif 'morwaqeco3d' in file_nc:
        stations_requested = ['NW', 'N', 'NE', 'W', 'Centre', 'E', 'SW', 'SE', 'Inlet_out_NW',
                              'Inlet_out_NE', 'Seegat', 'Inlet_in_1', 'Inlet_in_2', 'Inlet_in_3', 'Inlet_in_4', 'Inlet_in_5']
    elif 'sfincs' in file_nc:
        stations_requested = ['charleston_mouth', 'cooper_river', 'ashley_river']
        data_xr = data_xr.rename({'point_zs':'waterlevel','point_zb':'bedlevel'}) # for convenience in this script
    elif 'trih-cb2' in file_nc:
        stations_requested = ['Outer-south', 'inner-south', 'inner-middle']
        data_xr = data_xr.rename({'NOSTAT':'stations','ZWL':'waterlevel','DPS':'bedlevel'}) # for convenience in this script
    elif 'trih-thiery_002_coarse' in file_nc:
        stations_requested = ['ADCP1_final','ADCP2_final','KP1_016']
        data_xr = data_xr.rename({'NOSTAT':'stations','ZWL':'waterlevel','DPS':'bedlevel'}) # for convenience in this script
    
    #data_xr_indexlist = list(data_xr.indexes.keys()) #TODO: also add waterbalance as index?
    #data_xr_perdim = {dimname: dfmt.Dataset_varswithdim(data_xr,dimname=dimname) for dimname in data_xr.dims}
    statlist_pd = data_xr['stations'].to_dataframe() #alternatively use .to_series() for labels only, also possible for other indexed dimensions like cross_section and general_structures (list(data_xr.indexes.keys()))
    
    print('plot waterlevel from his')
    data_fromhis_xr = data_xr.waterlevel.sel(stations=stations_requested)
    fig, ax = plt.subplots(figsize=(10,6))
    data_fromhis_xr.plot.line('-',ax=ax,x='time')
    ax.legend(data_fromhis_xr.stations.to_series(),fontsize=9) #optional, to reduce legend font size
    data_fromhis_xr_dailymean = data_fromhis_xr.resample(time='D').mean(dim='time') #add daily mean values in the back
    data_fromhis_xr_dailymean.plot.line('-',ax=ax,x='time',add_legend=False,zorder=0,linewidth=.8,color='grey')
    fig.tight_layout()
    fig.savefig(os.path.join(dir_output,f'{basename}_waterlevel'))
    if 'RMM_dflowfm' in file_nc or 'moergestels_broek' in file_nc:
        continue
    
    print('plot bedlevel from his')
    data_fromhis_xr = data_xr.bedlevel.sel(stations=stations_requested)
    fig, ax = plt.subplots(figsize=(10,6))
    if 'time' in data_fromhis_xr.dims:
        data_fromhis_xr.plot.line('-',ax=ax,x='time')
    else:
        data_fromhis_xr.plot.line('-',ax=ax)
        stat_pd = data_fromhis_xr.stations.to_series()
        ax.set_xticks(range(len(stat_pd)))
        ax.set_xticklabels(stat_pd,rotation=45,ha='right') #optional, to rotate x-labels
    fig.tight_layout()
    fig.savefig(os.path.join(dir_output,f'{basename}_bedlevel'))
    if 'sfincs' in file_nc or 'morwaqeco3d' in file_nc or 'trih-' in file_nc:
        continue
    
    #at this stage there are only three hisfiles left where continue was not called for
    
    print('plot salinity from his') #TODO: retreive variable on fixed depth (like dfmt.get_mapdata_atfixedepth())
    data_fromhis_xr = data_xr.salinity.sel(stations=stations_requested).isel(laydim=20)
    fig, ax = plt.subplots(figsize=(10,6))
    data_fromhis_xr.plot.line('-',ax=ax,x='time')
    ax.legend(data_fromhis_xr.stations.to_series(),fontsize=9) #optional, to reduce legend font size
    fig.tight_layout()
    fig.savefig(os.path.join(dir_output,f'{basename}_salinity'))
    
    print('plot salinity over depth')
    #depth retrieval is probably wrong
    data_fromhis_xr = data_xr.salinity.sel(stations=stations_requested).isel(time=4)
    fig, ax = plt.subplots(figsize=(10,6))
    data_fromhis_xr.T.plot.line('-',ax=ax,y='zcoordinate_c')
    ax.legend(data_fromhis_xr.stations.to_series(),fontsize=9) #optional, to reduce legend font size
    fig.tight_layout()
    fig.savefig(os.path.join(dir_output,f'{basename}_salinityoverdepth'))
    
    print('zt temperature plot and wl')
    data_xr_selzt = data_xr.isel(stations=2).isel(time=slice(40,100))
    fig, ax = plt.subplots(1,1,figsize=(12,7))
    data_xr_selzt.waterlevel.plot.line(ax=ax,color='r') #waterlevel line
    pc = dfmt.plot_ztdata(data_xr_sel=data_xr_selzt, varname='temperature', ax=ax, cmap='jet') #temperature pcolormesh
    fig.colorbar(pc,ax=ax)
    CS = dfmt.plot_ztdata(data_xr_sel=data_xr_selzt, varname='temperature', ax=ax, only_contour=True, levels=6, colors='k', linewidths=0.8, linestyles='solid') #temperature contour
    ax.clabel(CS, fontsize=10)
    fig.tight_layout()
    fig.savefig(os.path.join(dir_output,f'{basename}_zt_temp'))
    ax.set_ylim(-2,0.5)
    fig.savefig(os.path.join(dir_output,f'{basename}_zt_temp_zoomwl'))

