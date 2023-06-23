# -*- coding: utf-8 -*-
"""
Created on Wed Oct 27 21:56:54 2021

@author: veenstra


RMM testmodel: convert all waqua output to netcdf after completion of the run by adding this to siminp file:
    NETCDFOUTPUT
      MAPS
    	OUTPUTNAME = 'nc_map.nc'
      HISTORIES
    	OUTPUTNAME = 'nc_his.nc'
      OPTIONS
        MAPEXTRA = 'HZETA,VICO'
        HISEXTRA = 'Z0'
    
    waqpro.pl uses netcdfoutput.pm to write map-variables by calling getdata (appended with mapextra):
        getdata.pl -f SDS-haven -o netcdf -d nc_map.nc -v xzeta,yzeta,xdep,ydep,sep,velu,velv,h
    waqpro.pl uses netcdfoutput.pm to write his-variables by calling getdata (appended with hisextra):
        getdata.pl -f SDS-haven -o netcdf -d nc_his.nc -v zwl,namwl,mwl,nwl,xzeta,yzeta,itdate,zcur,zcuru,zcurv,namc,ctr,fltr,namtra,ctrv,fltrv,namtrv
    get a list of all available variables in the SDS file:
        getdata.pl -f SDS-haven -v l

convert entire SDS file with getdata.pl (this gives you a list of variables but does not work yet):
    module load simona
    #mapvars_raw=$(getdata.pl -f SDS-haven -v l | grep 'TIME DEP' | grep MNMAXK | sed 's/\tREAL.*//' | sed 's/\tINT.*//' | tr '\n' ',')
    #hisvars_raw=$(getdata.pl -f SDS-haven -v l | grep 'TIME DEP' | grep -v MNMAXK | sed 's/\tREAL.*//' | sed 's/\tINT.*//' | tr '\n' ',')
    mapvars_raw=$(getdata.pl -f SDS-haven -v l | grep 'TIME DEP' | grep MNMAXK | grep -v "V\S*INT" | sed 's/\tREAL.*//' | sed 's/\tINT.*//' | tr '\n' ',')
    hisvars_raw=$(getdata.pl -f SDS-haven -v l | grep 'TIME DEP' | grep NO | sed 's/\tREAL.*//' | sed 's/\tINT.*//' | tr '\n' ',')
    getdata.pl -f SDS-haven -o netcdf -d nc_map.nc -v ${mapvars_raw%?}
    getdata.pl -f SDS-haven -o netcdf -d nc_his.nc -v ${hisvars_raw%?}

DCSM: convert existing waqua output to netcdf files via putty with:
    module load simona
    cd /p/1204257-dcsmzuno/2019/DCSMv6/A01
    getdata.pl -f SDS-A01 -v l
    getdata.pl -f SDS-A01 -v SEP,VELU,VELV,YZETA,XZETA -o netcdf -d SDS-A01_map
    getdata.pl -f SDS-A01 -v ZWL,ZCURU,ZCURV,NAMWL,NAMC -o netcdf -d SDS-A01_his
    http://simona.deltares.nl/release/doc/usedoc/getdata/getdata.pdf

RMM: convert existing waqua output to netcdf files via putty with:
    module load simona
    cd /p/11205258-006-kpp2020_rmm-g6/C_Work/07_WAQUAresultaten/j15
    getdata.pl -f SDS-riv_tba -v l
    getdata.pl -f SDS-riv_tba -v SEP,VELU,VELV,YZETA,XZETA -o netcdf -d SDS-riv_tba_map
    getdata.pl -f SDS-riv_tba -v ZWL,ZCURU,ZCURV,NAMWL,NAMC -o netcdf -d SDS-riv_tba_his
    http://simona.deltares.nl/release/doc/usedoc/getdata/getdata.pdf
    

"""

import os
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
plt.close('all')
import dfm_tools as dfmt

dir_testinput = r'c:\DATA\dfm_tools_testdata'
dir_output = '.'


file_nc = r'p:\1204257-dcsmzuno\2019\DCSMv6\A01\SDS-A01_map.nc'
file_nc = r'p:\archivedprojects\11205258-006-kpp2020_rmm-g6\C_Work\07_WAQUAresultaten\j15\SDS-riv_tba_map.nc'
file_nc = r'c:\DATA\dfm_tools_testdata\waqua_netcdf\SDS-haven_map.nc'



#MAP ZUNO
file_nc = r'p:\1204257-dcsmzuno\2019\DCSMv6\A01\SDS-A01_map.nc'
#vars_pd = get_ncvarproperties(file_nc=file_nc)

timestep = 1

file_nc = r'c:\DATA\dfm_tools_testdata\waqua_netcdf\SDS-haven_map.nc'
#file_nc = r'p:\archivedprojects\11205258-006-kpp2020_rmm-g6\C_Work\07_WAQUAresultaten\j15\SDS-riv_tba_map.nc'

if 1:
    ds = xr.open_dataset(file_nc)
    uds = dfmt.open_dataset_curvilinear(file_nc,varn_vert_lon='grid_x', varn_vert_lat='grid_y', ij_dims=['N','M'])
    fig,ax = plt.subplots()#figsize=figsize)
    uds.isel(TIME=timestep,LAYER=0).SEP.ugrid.plot(ax=ax, center=False, cmap='jet')
    #uds.grid.plot(ax=ax,linewidth=0.2)
    ax.set_aspect('equal')
    
else:
    data_xr = xr.open_dataset(file_nc)
    data_xr['grid_x_cen'] = data_xr['grid_x'].mean(dim='bounds').bfill(dim='M').bfill(dim='N').ffill(dim='M').ffill(dim='N') #bfill/ffill to avoid nans
    data_xr['grid_y_cen'] = data_xr['grid_y'].mean(dim='bounds').bfill(dim='M').bfill(dim='N').ffill(dim='M').ffill(dim='N') #bfill/ffill to avoid nans
    data_xr = data_xr.set_coords(['grid_x_cen','grid_y_cen'])
    data_nc_SEP = data_xr['SEP'].isel(TIME=timestep)
    fig, ax = plt.subplots()
    pc = data_nc_SEP.plot.pcolormesh(x='grid_x_cen',y='grid_y_cen',cmap='jet')


breakit
pc.set_clim([-0.1,0.1])
fig.tight_layout()
plt.savefig(os.path.join(dir_output,'waqua_DCSM_map_wl'))


#HIS ZUNO
file_nc = r'p:\1204257-dcsmzuno\2019\DCSMv6\A01\SDS-A01_his.nc'
#vars_pd = get_ncvarproperties(file_nc=file_nc)
data_xr = xr.open_dataset(file_nc)
stations_pd = data_xr.NAMWL.astype(str).to_pandas().str.strip()

fig, ax = plt.subplots(figsize=(16,7))
for iS in range(10):
    data_nc_ZWL = data_xr.ZWL.isel(STATION=iS).sel(TIME=slice('2018-12-22','2019-01-05'))
    ax.plot(data_nc_ZWL.TIME,data_nc_ZWL,label=stations_pd.iloc[iS], linewidth=1)
ax.legend()
ax.set_ylabel('%s (%s)'%(data_nc_ZWL.attrs['long_name'], data_nc_ZWL.attrs['units']))
time_ext = data_nc_ZWL.TIME[[0,-1]].to_numpy()
ax.set_xlim(time_ext)
plt.savefig(os.path.join(dir_output,'waqua_DCSM_his_ZWL'))




for RMM_name in ['RMM','RMMtestmodel']:
    if RMM_name=='RMM':
        file_nc_map = r'p:\archivedprojects\11205258-006-kpp2020_rmm-g6\C_Work\07_WAQUAresultaten\j15\SDS-riv_tba_map.nc'
        file_nc_his = r'p:\archivedprojects\11205258-006-kpp2020_rmm-g6\C_Work\07_WAQUAresultaten\j15\SDS-riv_tba_his.nc'
        timestep = 1
    elif RMM_name == 'RMMtestmodel':
        file_nc_map = r'c:\DATA\dfm_tools_testdata\waqua_netcdf\SDS-haven_map.nc'
        file_nc_his = r'c:\DATA\dfm_tools_testdata\waqua_netcdf\SDS-haven_his.nc'
        timestep = 10

    #MAP RMM
    file_nc = file_nc_map
    #vars_pd = get_ncvarproperties(file_nc=file_nc)
    
    data_xr = xr.open_dataset(file_nc)
    data_xr['grid_x_cen'] = data_xr['grid_x'].mean(dim='bounds').bfill(dim='M').bfill(dim='N').ffill(dim='M').ffill(dim='N') #bfill/ffill to avoid nans
    data_xr['grid_y_cen'] = data_xr['grid_y'].mean(dim='bounds').bfill(dim='M').bfill(dim='N').ffill(dim='M').ffill(dim='N') #bfill/ffill to avoid nans
    data_xr = data_xr.set_coords(['grid_x_cen','grid_y_cen'])
    data_nc_SEP = data_xr['SEP'].isel(TIME=timestep)
    data_nc_VELU = data_xr['VELU'].isel(TIME=timestep,LAYER=0)
    data_nc_VELV = data_xr['VELV'].isel(TIME=timestep,LAYER=0)
    
    
    fig, ax = plt.subplots(figsize=(16,7))
    pc = data_nc_SEP.plot.pcolormesh(x='grid_x_cen',y='grid_y_cen',cmap='jet')
    pc.set_clim([0,3])
    ax.set_aspect('equal')
    fig.tight_layout()
    plt.savefig(os.path.join(dir_output,'waqua_%s_map_wl'%(RMM_name)))

    fig, ax = plt.subplots(figsize=(16,7))
    data_nc_VELmagn = np.sqrt(data_nc_VELU**2 + data_nc_VELV**2)
    pc = data_nc_VELmagn.plot.pcolormesh(x='grid_x_cen',y='grid_y_cen',cmap='jet')
    pc.set_clim([0,1])
    if RMM_name=='RMM':
        thinning = 10
    else:
        thinning = 1
    data_nc_VELU_thin = data_nc_VELU.loc[::thinning,::thinning]
    data_nc_VELV_thin = data_nc_VELV.loc[::thinning,::thinning]
    ax.quiver(data_nc_VELU_thin.grid_x_cen, data_nc_VELU_thin.grid_y_cen, data_nc_VELU_thin, data_nc_VELV_thin, 
              color='w',scale=10)#,width=0.005)#, cmap='jet')
    if RMM_name=='RMM':
        ax.set_xlim([61000, 72000])
        ax.set_ylim([438000, 446000])
    ax.set_aspect('equal')
    fig.tight_layout()
    plt.savefig(os.path.join(dir_output,'waqua_%s_map_vel'%(RMM_name)))
    
    #HIS RMM
    file_nc = file_nc_his
    #vars_pd = get_ncvarproperties(file_nc=file_nc)
    
    data_xr = xr.open_dataset(file_nc)
    stations_pd = data_xr.NAMWL.astype(str).to_pandas().str.strip()
    
    fig, ax = plt.subplots(figsize=(16,7))
    for iS in range(10):
        data_nc_ZWL = data_xr.ZWL.isel(STATION=iS)
        ax.plot(data_nc_ZWL.TIME,data_nc_ZWL,label=stations_pd.iloc[iS], linewidth=1)
    ax.legend()
    ax.set_ylabel('%s (%s)'%(data_nc_ZWL.attrs['long_name'], data_nc_ZWL.attrs['units']))
    time_ext = data_nc_ZWL.TIME[[0,-1]].to_numpy()
    ax.set_xlim(time_ext)
    plt.savefig(os.path.join(dir_output,'waqua_%s_his_ZWL'%(RMM_name)))



