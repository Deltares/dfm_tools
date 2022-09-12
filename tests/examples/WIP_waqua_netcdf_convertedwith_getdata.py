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
    
OSR: convert existing waqua output to netcdf files via putty with:
    module load simona
    cd /p/archivedprojects/1230049-zoutlastbeperking/Gaten_langsdam/Simulaties/OSR-model_GatenLangsdam/berekeningen/run7
    getdata.pl -f SDS-nsctri -v l
    getdata.pl -f SDS-nsctri -v SEP,VELU,VELV -o netcdf -d SDS-nsctri_map
    getdata.pl -f SDS-nsctri -v ZWL,ZCURU,ZCURV,NAMWL,NAMC -o netcdf -d SDS-nsctri_his
    http://simona.deltares.nl/release/doc/usedoc/getdata/getdata.pdf
    #this file should be recreated with YZETA,XZETA added to map

RMM: convert existing waqua output to netcdf files via putty with:
    module load simona
    cd /p/11205258-006-kpp2020_rmm-g6/C_Work/07_WAQUAresultaten/j15
    getdata.pl -f SDS-riv_tba -v l
    getdata.pl -f SDS-riv_tba -v SEP,VELU,VELV,YZETA,XZETA -o netcdf -d SDS-riv_tba_map
    getdata.pl -f SDS-riv_tba -v ZWL,ZCURU,ZCURV,NAMWL,NAMC -o netcdf -d SDS-riv_tba_his
    http://simona.deltares.nl/release/doc/usedoc/getdata/getdata.pdf
    

"""

import os
import datetime as dt
import numpy as np
import matplotlib.pyplot as plt
plt.close('all')

from dfm_tools.get_nc import get_ncmodeldata#, get_netdata, plot_netmapdata
from dfm_tools.get_nc_helpers import get_ncvardimlist, get_hisstationlist#, get_varname_fromnc
import xarray as xr

dir_testinput = r'c:\DATA\dfm_tools_testdata'
dir_output = '.'


#MAP ZUNO
file_nc = r'p:\1204257-dcsmzuno\2019\DCSMv6\A01\SDS-A01_map.nc'
vars_pd, dims_pd = get_ncvardimlist(file_nc=file_nc)

data_nc_x = get_ncmodeldata(file_nc=file_nc, varname='grid_x')
data_nc_y = get_ncmodeldata(file_nc=file_nc, varname='grid_y')
data_nc_SEP = get_ncmodeldata(file_nc=file_nc, varname='SEP',timestep='all')
#data_nc_VELU = get_ncmodeldata(file_nc=file_nc, varname='VELU',timestep='all')
#data_nc_VELV = get_ncmodeldata(file_nc=file_nc, varname='VELV',timestep='all')
data_nc_xcen = np.mean(data_nc_x, axis=2)
data_nc_ycen = np.mean(data_nc_y, axis=2)

fig, ax = plt.subplots()
#vel_x, vel_y, vel_magn, direction_naut_deg = uva2xymagdeg(u=data_nc_U1[timestep,90,:,:],v=data_nc_V1[timestep,90,:,:],alpha=data_nc_ALFAS)
#pc = ax.pcolor(data_nc_XZ,data_nc_YZ,direction_naut_deg,cmap='jet')
#pc.set_clim([0,360])
timestep=0
pc = ax.pcolor(data_nc_xcen,data_nc_ycen,data_nc_SEP[timestep,:,:],cmap='jet')
pc.set_clim([-0.1,0.1])
cbar = fig.colorbar(pc, ax=ax)
cbar.set_label('%s (%s)'%(data_nc_SEP.var_varname, data_nc_SEP.var_ncattrs['units']))
ax.set_title('t=%d (%s)'%(timestep, data_nc_SEP.var_times.iloc[timestep]))
ax.set_aspect('equal')
#ax.quiver(data_nc_XZ[::2,::2], data_nc_YZ[::2,::2], vel_x[::2,::2], vel_y[::2,::2], 
#          scale=3,color='w',width=0.005)#, edgecolor='face', cmap='jet')
fig.tight_layout()
plt.savefig(os.path.join(dir_output,'waqua_DCSM_map_wl'))


#HIS ZUNO
file_nc = r'p:\1204257-dcsmzuno\2019\DCSMv6\A01\SDS-A01_his.nc'
vars_pd, dims_pd = get_ncvardimlist(file_nc=file_nc)
data_xr = xr.open_dataset(file_nc)

#data_nc_NAMWL = get_hisstationlist(file_nc=file_nc, varname='NAMWL')
#data_nc_NAMC = get_hisstationlist(file_nc=file_nc, varname='NAMC')
#data_nc_ZWL = get_ncmodeldata(file_nc=file_nc, varname='ZWL',timestep='all',station='all')
#data_nc_ZCURU = get_ncmodeldata(file_nc=file_nc, varname='ZCURU',timestep='all',station='all')
#data_nc_ZCURV = get_ncmodeldata(file_nc=file_nc, varname='ZCURV',timestep='all',station='all')

stations_pd = data_xr.NAMWL.astype(str).to_pandas().str.strip()

fig, ax = plt.subplots(figsize=(16,7))
for iS in range(10):
    data_nc_ZWL = data_xr.ZWL.isel(STATION=iS).sel(TIME=slice('2018-12-22','2019-01-05'))
    ax.plot(data_nc_ZWL.TIME,data_nc_ZWL,label=stations_pd.iloc[iS], linewidth=1)
ax.legend()
ax.set_ylabel('%s (%s)'%(data_nc_ZWL.attrs['long_name'], data_nc_ZWL.attrs['units']))
time_ext = data_nc_ZWL.TIME[[0,-1]].to_numpy()
ax.set_xlim(time_ext)



"""
#MAP OSR
file_nc = r'p:\11205258-006-kpp2020_rmm-g6\C_Work\ZZ_Jelmer\SDS-nsctri_map.nc'
vars_pd, dims_pd = get_ncvardimlist(file_nc=file_nc)

data_nc_x = get_ncmodeldata(file_nc=file_nc, varname='grid_x')
data_nc_y = get_ncmodeldata(file_nc=file_nc, varname='grid_y')
data_nc_xcen = np.mean(data_nc_x, axis=2)
data_nc_ycen = np.mean(data_nc_y, axis=2)

timestep = 10
data_nc_SEP = get_ncmodeldata(file_nc=file_nc, varname='SEP',timestep=timestep)
data_nc_VELU = get_ncmodeldata(file_nc=file_nc, varname='VELU',timestep=timestep, layer=9)
data_nc_VELV = get_ncmodeldata(file_nc=file_nc, varname='VELV',timestep=timestep, layer=9)

fig, ax = plt.subplots()
pc = ax.pcolor(data_nc_xcen,data_nc_ycen,data_nc_SEP[0,:,:],cmap='jet')
pc.set_clim([0,2])
cbar = fig.colorbar(pc, ax=ax)
cbar.set_label('%s (%s)'%(data_nc_SEP.var_varname, data_nc_SEP.var_ncattrs['units']))
ax.set_title('t=%d (%s)'%(timestep, data_nc_SEP.var_times.loc[timestep]))
ax.set_aspect('equal')
fig.tight_layout()
plt.savefig(os.path.join(dir_output,'waqua_OSR_map_wl'))

fig, ax = plt.subplots()
vel_magn = np.sqrt(data_nc_VELU**2 + data_nc_VELV**2)
pc = ax.pcolor(data_nc_xcen,data_nc_ycen,vel_magn[0,:,:,0],cmap='jet')
pc.set_clim([0,1])
cbar = fig.colorbar(pc, ax=ax)
cbar.set_label('velocity magnitude (%s)'%(data_nc_VELU.var_ncattrs['units']))
ax.set_title('t=%d (%s)'%(timestep, data_nc_VELU.var_times.loc[timestep]))
ax.set_aspect('equal')
thinning = 10
ax.quiver(data_nc_xcen[::thinning,::thinning], data_nc_ycen[::thinning,::thinning], data_nc_VELU[0,::thinning,::thinning,0], data_nc_VELV[0,::thinning,::thinning,0], 
          color='w',scale=10)#,width=0.005)#, edgecolor='face', cmap='jet')
ax.set_xlim([58000, 66000])
ax.set_ylim([442000, 448000])
fig.tight_layout()
plt.savefig(os.path.join(dir_output,'waqua_OSR_map_vel'))

#HIS OSR
file_nc = r'p:\11205258-006-kpp2020_rmm-g6\C_Work\ZZ_Jelmer\SDS-nsctri_his.nc'
vars_pd, dims_pd = get_ncvardimlist(file_nc=file_nc)
data_nc_NAMWL = get_hisstationlist(file_nc=file_nc, varname='NAMWL')
#data_nc_NAMC = get_hisstationlist(file_nc=file_nc, varname='NAMC')
data_nc_ZWL = get_ncmodeldata(file_nc=file_nc, varname='ZWL',timestep='all',station='all')
#data_nc_ZCURU = get_ncmodeldata(file_nc=file_nc, varname='ZCURU',timestep='all',station='all',layer='all')
#data_nc_ZCURV = get_ncmodeldata(file_nc=file_nc, varname='ZCURV',timestep='all',station='all',layer='all')

fig, ax = plt.subplots(figsize=(16,7))
for iS in range(10):
    ax.plot(data_nc_ZWL.var_times,data_nc_ZWL[:,iS],label=data_nc_NAMWL['NAMWL'].iloc[iS], linewidth=1)
ax.legend()
ax.set_ylabel('%s (%s)'%(data_nc_ZWL.var_varname, data_nc_ZWL.var_ncattrs['units']))
ax.set_xlim([data_nc_ZWL.var_times[0],data_nc_ZWL.var_times[0]+dt.timedelta(days=14)])
plt.savefig(os.path.join(dir_output,'waqua_OSR_his_ZWL'))
"""


#MAP RMM
RMM_names = ['RMM','RMMtestmodel']
for RMM_name in RMM_names:
    if RMM_name=='RMM':
        file_nc_map = r'p:\11205258-006-kpp2020_rmm-g6\C_Work\07_WAQUAresultaten\j15\SDS-riv_tba_map.nc'
        file_nc_his = r'p:\11205258-006-kpp2020_rmm-g6\C_Work\07_WAQUAresultaten\j15\SDS-riv_tba_his.nc'
        timestep = 1
    elif RMM_name == 'RMMtestmodel':
        file_nc_map = r'c:\DATA\dfm_tools_testdata\waqua_netcdf\SDS-haven_map.nc'
        file_nc_his = r'c:\DATA\dfm_tools_testdata\waqua_netcdf\SDS-haven_his.nc'
        timestep = 10

    file_nc = file_nc_map
    vars_pd, dims_pd = get_ncvardimlist(file_nc=file_nc)
    
    data_nc_x = get_ncmodeldata(file_nc=file_nc, varname='grid_x')
    data_nc_y = get_ncmodeldata(file_nc=file_nc, varname='grid_y')
    data_nc_xcen = np.mean(data_nc_x, axis=2)
    data_nc_ycen = np.mean(data_nc_y, axis=2)
    
    data_nc_SEP = get_ncmodeldata(file_nc=file_nc, varname='SEP',timestep=timestep)
    data_nc_VELU = get_ncmodeldata(file_nc=file_nc, varname='VELU',timestep=timestep, layer=0)
    data_nc_VELV = get_ncmodeldata(file_nc=file_nc, varname='VELV',timestep=timestep, layer=0)
    
    fig, ax = plt.subplots(figsize=(16,7))
    pc = ax.pcolor(data_nc_xcen,data_nc_ycen,data_nc_SEP[0,:,:],cmap='jet')
    pc.set_clim([0,3])
    cbar = fig.colorbar(pc, ax=ax)
    cbar.set_label('%s (%s)'%(data_nc_SEP.var_varname, data_nc_SEP.var_ncattrs['units']))
    ax.set_title('t=%d (%s)'%(timestep, data_nc_SEP.var_times.loc[timestep]))
    ax.set_aspect('equal')
    fig.tight_layout()
    plt.savefig(os.path.join(dir_output,'waqua_%s_map_wl'%(RMM_name)))

    if RMM_name=='RMM':
        fig, ax = plt.subplots()
    else:
        fig, ax = plt.subplots(figsize=(16,7))
    vel_magn = np.sqrt(data_nc_VELU**2 + data_nc_VELV**2)
    pc = ax.pcolor(data_nc_xcen,data_nc_ycen,vel_magn[0,:,:,0],cmap='jet')
    pc.set_clim([0,1])
    cbar = fig.colorbar(pc, ax=ax)
    cbar.set_label('velocity magnitude (%s)'%(data_nc_VELU.var_ncattrs['units']))
    ax.set_title('t=%d (%s)'%(timestep, data_nc_VELU.var_times.loc[timestep]))
    ax.set_aspect('equal')
    if RMM_name=='RMM':
        thinning = 10
    else:
        thinning = 1
    ax.quiver(data_nc_xcen[::thinning,::thinning], data_nc_ycen[::thinning,::thinning], data_nc_VELU[0,::thinning,::thinning,0], data_nc_VELV[0,::thinning,::thinning,0], 
              color='w',scale=10)#,width=0.005)#, edgecolor='face', cmap='jet')
    if RMM_name=='RMM':
        ax.set_xlim([61000, 72000])
        ax.set_ylim([438000, 446000])
    fig.tight_layout()
    plt.savefig(os.path.join(dir_output,'waqua_%s_map_vel'%(RMM_name)))
    
    #HIS RMM
    file_nc = file_nc_his
    vars_pd, dims_pd = get_ncvardimlist(file_nc=file_nc)
    #data_nc_NAMWL = get_hisstationlist(file_nc=file_nc, varname='NAMWL')
    #data_nc_NAMC = get_hisstationlist(file_nc=file_nc, varname='NAMC')
    #data_nc_ZWL = get_ncmodeldata(file_nc=file_nc, varname='ZWL',timestep='all',station='all')
    #data_nc_ZCURU = get_ncmodeldata(file_nc=file_nc, varname='ZCURU',timestep='all',station='all',layer='all')
    #data_nc_ZCURV = get_ncmodeldata(file_nc=file_nc, varname='ZCURV',timestep='all',station='all',layer='all')
    
    data_xr = xr.open_dataset(file_nc)
    stations_pd = data_xr.NAMWL.astype(str).to_pandas().str.strip()
    
    fig, ax = plt.subplots(figsize=(16,7))
    for iS in range(10):
        data_nc_ZWL = data_xr.ZWL.isel(STATION=iS)
        if RMM_name=='RMM':
            data_nc_ZWL = data_nc_ZWL.sel(TIME=slice('2014-01-15','2014-01-29'))
        ax.plot(data_nc_ZWL.TIME,data_nc_ZWL,label=stations_pd.iloc[iS], linewidth=1)
    ax.legend()
    ax.set_ylabel('%s (%s)'%(data_nc_ZWL.attrs['long_name'], data_nc_ZWL.attrs['units']))
    time_ext = data_nc_ZWL.TIME[[0,-1]].to_numpy()
    ax.set_xlim(time_ext)
    plt.savefig(os.path.join(dir_output,'waqua_%s_his_ZWL'%(RMM_name)))



