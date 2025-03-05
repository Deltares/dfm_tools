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

dir_output = '.'

file_nc_map = r'p:\archivedprojects\11205258-006-kpp2020_rmm-g6\C_Work\07_WAQUAresultaten\j15\SDS-riv_tba_map.nc'
file_nc_his = r'p:\archivedprojects\11205258-006-kpp2020_rmm-g6\C_Work\07_WAQUAresultaten\j15\SDS-riv_tba_his.nc'
timestep = 1
resolution = 300
figsize = (16,7)

#MAP RMM
uds = dfmt.open_dataset_curvilinear(
    file_nc_map, 
    x_dim='N',
    y_dim='M',
    x_bounds='grid_x',
    y_bounds='grid_y', 
    )
uds_sel = uds.isel(TIME=timestep,LAYER=0)
uds_rastered = dfmt.rasterize_ugrid(uds_sel,resolution=resolution)

fig,ax = plt.subplots(figsize=figsize)
pc = uds_sel.SEP.ugrid.plot(ax=ax, center=False, cmap='jet')
pc.set_clim([0,3])
ax.set_aspect('equal')
fig.tight_layout()
fig.savefig(os.path.join(dir_output,'waqua_RMM_map_wl'))

fig, ax = plt.subplots(figsize=figsize)
data_nc_VELmagn = np.sqrt(uds_sel.VELU**2 + uds_sel.VELV**2)
pc = data_nc_VELmagn.ugrid.plot(ax=ax, center=False, cmap='jet')
pc.set_clim([0,1])
ax.quiver(uds_rastered.x, uds_rastered.y, uds_rastered.VELU, uds_rastered.VELV, 
          color='w',scale=10)
ax.set_xlim(61000, 72000)
ax.set_ylim(438000, 446000)
ax.set_aspect('equal')
fig.tight_layout()
fig.savefig(os.path.join(dir_output,'waqua_RMM_map_vel'))

#HIS RMM
data_xr = xr.open_dataset(file_nc_his)
stations_pd = data_xr.NAMWL.astype(str).to_pandas().str.strip()

fig, ax = plt.subplots(figsize=figsize)
for iS in range(10):
    data_nc_ZWL = data_xr.ZWL.isel(STATION=iS)
    ax.plot(data_nc_ZWL.TIME,data_nc_ZWL,label=stations_pd.iloc[iS], linewidth=1)
ax.legend()
ax.set_ylabel('%s (%s)'%(data_nc_ZWL.attrs['long_name'], data_nc_ZWL.attrs['units']))
time_ext = data_nc_ZWL.TIME[[0,-1]].to_numpy()
ax.set_xlim(time_ext)
plt.savefig(os.path.join(dir_output,'waqua_RMM_his_ZWL'))


