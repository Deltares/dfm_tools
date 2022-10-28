# -*- coding: utf-8 -*-
"""
Created on Fri Oct 28 06:40:03 2022

@author: veenstra
"""


import os
import matplotlib.pyplot as plt
plt.close('all')
import xarray as xr
import pandas as pd
import datetime as dt
import glob

from dfm_tools.get_nc import get_netdata, get_ncmodeldata, plot_netmapdata
from dfm_tools.regulargrid import scatter_to_regulargrid
from dfm_tools.get_nc_helpers import get_ncvarproperties

overwrite = False

dir_testinput = r'c:\DATA\dfm_tools_testdata'
file_nc_star = os.path.join(dir_testinput,r'DFM_3D_z_Grevelingen\computations\run01\DFM_OUTPUT_Grevelingen-FM\Grevelingen-FM_00*_map.nc')
#file_nc_star = r'p:\11206813-006-kpp2021_rmm-2d\C_Work\31_RMM_FMmodel\computations\model_setup\run_206\results\RMM_dflowfm_00*_map.nc'

file_nc_list = glob.glob(file_nc_star)#[:1]

for file_nc in file_nc_list:
    dir_output = os.path.join(os.path.dirname(file_nc),'monthlymeans')
    if not os.path.exists(dir_output):
        os.mkdir(dir_output)
    postfix = '_'.join(os.path.basename(file_nc).split('_')[-2:])
    file_out = os.path.join(dir_output,os.path.basename(file_nc).replace(postfix,f'monthlymean_{postfix}'))
    if not overwrite:
        continue
    
    print('opening dataset')
    dtstart = dt.datetime.now()
    data_xr = xr.open_dataset(file_nc, mask_and_scale=False) #TODO: mask_and_scale is necessary to avoid nan (toohigh) quadrilateral values
    for varname in data_xr.data_vars: #TODO: this should be avoided, can also be done with mask_and_scale argument, but is that desireable?
        if (data_xr[varname].encoding['dtype']==int) & (data_xr[varname].dtype!=int):
            print(f'converting {varname} back to int dtype')
            data_xr[varname] = data_xr[varname].astype(int)
    data_xr_var = data_xr['mesh2d_sa1']
    time_passed = (dt.datetime.now()-dtstart).total_seconds()
    print(f'>>time passed: {time_passed:.2f} sec')
    
    print('computing monthly means')
    dtstart = dt.datetime.now()
    if False: #on unique month numbers
        #data_xr_monthmean = data_xr_var.groupby('time.month').mean()
        data_xr_monthmean = data_xr.groupby('time.month').mean()
    else: # on unique year+month combinations
        data_xr_monthmean = data_xr.copy()
        #data_xr_monthmean = data_xr_var.copy()
        data_xr_monthmean = data_xr_monthmean.resample(time='MS').mean(dim='time')
    #data_xr_monthmean.time.encoding = data_xr.time.encoding
    time_passed = (dt.datetime.now()-dtstart).total_seconds()
    print(f'>>time passed: {time_passed:.2f} sec')

    for varname in data_xr.data_vars: #TODO: avoid adding time dimension to time-independent variables: https://github.com/pydata/xarray/issues/2145
        if data_xr[varname].dims != data_xr_monthmean[varname].dims:
            #print(f'{varname} dims changed from {data_xr[varname].dims} to {data_xr_monthmean[varname].dims}')
            data_xr_monthmean[varname] = data_xr_monthmean[varname].isel(time=0) #selecting first time of variables that did not have time dimension before
            if data_xr[varname].dims != data_xr_monthmean[varname].dims:
                print(f'{varname} dims changed from {data_xr[varname].dims} to {data_xr_monthmean[varname].dims}')
    
    data_xr_monthmean.to_netcdf(file_out)

timestep = 0
layer = 33
clim_bl = None
clim_wl = [-0.5,1]
clim_sal = [28,30.2]
clim_tem = [4,10]
crs = "EPSG:28992"

ugrid_all = get_netdata(file_out)


print('plot grid and values from mapdata (salinity on layer, 3dim, on cell centers)')
try:
    data_frommap = get_ncmodeldata(file_nc=file_out, varname='mesh2d_sa1', timestep=timestep, layer=layer)#, multipart=False)
except:
    data_frommap = get_ncmodeldata(file_nc=file_out, varname='mesh2d_sa1', timestep=timestep)#, multipart=False)
fig, ax = plt.subplots()
pc = plot_netmapdata(ugrid_all.verts, values=data_frommap, ax=None, linewidth=0.5, cmap="jet", edgecolor='face')
pc.set_clim(clim_sal)
cbar = fig.colorbar(pc, ax=ax)
cbar.set_label('%s [%s]'%(data_frommap.var_ncattrs['long_name'], data_frommap.var_ncattrs['units']))
#ax.set_xlabel('x [%s]'%(gridunits))
#ax.set_ylabel('y [%s]'%(gridunits))
ax.set_aspect('equal')
fig.tight_layout()
#fig.savefig(os.path.join(dir_output,'%s_mesh2d_sa1'%(os.path.basename(file_nc).replace('.',''))))

