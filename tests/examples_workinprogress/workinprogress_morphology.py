# -*- coding: utf-8 -*-
"""
Created on Wed Oct 27 22:01:16 2021

@author: veenstra
"""

import os
import matplotlib.pyplot as plt
plt.close('all')
import xarray as xr
import dfm_tools as dfmt

dir_testinput = r'c:\DATA\dfm_tools_testdata'
dir_output = '.'


#MAPFILE
file_nc = r'p:\archivedprojects\11203869-morwaqeco3d\05-Tidal_inlet\02_FM_201910\FM_MF10_Max_30s\fm\DFM_OUTPUT_inlet\inlet_map.nc'
data_frommap_merged = dfmt.open_partitioned_dataset(file_nc)

vars_pd = dfmt.get_ncvarproperties(file_nc=file_nc)
vars_pd_sel = vars_pd[vars_pd['long_name'].str.contains('transport')]
#vars_pd_sel = vars_pd[vars_pd['dimensions'].str.contains('mesh2d_nFaces') & vars_pd['long_name'].str.contains('wave')]

varname = 'mesh2d_mor_bl'
var_clims = [-50,0]
var_longname = vars_pd['long_name'][vars_pd.index==varname].iloc[0]
fig, axs = plt.subplots(3,1, figsize=(6,9))

ax = axs[0]
data_frommap_0 = data_frommap_merged[varname].isel(time=0)
pc = data_frommap_0.ugrid.plot(ax=ax, linewidth=0.5, cmap='jet', clim=var_clims)

ax = axs[1]
data_frommap_end = data_frommap_merged[varname].isel(time=-1)
pc = data_frommap_end.ugrid.plot(ax=ax, linewidth=0.5, cmap='jet', clim=var_clims)

ax = axs[2]
data_diff = data_frommap_end-data_frommap_0
pc = data_diff.ugrid.plot(ax=ax, linewidth=0.5, cmap='jet', clim=[-3,3])

for ax in axs:
    ax.set_aspect('equal')
    #ax.set_ylim(val_ylim)
fig.tight_layout()
plt.savefig(os.path.join(dir_output,'%s_%s'%(os.path.basename(file_nc).replace('.',''), varname)))

varname = 'mesh2d_hwav'
var_longname = vars_pd['long_name'][vars_pd.index==varname].iloc[0]
fig, ax = plt.subplots(1,1)
data_frommap = data_frommap_merged[varname].isel(time=-1)
pc = data_frommap.ugrid.plot(ax=ax, linewidth=0.5, cmap='jet')
ax.set_aspect('equal')

fig.tight_layout()
plt.savefig(os.path.join(dir_output,'%s_%s'%(os.path.basename(file_nc).replace('.',''), varname)))



#WAVM FILE
file_nc = r'p:\archivedprojects\11203869-morwaqeco3d\05-Tidal_inlet\02_FM_201910\FM_MF10_Max_30s\wave\wavm-inlet.nc'
vars_pd = dfmt.get_ncvarproperties(file_nc=file_nc)
vars_pd_sel = vars_pd[vars_pd['long_name'].str.contains('dissi')]
#vars_pd_sel = vars_pd[vars_pd['dimensions'].str.contains('mesh2d_nFaces') & vars_pd['long_name'].str.contains('wave')]

data_xr = xr.open_dataset(file_nc)

#plt.close('all')
varname_list = ['hsign', 'dir', 'period', 'dspr']
var_clim = [[0,2], [0,360], [0,7.5], [0,35], [0,20]]
for iV, varname in enumerate(varname_list):
    var_longname = vars_pd['long_name'][vars_pd.index==varname].iloc[0]
    
    fig, axs = plt.subplots(1,2, figsize=(12,7))
    fig.suptitle('%s (%s)'%(varname, var_longname))

    timestep = 10
    ax = axs[0]
    pc = data_xr[varname].isel(time=timestep).plot(ax=ax, cmap='jet')
    pc.set_clim(var_clim[iV])
    ax.set_aspect('equal')
    
    timestep = -1
    ax = axs[1]
    pc = data_xr[varname].isel(time=timestep).plot(ax=ax, cmap='jet')
    pc.set_clim(var_clim[iV])
    ax.set_aspect('equal')
    
    fig.tight_layout()
    plt.savefig(os.path.join(dir_output,'%s_%s'%(os.path.basename(file_nc).replace('.',''), varname)))
    


