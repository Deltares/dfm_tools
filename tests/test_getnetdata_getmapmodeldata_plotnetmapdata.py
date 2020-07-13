# -*- coding: utf-8 -*-
"""
Created on Mon May 18 13:49:29 2020

@author: laan_st
"""
import os

dir_testinput = r'z:\OneDrive - Stichting Deltares\_Stendert\Projects\dfm_tools\testData'
file_nc = os.path.join(dir_testinput,r'DFM_sigma_curved_bend\DFM_OUTPUT_cb_3d\cb_3d_map.nc')
#file_nc = os.path.join(dir_testinput,r'DFM_3D_z_Grevelingen\computations\run01\DFM_OUTPUT_Grevelingen-FM\Grevelingen-FM_0000_map.nc')
#file_nc = r'p:\11205258-006-kpp2020_rmm-g6\C_Work\08_RMM_FMmodel\computations\run_156\DFM_OUTPUT_RMM_dflowfm\RMM_dflowfm_0000_map.nc'

dir_output = r'Z:/OneDrive - Stichting Deltares/_Stendert/Projects/dfm_tools/gitHubClone/tests/output'

"""
this test retrieves grid data, retrieves map data, and plots it
file_nc = os.path.join(dir_testinput,'DFM_3D_z_Grevelingen','computations','run01','DFM_OUTPUT_Grevelingen-FM','Grevelingen-FM_0000_map.nc')
dir_output = './test_output'
"""

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
plt.close('all')

from dfm_tools.get_nc import get_netdata, get_ncmodeldata, plot_netmapdata

if 'cb_3d_map' in file_nc:
    timestep = 3
    layer = 5
    clim_bl = None
    clim_wl = [-0.5,1]
    clim_sal = None
    clim_tem = None
elif 'Grevelingen-FM_0000_map' in file_nc:
    timestep = 3
    layer = 33
    clim_bl = None
    clim_wl = [-0.5,1]
    clim_sal = [28,30.2]
    clim_tem = [4,10]
elif 'RMM_dflowfm_0000_map' in file_nc:
    timestep = 50
    layer = None
    clim_bl = [-10,10]
    clim_wl = [-2,2]
    clim_sal = None
    clim_tem = None
else:
    raise Exception('ERROR: no settings provided for this mapfile')
    

#PLOT GRID
print('plot only grid from mapdata')
ugrid_all = get_netdata(file_nc=file_nc)#,multipart=False)
fig, ax = plt.subplots()
pc = plot_netmapdata(ugrid_all.verts, values=None, ax=None, linewidth=0.5, color="crimson", facecolor="None")
ax.set_xlabel('x-direction')
ax.set_ylabel('y-direction')
ax.set_aspect('equal')
plt.tight_layout()
plt.savefig(os.path.join(dir_output,'%s_grid'%(os.path.basename(file_nc).replace('.',''))))


#PLOT bedlevel
if not 'cb_3d_map' in file_nc:
    print('plot grid and values from mapdata (constantvalue, 1 dim)')
    ugrid = get_netdata(file_nc=file_nc)#,multipart=False)
    #iT = 3 #for iT in range(10):
    data_frommap = get_ncmodeldata(file_nc=file_nc, varname='mesh2d_flowelem_bl')#, multipart=False)
    data_frommap_flat = data_frommap.flatten()
    fig, ax = plt.subplots()
    pc = plot_netmapdata(ugrid.verts, values=data_frommap_flat, ax=None, linewidth=0.5, cmap="jet")
    pc.set_clim(clim_bl)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = fig.colorbar(pc, cax=cax)
    cbar.set_label('%s [%s]'%(data_frommap.var_ncvarobject.long_name, data_frommap.var_ncvarobject.units))
    ax.set_xlabel('%s [%s]'%(data_frommap.var_ncobject.variables.get(data_frommap.var_ncvarobject.coordinates.split()[0]).long_name,data_frommap.var_ncobject.variables.get(data_frommap.var_ncvarobject.coordinates.split()[0]).units))
    ax.set_ylabel('%s [%s]'%(data_frommap.var_ncobject.variables.get(data_frommap.var_ncvarobject.coordinates.split()[1]).long_name,data_frommap.var_ncobject.variables.get(data_frommap.var_ncvarobject.coordinates.split()[1]).units))
    ax.set_aspect('equal')
    plt.tight_layout()
    plt.savefig(os.path.join(dir_output,'%s_%s'%(os.path.basename(file_nc).replace('.',''),data_frommap.var_varname)))
    

#PLOT water level on map
print('plot grid and values from mapdata (waterlevel, 2dim)')
data_frommap = get_ncmodeldata(file_nc=file_nc, varname='mesh2d_s1', timestep=timestep)#, multipart=False)
data_frommap_flat = data_frommap.flatten()
fig, ax = plt.subplots()
pc = plot_netmapdata(ugrid_all.verts, values=data_frommap_flat, ax=None, linewidth=0.5, cmap="jet")
pc.set_clim(clim_wl)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
cbar = fig.colorbar(pc, cax=cax)
cbar.set_label('%s [%s]'%(data_frommap.var_ncvarobject.long_name, data_frommap.var_ncvarobject.units))
ax.set_xlabel('%s [%s]'%(data_frommap.var_ncobject.variables.get(data_frommap.var_ncvarobject.coordinates.split()[0]).long_name,data_frommap.var_ncobject.variables.get(data_frommap.var_ncvarobject.coordinates.split()[0]).units))
ax.set_ylabel('%s [%s]'%(data_frommap.var_ncobject.variables.get(data_frommap.var_ncvarobject.coordinates.split()[1]).long_name,data_frommap.var_ncobject.variables.get(data_frommap.var_ncvarobject.coordinates.split()[1]).units))
ax.set_aspect('equal')
plt.tight_layout()
plt.savefig(os.path.join(dir_output,'%s_%s'%(os.path.basename(file_nc).replace('.',''),data_frommap.var_varname)))

#PLOT var layer on map
if not 'RMM_dflowfm_0000_map' in file_nc:
    print('plot grid and values from mapdata (salinity on layer, 3dim)')
    data_frommap = get_ncmodeldata(file_nc=file_nc, varname='mesh2d_sa1', timestep=timestep, layer=layer)#, multipart=False)
    data_frommap_flat = data_frommap.flatten()
    fig, ax = plt.subplots()
    pc = plot_netmapdata(ugrid_all.verts, values=data_frommap_flat, ax=None, linewidth=0.5, cmap="jet")
    pc.set_clim(clim_sal)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = fig.colorbar(pc, cax=cax)
    cbar.set_label('%s [%s]'%(data_frommap.var_ncvarobject.long_name, data_frommap.var_ncvarobject.units))
    ax.set_xlabel('%s [%s]'%(data_frommap.var_ncobject.variables.get(data_frommap.var_ncvarobject.coordinates.split()[0]).long_name,data_frommap.var_ncobject.variables.get(data_frommap.var_ncvarobject.coordinates.split()[0]).units))
    ax.set_ylabel('%s [%s]'%(data_frommap.var_ncobject.variables.get(data_frommap.var_ncvarobject.coordinates.split()[1]).long_name,data_frommap.var_ncobject.variables.get(data_frommap.var_ncvarobject.coordinates.split()[1]).units))
    ax.set_aspect('equal')
    plt.tight_layout()
    plt.savefig(os.path.join(dir_output,'%s_%s'%(os.path.basename(file_nc).replace('.',''),data_frommap.var_varname)))

    print('plot grid and values from mapdata (temperature on layer, 3dim)')
    data_frommap = get_ncmodeldata(file_nc=file_nc, varname='mesh2d_tem1', timestep=timestep, layer=layer)#, multipart=False)
    data_frommap_flat = data_frommap.flatten()
    fig, ax = plt.subplots()
    pc = plot_netmapdata(ugrid_all.verts, values=data_frommap_flat, ax=None, linewidth=0.5, cmap="jet")
    pc.set_clim(clim_tem)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = fig.colorbar(pc, cax=cax)
    cbar.set_label('%s [%s]'%(data_frommap.var_ncvarobject.long_name, data_frommap.var_ncvarobject.units))
    ax.set_xlabel('%s [%s]'%(data_frommap.var_ncobject.variables.get(data_frommap.var_ncvarobject.coordinates.split()[0]).long_name,data_frommap.var_ncobject.variables.get(data_frommap.var_ncvarobject.coordinates.split()[0]).units))
    ax.set_ylabel('%s [%s]'%(data_frommap.var_ncobject.variables.get(data_frommap.var_ncvarobject.coordinates.split()[1]).long_name,data_frommap.var_ncobject.variables.get(data_frommap.var_ncvarobject.coordinates.split()[1]).units))
    ax.set_aspect('equal')
    plt.tight_layout()
    plt.savefig(os.path.join(dir_output,'%s_%s'%(os.path.basename(file_nc).replace('.',''),data_frommap.var_varname)))