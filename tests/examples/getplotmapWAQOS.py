# -*- coding: utf-8 -*-
"""
Created on Wed Oct 27 14:21:28 2021

@author: veenstra
"""

import os
import matplotlib.pyplot as plt
plt.close('all')

from dfm_tools.get_nc import get_netdata, get_ncmodeldata, plot_netmapdata

dir_testinput = r'c:\DATA\dfm_tools_testdata'
file_nc = os.path.join(dir_testinput,'oost_tracer_2_map.nc')

ugrid = get_netdata(file_nc=file_nc)

print('plot grid and values from mapdata (constantvalue, 1 dim)')
if 'oost_tracer_map' in file_nc:
    var_names = ['FColi1','HIWAI1','mspaf1','Pharma1'] #nieuwe file, te veel dimensies
    var_clims = [None,[0,100000000000],None,[0,10000]]
elif 'oost_tracer_2_map' in file_nc:
    var_names = ['mesh2d_FColi','mesh2d_HIWAI','mesh2d_Pharma'] #oude file
    var_clims = [None,[0,100000000000],[0,10000]]
else:
    raise Exception('ERROR: no settings provided for this mapfile')

for var_name, var_clim in zip(var_names, var_clims):
    fig, ax = plt.subplots()
    data_frommap = get_ncmodeldata(file_nc=file_nc, varname=var_name)#, multipart=False)
    pc = plot_netmapdata(ugrid.verts, values=data_frommap, ax=None, linewidth=0.5, cmap="jet")
    if var_clim != None:
        pc.set_clim(var_clim)
    fig.colorbar(pc, ax=ax)
    ax.set_aspect('equal')
    ax.set_xlabel(var_name)
    plt.savefig('%s_%s'%(os.path.basename(file_nc).replace('.',''),var_name))
