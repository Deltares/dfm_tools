# -*- coding: utf-8 -*-
"""
Created on Wed Oct 27 14:50:46 2021

@author: veenstra
"""

import os
import matplotlib.pyplot as plt
plt.close('all')

from dfm_tools.get_nc import get_netdata, plot_netmapdata

dir_testinput = r'c:\DATA\dfm_tools_testdata'
dir_output = '.'

file_nc_list = [os.path.join(dir_testinput,r'DFM_3D_z_Grevelingen\computations\run01\Grevelingen_FM_grid_20190603_net.nc'),
                os.path.join(dir_testinput,'vanNithin','myortho3_RGFGRID_net.nc'),
                os.path.join(dir_testinput,'DFM_fou_RMM','rooster_rmm_v1p5_net.nc'),
                ]

for file_nc in file_nc_list:

    print('plot only grid from net.nc')
    ugrid = get_netdata(file_nc=file_nc)
    fig, ax = plt.subplots()
    plot_netmapdata(ugrid.verts, values=None, ax=None, linewidth=0.5, color="crimson", facecolor="None")
    ax.set_aspect('equal')
    plt.savefig(os.path.join(dir_output,os.path.basename(file_nc).replace('.','')))
    
    
