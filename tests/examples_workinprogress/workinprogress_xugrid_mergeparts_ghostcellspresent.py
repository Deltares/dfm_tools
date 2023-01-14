# -*- coding: utf-8 -*-
"""
Created on Sat Jan 14 20:22:20 2023

@author: veenstra
"""


import os
import matplotlib.pyplot as plt
plt.close('all')
import glob
import dfm_tools as dfmt

dir_testinput = r'c:\DATA\dfm_tools_testdata'

file_nc = os.path.join(dir_testinput,'DFM_3D_z_Grevelingen','computations','run01','DFM_OUTPUT_Grevelingen-FM','Grevelingen-FM_0*_map.nc')
#file_nc = glob.glob(file_nc)[:4]
file_nc = [file_nc.replace('_0*_',f'_{idno:04d}_') for idno in [2,3]] #most leftward partitions

varn = 'mesh2d_flowelem_domain' # 'mesh2d_flowelem_bl' 'mesh2d_flowelem_domain'
alpha = 1.0
edgecolor = 'grey'

data_frommap_merged = dfmt.open_partitioned_dataset(file_nc, merge_xugrid=False)
fig, ax_input = plt.subplots()
pc = data_frommap_merged[varn].ugrid.plot(edgecolor=edgecolor,cmap='jet', alpha=alpha, linewidth=0.1)

data_frommap_merged = dfmt.open_partitioned_dataset(file_nc, merge_xugrid=True)
fig, ax_input = plt.subplots()
pc = data_frommap_merged[varn].ugrid.plot(edgecolor=edgecolor,cmap='jet', alpha=alpha, linewidth=0.1)

