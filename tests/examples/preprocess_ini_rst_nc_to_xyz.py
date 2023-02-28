# -*- coding: utf-8 -*-
"""
Created on Tue Jun 15 15:32:15 2021

@author: buckman
"""

import os
import numpy as np
import dfm_tools as dfmt
import xarray as xr
import matplotlib.pyplot as plt
plt.close('all')

dir_output = '.' #r'p:\11206811-002-kpp-veerse-meer\model\initial_conditions\VM_WQ_3D_run9c'

if not os.path.exists(dir_output):
    os.makedirs(dir_output)

file_nc = r'p:\11206811-002-kpp-veerse-meer\model\runs_2011-2012\VM_WQ_3D_run9_c\DFM_OUTPUT_VM_WQ_3D\VM_WQ_3D_0000_20130101_000000_rst.nc' #mf1_rstfile (without topology var)

ds = xr.open_dataset(file_nc)
ds = ds.set_coords(['FlowElem_xzw','FlowElem_yzw'])
#uds = dfmt.open_partitioned_dataset(file_nc) #TODO: fails since no grids are present (no variable with cf_role:mesh_topology attr, which can be reconstructed but also no node_coordinates present in dataset)

# define map variables for VM
subs = ['DetCS1','DetNS1','DetPS1','DetSiS1']

for sub in subs:
    #var_coords = ds[sub].attrs['coordinates'].split()
    #ds = ds.set_coords(var_coords) #TODO: this does not work since FlowElem_xcc/FlowElem_ycc are not in file. using FlowElem_xzw/FlowElem_yzw instead
    
    #get data to plot
    data_frommap = ds[sub].isel(time=0)
    xyz = data_frommap.to_dataframe().drop('time',axis=1).values
    
    np.savetxt(os.path.join(sub+'.xyz'),xyz,delimiter='\t')
