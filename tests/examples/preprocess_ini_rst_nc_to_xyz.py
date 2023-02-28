# -*- coding: utf-8 -*-
"""
Created on Tue Jun 15 15:32:15 2021

@author: buckman
"""

import os
import numpy as np
import dfm_tools as dfmt
import xarray as xr
import xugrid as xu
import matplotlib.pyplot as plt
plt.close('all')

dir_output = '.' #r'p:\11206811-002-kpp-veerse-meer\model\initial_conditions\VM_WQ_3D_run9c'

if not os.path.exists(dir_output):
    os.makedirs(dir_output)


file_nc = r'p:\11206811-002-kpp-veerse-meer\model\runs_2011-2012\VM_WQ_3D_run9_c\DFM_OUTPUT_VM_WQ_3D\VM_WQ_3D_0000_20130101_000000_rst.nc' #mf1_rstfile (without topology var)

ds = xr.open_dataset(file_nc)
uds = dfmt.open_partitioned_dataset(file_nc) #TODO: fails since no grids are present (no variable with cf_role:mesh_topology attr, which can be reconstructed but also no node_coordinates present in dataset)

# define map variables for VM
subs = ['DetCS1','DetNS1','DetPS1','DetSiS1']
x_coords = dfmt.get_ncmodeldata(file_nc=file_nc, varname='FlowElem_xzw', multipart=True)
y_coords = dfmt.get_ncmodeldata(file_nc=file_nc, varname='FlowElem_yzw', multipart=True)

#data_frommap_merged = dfmt.open_partitioned_dataset(file_nc.replace('_0000_','_0*_')) #TODO: make starred default, but not supported by older code #TODO: 0000 not at end of filename and no mesh2d variable in the file

for sub in subs:
    #get data to plot
    data_frommap = dfmt.get_ncmodeldata(file_nc=file_nc, varname=sub, timestep=0, layer=None, multipart=True)
    xyz = np.stack([x_coords,y_coords,data_frommap.transpose().reshape(-1)])
    
    xyz_out = xyz.transpose()
    np.savetxt(os.path.join(sub+'.xyz'),xyz_out,delimiter='\t')
