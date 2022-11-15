# -*- coding: utf-8 -*-
"""
Created on Tue Jun 15 15:32:15 2021

@author: buckman
"""
import os

import numpy as np
import dfm_tools as dfmt

file_nc = r'p:\11206811-002-kpp-veerse-meer\model\runs_2011-2012\VM_WQ_3D_run9_c\DFM_OUTPUT_VM_WQ_3D\VM_WQ_3D_0000_20130101_000000_rst.nc'
dir_output = r'p:\11206811-002-kpp-veerse-meer\model\initial_conditions\VM_WQ_3D_run9c'

if not os.path.exists(dir_output):
    os.makedirs(dir_output)

# define map variables for VM
times_pd = dfmt.get_timesfromnc(file_nc=file_nc)
vars_pd = dfmt.get_ncvarproperties(file_nc=file_nc)
subs = ['DetCS1','DetNS1','DetPS1','DetSiS1']
x_coords = dfmt.get_ncmodeldata(file_nc=file_nc, varname='FlowElem_xzw', multipart=True)
y_coords = dfmt.get_ncmodeldata(file_nc=file_nc, varname='FlowElem_yzw', multipart=True)

for sub in subs:    
    #get data to plot
    data_frommap = dfmt.get_ncmodeldata(file_nc=file_nc, varname=sub, timestep=0, layer=None, multipart=True)
    xyz = np.stack([x_coords,y_coords,data_frommap.transpose().reshape(-1)])
    
    xyz_out = xyz.transpose()
    np.savetxt(os.path.join(dir_output,sub+'.xyz'),xyz_out,delimiter='\t')
