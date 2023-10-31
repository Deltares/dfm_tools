# -*- coding: utf-8 -*-
"""
Created on Tue Jun 15 15:32:15 2021

@author: buckman
"""

import os
import numpy as np
import dfm_tools as dfmt
import matplotlib.pyplot as plt
plt.close('all')

# mf1_rstfile (without topology var)
file_nc = r'p:\archivedprojects\11206811-002-kpp-veerse-meer\model\runs_2011-2012\VM_WQ_3D_run9_c\DFM_OUTPUT_VM_WQ_3D\VM_WQ_3D_0000_20130101_000000_rst.nc'

uds = dfmt.open_partitioned_dataset(file_nc, preprocess=dfmt.enrich_rst_with_map)

# define map variables for VM
subs = ['DetCS1','DetNS1','DetPS1','DetSiS1']

for sub in subs:
    # uds[sub].isel(time=0).ugrid.plot()
    
    #get data to plot
    data_frommap = uds[sub].isel(time=0)
    xyz = data_frommap.to_dataframe().drop('time',axis=1).values
    
    np.savetxt(os.path.join(sub+'.xyz'),xyz,delimiter='\t')
