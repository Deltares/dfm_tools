# -*- coding: utf-8 -*-
"""
Created on Wed Mar 18 16:33:52 2020

@author: veenstra
"""

import pytest
import inspect
import os

dir_testinput = os.path.join(r'c:/DATA','dfm_tools_testdata')
from tests.TestUtils import getmakeoutputdir


@pytest.mark.unittest    
def test_delft3d4():
    import matplotlib.pyplot as plt
    plt.close('all')
    #import numpy as np
    
    #from dfm_tools.get_nc import plot_netmapdata
    from dfm_tools.delft3d4 import grid
    from dfm_tools.delft3d4 import dep
    
    dir_output = getmakeoutputdir(__file__,inspect.currentframe().f_code.co_name)
    #dir_output = '.'

    #file_d3d_grid = r'c:\DATA\werkmap\dfm_tools_testdata\D3D_3D_sigma_curved_bend\bendcurv.grd'
    file_d3d_grid = os.path.join(dir_testinput, 'brazil_patos_lagoon_52S_32E', 'lake_and_sea_5_xy.grd')
    file_d3d_dep = os.path.join(dir_testinput, 'brazil_patos_lagoon_52S_32E', 'dep_at_cor_triangulated_filled_corners.dep')
    #file_d3d_mdf = os.path.join(dir_testinput, 'brazil_patos_lagoon_52S_32E', '3d1.mdf')
    #file_d3d_mdf = os.path.join(dir_testinput, 'D3D_3D_sigma_curved_bend','cb2-sal-added-3d.mdf')
    
    data_grd = grid.Grid.fromfile(file_d3d_grid)
    grd_x = data_grd.x
    grd_y = data_grd.y
    grd_shape = data_grd.shape
    data_dep = dep.Dep.read(file_d3d_dep,grd_shape)
    
    fig, ax = plt.subplots(1,1)
    ax.plot(data_grd.x.transpose(), data_grd.y.transpose(), 'g')
    ax.plot(data_grd.x, data_grd.y, 'g')
    ax.scatter(grd_x,grd_y,5,c='b')
    ax.set_aspect('equal')
    plt.savefig(os.path.join(dir_output,'d3d_grd'))
    
    fig, ax = plt.subplots(1,1)
    ax.pcolor(data_grd.x,data_grd.y,data_dep.val[0:-1,0:-1])
    ax.set_aspect('equal')
    plt.savefig(os.path.join(dir_output,'d3d_dep'))
    
    #data_mdf = mdf.read(file_d3d_mdf)


