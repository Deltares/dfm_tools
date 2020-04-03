# -*- coding: utf-8 -*-
"""
Created on Wed Mar 18 16:33:52 2020

@author: veenstra
"""

import pytest
import inspect
import os

dir_testinput = os.path.join(r'c:/DATA','dfm_tools_testdata')
from dfm_tools.testutils import getmakeoutputdir


@pytest.mark.unittest    
def test_delft3d4():
    import matplotlib.pyplot as plt
    plt.close('all')
    #import numpy as np
    
    #from dfm_tools.get_nc import plot_netmapdata
    from dfm_tools.io import grd
    from dfm_tools.io import dep
    #from dfm_tools.io import mdf #does not work yet
    
    dir_output = getmakeoutputdir(__file__,inspect.currentframe().f_code.co_name)
    #dir_output = './test_output'

    #file_d3d_grid = r'c:\DATA\werkmap\dfm_tools_testdata\D3D_3D_sigma_curved_bend\bendcurv.grd'
    file_d3d_grid = os.path.join(dir_testinput, 'brazil_patos_lagoon_52S_32E', 'lake_and_sea_5_xy.grd')
    file_d3d_dep = os.path.join(dir_testinput, 'brazil_patos_lagoon_52S_32E', 'dep_at_cor_triangulated_filled_corners.dep')
    #file_d3d_mdf = os.path.join(dir_testinput, 'brazil_patos_lagoon_52S_32E', '3d1.mdf')
    #file_d3d_mdf = os.path.join(dir_testinput, 'D3D_3D_sigma_curved_bend','cb2-sal-added-3d.mdf')
    
    data_grd = grd.Grid.fromfile(file_d3d_grid)
    grd_shape = data_grd.shape
    data_dep = dep.Dep.read(file_d3d_dep,grd_shape)
    
    fig, ax = plt.subplots(1,1)
    ax.plot(data_grd.x.transpose(), data_grd.y.transpose(), 'g', linewidth=0.5)
    ax.plot(data_grd.x, data_grd.y, 'g', linewidth=0.5)
    #ax.scatter(data_grd.x,data_grd.y,2,c=data_dep.val[0:-1,0:-1])
    ax.set_aspect('equal')
    plt.savefig(os.path.join(dir_output,'d3d_grd'))
    
    fig, ax = plt.subplots(1,1)
    ax.pcolor(data_grd.x,data_grd.y,data_dep.val[0:-1,0:-1])
    ax.set_aspect('equal')
    plt.savefig(os.path.join(dir_output,'d3d_deppcolor'))
    
    #data_mdf = mdf.read(file_d3d_mdf)






@pytest.mark.acceptance
def test_mdu():
    """
    tests whether mdu file can be imported and exported
    """

    dir_output = getmakeoutputdir(__file__,inspect.currentframe().f_code.co_name)
    #dir_output = './test_output'
    if not os.path.exists(dir_output):
        os.mkdir(dir_output)

    from dfm_tools.io import mdu
        
    try:
        filename_mdu = os.path.join(dir_testinput, r'DFM_3D_z_Grevelingen\computations\run01\Grevelingen-FM.mdu')
        data_mdu = mdu.read_deltares_ini(filename_mdu)
        print(data_mdu)
        import_success = True
    except:
        import_success = False

    try:
        filename_mdu_out = os.path.join(dir_output, 'Grevelingen-FM_out.mdu')
        mdu.write_deltares_ini(data_mdu, filename_mdu_out)
        export_success = True
    except:
        export_success = False
    
    assert import_success == True
    assert export_success == True






@pytest.mark.parametrize("file_pol", [pytest.param(os.path.join(dir_testinput,'DFM_3D_z_Grevelingen\\geometry\\structures\\Grevelingen-FM_BL_fxw.pliz'), id='Grevelingen pliz'),
                                      pytest.param(os.path.join(dir_testinput,'DFM_3D_z_Grevelingen\\geometry\\structures\\Grevelingen-FM_BL_fxw.pli'), id='Grevelingen pli'),
                                      pytest.param(os.path.join(dir_testinput,'world.ldb'), id='world'),
                                      pytest.param(os.path.join(dir_testinput,'Maeslant.tek'), id='Maeslant')])
@pytest.mark.unittest    
def test_readpolygon(file_pol):
    """
    this test tests if a netcdf varname can be retrieved from the 'dictionary' and if the variable can be retrieved from de netcdf
    file_pol = os.path.join(dir_testinput,'DFM_3D_z_Grevelingen\\geometry\\structures\\Grevelingen-FM_BL_fxw.pliz')
    file_pol = os.path.join(dir_testinput,'DFM_3D_z_Grevelingen\\geometry\\structures\\Grevelingen-FM_BL_fxw.pli')
    file_pol = os.path.join(dir_testinput,'world.ldb')
    file_pol = 'p:\\11205258-006-kpp2020_rmm-g6\\C_Work\\08_RMM_FMmodel\\geometry_j13_6-w3\\rmm_v1p3_fixed_weirs.pliz'
    file_pol = 'p:\\11205258-006-kpp2020_rmm-g6\\C_Work\\08_RMM_FMmodel\\geometry_j13_6-w3\\structures\\rmm_v1p3_structures.pli'
    file_pol = 'p:\\11205258-006-kpp2020_rmm-g6\\C_Work\\04_randvoorwaarden\\keringen\\Maeslantkering\\Maeslant.tek'
    """
    
    dir_output = getmakeoutputdir(__file__,inspect.currentframe().f_code.co_name)
    #dir_output = './test_output'
    
    import matplotlib.pyplot as plt
    plt.close('all')
    
    from dfm_tools.io.polygon import Polygon
        
    pol_data_list, pol_name_list, pol_comment_list = Polygon.fromfile(file_pol, pd_output=False)
    pol_data_pd_list = Polygon.fromfile(file_pol, pd_output=True)
    
    fig, ax = plt.subplots()
    for iP, pol_data in enumerate(pol_data_list):
        if 'datetime' in pol_data_pd_list[iP].columns.tolist():
            ax.plot(pol_data_pd_list[iP].loc[:,'datetime'],pol_data[:,2:],'-',linewidth=0.5)
        else:
            ax.plot(pol_data[:,0],pol_data[:,1],'-',linewidth=0.5)
    plt.savefig(os.path.join(dir_output,os.path.basename(file_pol).replace('.','')))




