# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 19:27:21 2020

@author: veenstra
"""

import pytest
import inspect
import os

dir_tests = os.path.join(os.path.realpath(__file__), os.pardir)
dir_testoutput = os.path.join(dir_tests,'test_output')
if not os.path.exists(dir_testoutput):
    os.mkdir(dir_testoutput)
dir_testinput = os.path.join(r'c:/DATA','dfm_tools_testdata')

@pytest.mark.parametrize("file_pol", [pytest.param(os.path.join(dir_testinput,r'DFM_3D_z_Grevelingen\geometry\structures\Grevelingen-FM_BL_fxw.pliz'), id='Grevelingen pliz'),
                                      pytest.param(os.path.join(dir_testinput,r'DFM_3D_z_Grevelingen\geometry\structures\Grevelingen-FM_BL_fxw.pli'), id='Grevelingen pli'),
                                      pytest.param(os.path.join(dir_testinput,r'world.ldb'), id='world')])
@pytest.mark.unittest    
def test_readpolygon(file_pol):
    """
    this test tests if a netcdf varname can be retrieved from the 'dictionary' and if the variable can be retrieved from de netcdf
    """
    
    import matplotlib.pyplot as plt
    plt.close('all')
    
    from dfm_tools.polygon import Polygon
    
    """
    file_pol = os.path.join(dir_testinput,r'DFM_3D_z_Grevelingen\geometry\structures\Grevelingen-FM_BL_fxw.pliz')
    file_pol = os.path.join(dir_testinput,r'DFM_3D_z_Grevelingen\geometry\structures\Grevelingen-FM_BL_fxw.pli')
    file_pol = os.path.join(dir_testinput,r'world.ldb')
    file_pol = r'p:\11205258-006-kpp2020_rmm-g6\C_Work\08_RMM_FMmodel\geometry_j13_6-w3\rmm_v1p3_fixed_weirs.pliz'
    file_pol = r'p:\11205258-006-kpp2020_rmm-g6\C_Work\08_RMM_FMmodel\geometry_j13_6-w3\structures\rmm_v1p3_structures.pli'
    """
    
    pol_data_list, pol_name_list = Polygon.fromfile(file_pol)
    fig, ax = plt.subplots()
    for iP, pol_data in enumerate(pol_data_list):
        ax.plot(pol_data[:,0],pol_data[:,1],'-',linewidth=0.5)

    
