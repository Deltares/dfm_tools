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
                                      pytest.param(os.path.join(dir_testinput,r'world.ldb'), id='world'),
                                      pytest.param(os.path.join(dir_testinput,r'Maeslant.tek'), id='Maeslant')])
@pytest.mark.unittest    
def test_readpolygon(file_pol):
    """
    this test tests if a netcdf varname can be retrieved from the 'dictionary' and if the variable can be retrieved from de netcdf
    """
    
    dir_output = getmakeoutputdir(function_name=inspect.currentframe().f_code.co_name)
    #dir_output = dir_testoutput
    
    import matplotlib.pyplot as plt
    plt.close('all')
    
    from dfm_tools.polygon import Polygon
    
    """
    file_pol = os.path.join(dir_testinput,r'DFM_3D_z_Grevelingen\geometry\structures\Grevelingen-FM_BL_fxw.pliz')
    file_pol = os.path.join(dir_testinput,r'DFM_3D_z_Grevelingen\geometry\structures\Grevelingen-FM_BL_fxw.pli')
    file_pol = os.path.join(dir_testinput,r'world.ldb')
    file_pol = r'p:\11205258-006-kpp2020_rmm-g6\C_Work\08_RMM_FMmodel\geometry_j13_6-w3\rmm_v1p3_fixed_weirs.pliz'
    file_pol = r'p:\11205258-006-kpp2020_rmm-g6\C_Work\08_RMM_FMmodel\geometry_j13_6-w3\structures\rmm_v1p3_structures.pli'
    file_pol = r'p:\11205258-006-kpp2020_rmm-g6\C_Work\04_randvoorwaarden\keringen\Maeslantkering\Maeslant.tek'
    """
    
    pol_data_list, pol_name_list, pol_comment_list = Polygon.fromfile(file_pol, pd_output=False)
    pol_data_pd_list = Polygon.fromfile(file_pol, pd_output=True)
    
    fig, ax = plt.subplots()
    for iP, pol_data in enumerate(pol_data_list):
        if pol_data.shape[1]>2:
            ax.plot(pol_data_pd_list[iP].loc[:,'datetime'],pol_data[:,2:],'-',linewidth=0.5)
        else:
            ax.plot(pol_data[:,0],pol_data[:,1],'-',linewidth=0.5)
    plt.savefig(os.path.join(dir_output,os.path.basename(file_pol).split('.')[0]))



# region // Helpers

def getmakeoutputdir(function_name):
    dir_output = os.path.join(dir_testoutput,function_name)
    if not os.path.exists(dir_output):
        os.mkdir(dir_output)
    return dir_output

# endregion


