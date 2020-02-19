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
dir_testinput = os.path.join(r'c:/DATA/werkmap','dfm_tools_testdata')

@pytest.mark.parametrize("file_pol", [pytest.param(os.path.join(dir_testinput,r'DFM_3D_z_Grevelingen\geometry\structures\Grevelingen-FM_BL_fxw.pliz'), id='Grevelingen'),
                                      pytest.param(os.path.join(dir_testinput,r'world.ldb'), id='world')])
@pytest.mark.unittest    
def test_readpolygon(file_pol):
    """
    this test tests if a netcdf varname can be retrieved from the 'dictionary' and if the variable can be retrieved from de netcdf
    """
    
    from dfm_tools.polygon import Polygon
    
    #file_pol = os.path.join(dir_testinput,r'DFM_3D_z_Grevelingen\geometry\structures\Grevelingen-FM_BL_fxw.pli')
    #file_pol = os.path.join(dir_testinput,r'world.ldb')
    
    data_pol = Polygon.fromfile(file_pol)
    data_pol_x = data_pol.x
    data_pol_y = data_pol.y
    data_pol_names = data_pol.name
    data_pol_linearray = data_pol.line_array

    
    