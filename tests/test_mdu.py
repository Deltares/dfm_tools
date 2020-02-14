# -*- coding: utf-8 -*-
"""
Created on Thu Feb 13 23:05:01 2020

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


@pytest.mark.acceptance
def test_mdu():
    """
    tests whether mdu file can be imported and exported
    """

    this_function_name = inspect.currentframe().f_code.co_name
    dir_output = os.path.join(dir_testoutput,this_function_name)
    if not os.path.exists(dir_output):
        os.mkdir(dir_output)

    from dfm_tools.mdu import read_deltares_ini, write_deltares_ini    
        
    try:
        filename_mdu = os.path.join(dir_testinput, r'DFM_3D_z_Grevelingen\computations\run01\Grevelingen-FM.mdu')
        data_mdu = read_deltares_ini(filename_mdu)
        print(data_mdu)
        import_success = True
    except:
        import_success = False

    try:
        filename_mdu_out = os.path.join(dir_output, 'Grevelingen-FM_out.mdu')
        write_deltares_ini(data_mdu, filename_mdu_out)
        export_success = True
    except:
        export_success = False
    
    assert import_success == True
    assert export_success == True

