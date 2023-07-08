# -*- coding: utf-8 -*-
"""
Created on Sat Jul  8 11:04:27 2023

@author: veenstra
"""

import pytest
import os
import glob

# ACCEPTANCE TESTS VIA EXAMPLE SCRIPTS, these are the ones who are only meant to generate output files

dir_tests = os.path.dirname(__file__) #F9 doesnt work, only F5 (F5 also only method to reload external definition scripts)
list_configfiles = glob.glob(os.path.join(dir_tests,'examples','*.py')) + glob.glob(os.path.join(dir_tests,'examples_workinprogress','*.py'))
list_configfiles = [x for x in list_configfiles if 'workinprogress_xarray_performance' not in x] #ignore this slow script
dir_output_general = os.path.join(dir_tests,'examples_output')
os.makedirs(dir_output_general, exist_ok=True)

@pytest.mark.requireslocaldata
@pytest.mark.acceptance
@pytest.mark.parametrize("file_config", [pytest.param(file_config, id=os.path.basename(file_config).replace('.py','')) for file_config in list_configfiles])
def test_run_examples(file_config):
    # 1. Set up test data
    dir_output = os.path.join(dir_output_general,os.path.basename(file_config).replace('.py',''))
    if not os.path.exists(dir_output):
        os.mkdir(dir_output)
    os.chdir(dir_output)
    test = os.system(f'python {file_config}')#+ " & pause")
    
    if test:
        raise OSError('execution did not finish properly')