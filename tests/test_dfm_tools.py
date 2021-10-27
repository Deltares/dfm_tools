#!/usr/bin/env python

"""Tests for dfm_tools package environment"""

import pytest
#import inspect
import os
import glob

from dfm_tools.testutils import getmakeoutputdir, gettestinputdir
dir_testinput = os.path.join(r'c:\DATA','dfm_tools_testdata')




modulename_list = ['os','sys','glob','shutil','scipy','numpy','datetime','pandas','matplotlib','netCDF4','click','shapely','shapely.geometry','cartopy','pyepsg']
@pytest.mark.parametrize("modulename", [pytest.param('%s'%(stat), id='%s'%(stat)) for stat in modulename_list])
@pytest.mark.unittest
def test_import_libraries(modulename):
    """
    tests whether shapely can be imported successfully, this is a problem in some environments
    in that case 'import shapely' works, but import 'shapely.geometry' fails
    """
    from dfm_tools.testutils import try_importmodule
    
    try_importmodule(modulename=modulename)



dir_scriptfile = os.path.realpath(__file__) #F9 doesnt work, only F5 (F5 also only method to reload external definition scripts)
dir_tests = os.path.abspath(os.path.join(dir_scriptfile,os.pardir))  #1 level up from dir_scripts

dir_output_general = os.path.join(dir_tests,'examples_output')
if not os.path.exists(dir_output_general):
    os.mkdir(dir_output_general)
#os.chdir(dir_output)


# High level acceptance tests, these are the ones who are only meant to generate output files
# for the testers to verify (in Teamcity) whether the runs generate the expected files or not.
""" Run hatyan_main.py with test-configfiles as input """
list_configfiles = glob.glob(os.path.join(dir_tests,'examples','*.py'))
#list_configfiles = ['predictie_2019_b02ex2_19Ycomp4Ydia_CUXHVN_test.py']


@pytest.mark.acceptance
@pytest.mark.parametrize("file_config", [pytest.param(file_config, id=os.path.basename(file_config).replace('.py','')) for file_config in list_configfiles])
def test_run_examples(file_config):
    """
    file_config = os.path.join(dir_tests,'configfiles','predictie_2019_b02ex2_19Ycomp4Ydia_CUXHVN_test.py')
    """
    # 1. Set up test data
    dir_output = os.path.join(dir_output_general,os.path.basename(file_config).replace('.py',''))
    if not os.path.exists(dir_output):
        os.mkdir(dir_output)
    os.chdir(dir_output)
    #test = os.system("python {0} > {1}/FILE_DIAGNOSTICS.txt 2>&1".format(file_config, dir_output))#+ " & pause")
    test = os.system("python {0}".format(file_config))#+ " & pause")
    
    if test:
        raise Exception('execution did not finish properly')


