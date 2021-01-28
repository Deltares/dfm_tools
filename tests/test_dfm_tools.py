#!/usr/bin/env python

"""Tests for dfm_tools package environment"""

import pytest
#import inspect
import os

if 'TEAMCITY_VERSION' in os.environ.keys(): #teamcity path
    dir_testinput = r'\\dfs-trusted.directory.intra\dfs\Teamcity\Testdata\dfm_tools'
else: #default to this path
    dir_testinput = os.path.join(r'c:/DATA','dfm_tools_testdata')
#from dfm_tools.testutils import getmakeoutputdir


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

        
