#!/usr/bin/env python

"""Tests for dfm_tools package environment"""

import pytest
#import inspect
import os

from dfm_tools.testutils import getmakeoutputdir, gettestinputdir
dir_testinput = gettestinputdir()



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

        
