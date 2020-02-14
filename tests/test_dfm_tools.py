#!/usr/bin/env python

"""Tests for `dfm_tools` package environment"""

import pytest
#import inspect
import os

dir_tests = os.path.join(os.path.realpath(__file__), os.pardir)
dir_testoutput = os.path.join(dir_tests,'test_output')
if not os.path.exists(dir_testoutput):
    os.mkdir(dir_testoutput)
dir_testinput = os.path.join(r'c:/DATA/werkmap','dfm_tools_testdata')


@pytest.mark.unittest
def test_import_shapely():
    """
    tests whether shapely can be imported successfully, this is a problem in some environments
    in that case 'import shapely' works, but import 'shapely.geometry' fails
    """
    
    try:
        import shapely
        import_success1 = True
        print(shapely)
    except:
        print('"import shapely" failed')
        import_success1 = False

    try:
        import shapely.geometry
        import_success2 = True
        print(shapely.geometry)
    except:
        print('"shapely.geometry" failed')
        import_success2 = False

    try:
        from shapely.geometry import Point
        import_success3 = True
        print(Point)
    except:
        print('"from shapely.geometry import Point" failed')
        import_success3 = False
    
    assert import_success1 == True
    assert import_success2 == True
    assert import_success3 == True

        
