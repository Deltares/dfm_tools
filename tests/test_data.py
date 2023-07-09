# -*- coding: utf-8 -*-
"""
Created on Sun Jul  9 22:17:16 2023

@author: veenstra
"""

import pytest
import xugrid as xu
import dfm_tools as dfmt


@pytest.mark.unittest
def test_fm_curvedbend_map():
    uds = dfmt.data.fm_curvedbend_map()
    assert isinstance(uds,xu.UgridDataset)


@pytest.mark.unittest
def test_gshhs_coastlines_shp():
    gshhs_dir = dfmt.data.gshhs_coastlines_shp()
    assert isinstance(gshhs_dir,str)
