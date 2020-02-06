# -*- coding: utf-8 -*-
"""
Created on Thu Feb  6 10:00:01 2020

@author: veenstra
"""

#!/usr/bin/env python

"""Tests for dflowutil form dfm_tools package."""

import pytest


from dflowutil import mesh



class Test_dflowutil_mesh:
    @pytest.mark.acceptance
    def mesh_plotncmap(self):
        file_map = r'n:\My Documents\werkmap\vanJulien_shortmodelfiles\DFM_3D_z_Grevelingen\computations\run01\DFM_OUTPUT_Grevelingen-FM\Grevelingen-FM_0000_map.nc'
        file_net = r'n:\My Documents\werkmap\vanJulien_shortmodelfiles\DFM_3D_z_Grevelingen\computations\run01\Grevelingen_FM_grid_20190603_net.nc'
        mesh.plot_nc_map(file_map,elem='waterlevel',time=1,layer=1)
        mesh.plot_net(file_net)
        

@pytest.fixture
def response():
    """Sample pytest fixture.

    See more at: http://doc.pytest.org/en/latest/fixture.html
    """
    # import requests
    # return requests.get('https://github.com/audreyr/cookiecutter-pypackage')


def test_content(response):
    """Sample pytest test function with the pytest fixture as an argument."""
    # from bs4 import BeautifulSoup
    # assert 'GitHub' in BeautifulSoup(response.content).title.string


