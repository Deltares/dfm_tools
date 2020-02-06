#!/usr/bin/env python

"""Tests for `dfm_tools` package."""

import pytest


from dfm_tools.grid import UGrid


file_map = r'n:\My Documents\werkmap\vanJulien_shortmodelfiles\DFM_3D_z_Grevelingen\computations\run01\DFM_OUTPUT_Grevelingen-FM\Grevelingen-FM_0000_map.nc'
file_net = r'n:\My Documents\werkmap\vanJulien_shortmodelfiles\DFM_3D_z_Grevelingen\computations\run01\Grevelingen_FM_grid_20190603_net.nc'

data_nc = UGrid.fromfile(file_net)


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
