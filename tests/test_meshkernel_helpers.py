# -*- coding: utf-8 -*-
"""
Created on Wed Jul  5 12:59:23 2023

@author: veenstra
"""

import pytest
import xugrid as xu
from dfm_tools.meshkernel_helpers import add_crs_to_dataset


@pytest.mark.unittest
def test_add_crs_to_dataset_cartesian():
    uds = xu.data.adh_san_diego()
    crs='EPSG:26946' # this is not the correct crs for this model, but that does not matter
    add_crs_to_dataset(uds,crs=crs,is_geographic=False)
    
    assert 'projected_coordinate_system' in uds.data_vars
    crs_attrs = uds.projected_coordinate_system.attrs
    assert crs_attrs['EPSG_code'] == 'EPSG:26946'
    assert crs_attrs['epsg'] == 26946
    assert crs_attrs['grid_mapping_name'] == 'Unknown projected'

def test_add_crs_to_dataset_spherical():
    uds = xu.data.adh_san_diego()
    crs='EPSG:4326' # this is not the correct crs for this model, but that does not matter
    add_crs_to_dataset(uds,crs=crs,is_geographic=True)
    
    assert 'wgs84' in uds.data_vars
    crs_attrs = uds.wgs84.attrs
    assert crs_attrs['EPSG_code'] == 'EPSG:4326'
    assert crs_attrs['epsg'] == 4326
    assert crs_attrs['grid_mapping_name'] == 'latitude_longitude'