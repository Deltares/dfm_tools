# -*- coding: utf-8 -*-
"""
Created on Wed Jul  5 12:59:23 2023

@author: veenstra
"""

import pytest
import xugrid as xu
from dfm_tools.meshkernel_helpers import add_crs_to_dataset, make_basegrid, meshkernel_check_geographic


@pytest.mark.unittest
def test_add_crs_to_dataset_cartesian():
    uds = xu.data.adh_san_diego()
    crs='EPSG:26946' # this is not the correct crs for this model, but that does not matter
    add_crs_to_dataset(uds,is_geographic=False,crs=crs)
    
    assert 'projected_coordinate_system' in uds.data_vars
    crs_attrs = uds.projected_coordinate_system.attrs
    assert crs_attrs['EPSG_code'] == 'EPSG:26946'
    assert crs_attrs['epsg'] == 26946
    assert crs_attrs['grid_mapping_name'] == 'Unknown projected'


@pytest.mark.unittest
def test_add_crs_to_dataset_spherical():
    uds = xu.data.adh_san_diego()
    crs='EPSG:4326' # this is not the correct crs for this model, but that does not matter
    add_crs_to_dataset(uds,is_geographic=True,crs=crs)
    
    assert 'wgs84' in uds.data_vars
    crs_attrs = uds.wgs84.attrs
    assert crs_attrs['EPSG_code'] == 'EPSG:4326'
    assert crs_attrs['epsg'] == 4326
    assert crs_attrs['grid_mapping_name'] == 'latitude_longitude'
    

@pytest.mark.systemtest
def test_meshkernel_check_geographic():
    """
    to check whether is_geographic can be correctly derived from the mk object, for cartesian as well spherical objects
    """
    lon_min, lon_max, lat_min, lat_max = -68.55, -67.9, 11.8, 12.6
    dxy = 0.05
    
    mk_cartesian = make_basegrid(lon_min, lon_max, lat_min, lat_max, dx=dxy, dy=dxy, is_geographic=False)
    mk_cartesian_geograph = meshkernel_check_geographic(mk_cartesian)
    mk_spherical = make_basegrid(lon_min, lon_max, lat_min, lat_max, dx=dxy, dy=dxy, is_geographic=True)
    mk_spherical_geograph = meshkernel_check_geographic(mk_spherical)
    
    assert mk_cartesian_geograph==False
    assert mk_spherical_geograph==True

