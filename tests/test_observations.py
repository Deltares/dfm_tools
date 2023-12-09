# -*- coding: utf-8 -*-
"""
Created on Sat Dec  9 17:46:57 2023

@author: veenstra
"""

import pytest
from dfm_tools.observations import (ssh_catalog_subset,
                                    ssc_sscid_from_otherid,
                                    ssc_ssh_subset_groups,
                                    )


@pytest.mark.unittest
def test_ssh_catalog_subset():
    
    lon_min, lon_max, lat_min, lat_max = -6, 5, 48, 50.5 # france
    # lon_min, lon_max, lat_min, lat_max = 123, 148, 23, 47 # japan
    # lon_min, lon_max, lat_min, lat_max = -20, 40, 25, 72
    # time_min, time_max = '2016-01-01','2016-06-01'
    time_min, time_max = '2020-01-01','2020-06-01'
    
    fields_expected = ["geometry", "source"]
    source_list_witime = ["cmems", "uhslc-fast", "uhslc-rqds", "psmsl-gnssir"] # "gesla3" not possible on github
    source_list_notime = ["ssc", "ioc"]
    for source in source_list_witime+source_list_notime:
        ssc_catalog_gpd = ssh_catalog_subset(source=source)
        if source in source_list_notime:
            ssc_catalog_gpd_sel = ssh_catalog_subset(source=source,
                                                     lon_min=lon_min, lon_max=lon_max, 
                                                     lat_min=lat_min, lat_max=lat_max, 
                                                     )
        else:
            ssc_catalog_gpd_sel = ssh_catalog_subset(source=source,
                                                      lon_min=lon_min, lon_max=lon_max, 
                                                      lat_min=lat_min, lat_max=lat_max, 
                                                      time_min=time_min, time_max=time_max,
                                                      )
        assert len(ssc_catalog_gpd) > len(ssc_catalog_gpd_sel)
        for field in fields_expected:
            assert field in ssc_catalog_gpd.columns


@pytest.mark.unittest
def test_ssc_sscid_from_otherid():
    sscid_from_uhslcid = ssc_sscid_from_otherid(group_id=347, groupname='uhslc')
    assert sscid_from_uhslcid=="SSC-abas"


@pytest.mark.unittest
def test_ssc_ssh_subset_groups():
    ssc_catalog_gpd_uhslc = ssc_ssh_subset_groups(groups='uhslc')
    ssc_catalog_gpd_ioc = ssc_ssh_subset_groups(groups='ioc')
    ssc_catalog_gpd_twogroups = ssc_ssh_subset_groups(groups=['ioc','uhslc'])
    assert len(ssc_catalog_gpd_uhslc) < len(ssc_catalog_gpd_twogroups)
    assert len(ssc_catalog_gpd_ioc) < len(ssc_catalog_gpd_twogroups)
    assert len(ssc_catalog_gpd_ioc) + len(ssc_catalog_gpd_uhslc) > len(ssc_catalog_gpd_twogroups)

