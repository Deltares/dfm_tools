# -*- coding: utf-8 -*-
"""
Created on Sat Dec  9 17:46:57 2023

@author: veenstra
"""

import os
import pytest
import dfm_tools as dfmt
from dfm_tools.observations import (ssc_sscid_from_otherid,
                                    ssc_ssh_subset_groups,
                                    )


@pytest.mark.unittest
def test_ssh_catalog_subset_expected_fields():
    fields_expected = ["geometry", "source", "country", "station_name_unique"]
    source_list = ["uhslc-fast", "uhslc-rqds", "psmsl-gnssir", "ssc", "ioc"]
    if os.path.exists(r"p:\1230882-emodnet_hrsm\data\GESLA3"):
        # not possible without p-drive connection
        source_list += ["gesla3"]
    if os.name=="nt":
        source_list += ["cmems"] # TODO: not possible on Github, due to missing credentials
    for source in source_list:
        ssc_catalog_gpd = dfmt.ssh_catalog_subset(source=source)
        for field in fields_expected:
            assert field in ssc_catalog_gpd.columns
        if source not in ["ssc", "psmsl-gnssir"]:
            assert "time_ndays" in ssc_catalog_gpd.columns


@pytest.mark.unittest
def test_ssh_catalog_subset():
    lon_min, lon_max, lat_min, lat_max = -6, 5, 48, 50.5 # france
    # lon_min, lon_max, lat_min, lat_max = 123, 148, 23, 47 # japan
    # lon_min, lon_max, lat_min, lat_max = -20, 40, 25, 72
    # time_min, time_max = '2016-01-01','2016-06-01'
    time_min, time_max = '2020-01-01','2020-06-01'
    
    source_list_witime = ["uhslc-fast", "uhslc-rqds", "psmsl-gnssir", "ioc"]
    if os.path.exists(r"p:\1230882-emodnet_hrsm\data\GESLA3"):
        # not possible without p-drive connection
        source_list_witime += ["gesla3"]
    if os.name=="nt":
        source_list_witime += ["cmems"] # TODO: not possible on Github, due to missing credentials
    source_list_notime = ["ssc"]
    for source in source_list_witime+source_list_notime:
        ssc_catalog_gpd = dfmt.ssh_catalog_subset(source=source)
        if source in source_list_notime:
            ssc_catalog_gpd_sel = dfmt.ssh_catalog_subset(source=source,
                                                          lon_min=lon_min, lon_max=lon_max, 
                                                          lat_min=lat_min, lat_max=lat_max)
        else:
            ssc_catalog_gpd_sel = dfmt.ssh_catalog_subset(source=source,
                                                          lon_min=lon_min, lon_max=lon_max, 
                                                          lat_min=lat_min, lat_max=lat_max, 
                                                          time_min=time_min, time_max=time_max)
        assert len(ssc_catalog_gpd) > len(ssc_catalog_gpd_sel)


@pytest.mark.unittest
def test_ssh_retrieve_data(tmp_path):
    time_min, time_max = '2020-01-01','2020-02-01'
    
    source_list = ["ioc", "uhslc-fast", "uhslc-rqds", "psmsl-gnssir"]
    if os.path.exists(r"p:\1230882-emodnet_hrsm"):
        # not possible without p-drive connection
        source_list += ["gesla3"]
    if os.name=="nt":
        source_list += ["cmems"] # TODO: not possible on Github, due to missing credentials
    for source in source_list:
        ssc_catalog_gpd = dfmt.ssh_catalog_subset(source=source)
        ssc_catalog_gpd_sel = ssc_catalog_gpd.iloc[:1]
        if source=="cmems": #TODO: remove this exception when the cmems API works for insitu data
            dfmt.ssh_retrieve_data(ssc_catalog_gpd_sel, dir_output=tmp_path)
        else:
            dfmt.ssh_retrieve_data(ssc_catalog_gpd_sel, dir_output=tmp_path, 
                                   time_min=time_min, time_max=time_max)


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


@pytest.mark.unittest
def test_ssh_catalog_toxynfile(tmp_path):
    ssc_catalog_gpd = dfmt.ssh_catalog_subset(source="ssc")
    file_xyn = tmp_path / 'test_ssc_obs.xyn'
    dfmt.ssh_catalog_toxynfile(ssc_catalog_gpd, file_xyn)
    assert os.path.isfile(file_xyn)

