# -*- coding: utf-8 -*-
"""
Created on Sat Dec  9 17:46:57 2023

@author: veenstra
"""

import os
import pytest
import glob
import dfm_tools as dfmt
from dfm_tools.observations import (ssc_sscid_from_otherid,
                                    ssc_ssh_subset_groups,
                                    )

source_list = ["uhslc-fast", "uhslc-rqds", "psmsl-gnssir", "ssc", "ioc", "rwsddl", 
               "cmems", "cmems-nrt"] # cmems requires credentials
if os.path.exists(r"p:\1230882-emodnet_hrsm\data\GESLA3"):
    # not possible without p-drive connection
    source_list += ["gesla3"]


@pytest.mark.requiressecrets
@pytest.mark.unittest
@pytest.mark.parametrize("source", source_list)
def test_ssh_catalog_subset_expected_fields(source):
    fields_expected = ["geometry", "source", "country", "station_name_unique"]
    
    ssc_catalog_gpd = dfmt.ssh_catalog_subset(source=source)
    for field in fields_expected:
        assert field in ssc_catalog_gpd.columns
    if source not in ["ssc", "psmsl-gnssir", "rwsddl"]:
        assert "time_ndays" in ssc_catalog_gpd.columns


@pytest.mark.requiressecrets
@pytest.mark.unittest
@pytest.mark.parametrize("source", source_list)
def test_ssh_catalog_subset(source):
    lon_min, lon_max, lat_min, lat_max = -6, 5, 48, 50.5 # france
    # lon_min, lon_max, lat_min, lat_max = 123, 148, 23, 47 # japan
    # lon_min, lon_max, lat_min, lat_max = -20, 40, 25, 72
    # time_min, time_max = '2016-01-01','2016-06-01'
    time_min, time_max = '2020-01-01','2020-06-01'
    
    source_list_notime = ["ssc","rwsddl"]
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


@pytest.mark.requiressecrets
@pytest.mark.unittest
@pytest.mark.parametrize("source", source_list)
def test_ssh_retrieve_data(source, tmp_path):
    # ssc does not contain data, only station locations, so early return
    if source=="ssc":
        return
    
    if source=="uhslc-rqds":
        time_min, time_max = '2018-01-01','2018-02-01'
    else:
        time_min, time_max = '2020-01-01','2020-02-01'
    
    ssc_catalog_gpd = dfmt.ssh_catalog_subset(source=source)
    
    index_dict = {"uhslc-fast":0, "uhslc-rqds":2, 
                  "psmsl-gnssir":0, "ioc":0, "rwsddl":6, 
                  "cmems":0, "cmems-nrt":0, # cmems requires credentials
                  "gesla3":0}
    index = index_dict[source]
    ssc_catalog_gpd_sel = ssc_catalog_gpd.iloc[index:index+1]
    
    dfmt.ssh_retrieve_data(ssc_catalog_gpd_sel, dir_output=tmp_path, 
                           time_min=time_min, time_max=time_max)
    nc_list = glob.glob(os.path.join(tmp_path, "*.nc"))
    assert len(nc_list)==1


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

