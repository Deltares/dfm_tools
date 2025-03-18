# -*- coding: utf-8 -*-
"""
Created on Sat Dec  9 17:46:57 2023

@author: veenstra
"""

import os
import pytest
import glob
import ddlpy
import dfm_tools as dfmt
from dfm_tools.observations import (ssc_sscid_from_otherid,
                                    ssc_ssh_subset_groups,
                                    )
from dfm_tools.observations import (gtsm3_era5_cds_ssh_read_catalog,
                                    gtsm3_era5_cds_ssh_retrieve_data,
                                    )

source_list = ["uhslc-fast", "uhslc-rqds", "psmsl-gnssir", "ssc", "ioc", "rwsddl", 
               "cmems", "cmems-nrt", # requires CMEMS credentials
               "gtsm3-era5-cds", # requires CDS credentials
               ] 
if os.path.exists(r"p:\1230882-emodnet_hrsm\data\GESLA3"):
    # not possible without p-drive connection
    source_list += ["gesla3"]


@pytest.mark.requiressecrets
@pytest.mark.unittest
@pytest.mark.parametrize("source", source_list)
def test_ssh_catalog_subset(source):
    lon_min, lon_max, lat_min, lat_max = -6, 5, 48, 50.5 # france
    # lon_min, lon_max, lat_min, lat_max = 123, 148, 23, 47 # japan
    # lon_min, lon_max, lat_min, lat_max = -20, 40, 25, 72
    # time_min, time_max = '2016-01-01','2016-06-01'
    time_min, time_max = '2020-01-01','2020-06-01'
    
    ssc_catalog_gpd = dfmt.ssh_catalog_subset(source=source)
    
    # check the returned fields
    fields_expected = ["geometry", "source", "country", 
                       "station_name", "station_id", "station_name_unique"]
    for field in fields_expected:
        assert field in ssc_catalog_gpd.columns
    if source not in ["ssc", "psmsl-gnssir", "rwsddl"]:
        assert "time_min" in ssc_catalog_gpd.columns
        assert "time_max" in ssc_catalog_gpd.columns
        assert "time_ndays" in ssc_catalog_gpd.columns
    assert ssc_catalog_gpd.crs.to_string()=='EPSG:4326'
    
    # do actual subset
    source_list_notime = ["ssc","rwsddl"]
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
        # 2020 not available in uhslc-rqds yet
        time_min, time_max = '2018-01-01','2018-02-01'
    else:
        time_min, time_max = '2020-01-01','2020-02-01'
    
    ssc_catalog_gpd = dfmt.ssh_catalog_subset(source=source)
    if source=="rwsddl":
        # order of rows in rwsddl locations dataframe is python-version-dependent
        # make sure we always test on the same hist station (no realtime data)
        bool_hoekvhld = ssc_catalog_gpd["Code"].isin(["HOEKVHLD"])
        ssc_catalog_gpd = ssc_catalog_gpd.loc[bool_hoekvhld]
    
    index_dict = {"uhslc-fast":0, "uhslc-rqds":2, 
                  "psmsl-gnssir":0, "ioc":0, "rwsddl":0, 
                  "cmems":0, "cmems-nrt":0, # requires CMEMS credentials
                  "gtsm3-era5-cds":0, # requires CDS credentials
                  "gesla3":0,
                  }
    index = index_dict[source]
    ssc_catalog_gpd_sel = ssc_catalog_gpd.iloc[index:index+1]
    
    dfmt.ssh_retrieve_data(ssc_catalog_gpd_sel, dir_output=tmp_path, 
                           time_min=time_min, time_max=time_max)
    nc_list = glob.glob(os.path.join(tmp_path, "*.nc"))
    assert len(nc_list)==1


@pytest.mark.unittest
def test_ssh_netcdf_overview(tmp_path):
    source = "rwsddl"
    time_min, time_max = '2020-01-01','2020-01-05'
    
    ssc_catalog_gpd = dfmt.ssh_catalog_subset(source=source)
    
    # order of rows in rwsddl locations dataframe is python-version-dependent
    # make sure we always test on the same hist station (no realtime data)
    bool_hoekvhld = ssc_catalog_gpd["Code"].isin(["HOEKVHLD"])
    ssc_catalog_gpd_sel = ssc_catalog_gpd.loc[bool_hoekvhld]
    
    dfmt.ssh_retrieve_data(ssc_catalog_gpd_sel, dir_output=tmp_path, 
                           time_min=time_min, time_max=time_max)
    
    dfmt.ssh_netcdf_overview(tmp_path)
    assert os.path.isdir(os.path.join(tmp_path, "overview"))
    assert os.path.isfile(os.path.join(tmp_path, "overview", "overview_availability_001_001.png"))
    assert os.path.isfile(os.path.join(tmp_path, "overview", "waterlevel_data_netcdf_overview.csv"))


@pytest.mark.unittest
def test_rwsddl_ssh_get_time_max():
    locations = ddlpy.locations()
    bool_hoedanigheid = locations['Hoedanigheid.Code'].isin(['NAP'])
    bool_stations = locations.index.isin(['HOEKVHLD', 'IJMDBTHVN','SCHEVNGN'])
    bool_grootheid = locations['Grootheid.Code'].isin(['WATHTE'])
    bool_groepering = locations['Groepering.Code'].isin(['NVT'])
    selected = locations.loc[bool_grootheid & bool_hoedanigheid & bool_groepering & bool_stations]
    selected_withtimemax = dfmt.observations.rwsddl_ssh_get_time_max(selected)
    assert "time_max" not in selected.columns
    assert "time_max" in selected_withtimemax.columns

    
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


@pytest.mark.unittest
def test_ssh_catalog_tokmlfile(tmp_path):
    ssc_catalog_gpd = dfmt.ssh_catalog_subset(source="ssc")
    file_xyn = tmp_path / 'test_ssc.kml'
    dfmt.ssh_catalog_tokmlfile(ssc_catalog_gpd, file_xyn)
    assert os.path.isfile(file_xyn)


@pytest.mark.unittest
def test_gtsm3_era5_cds_ssh_retrieve_data_invalidfreq_nonetimes():
    df = gtsm3_era5_cds_ssh_read_catalog()
    
    with pytest.raises(ValueError) as e:
        gtsm3_era5_cds_ssh_retrieve_data(
            row=df.iloc[0],
            dir_output=".",
            time_min=None,
            time_max=None,
            time_freq='10min',
            )
    assert "time frequency for retrieving gtsm3-era5-cds data should" in str(e.value)
