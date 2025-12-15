# -*- coding: utf-8 -*-
"""
Created on Sat Dec  9 17:46:57 2023

@author: veenstra
"""

import os
import pytest
import glob
import ddlpy
import numpy as np
import dfm_tools as dfmt
from dfm_tools.observations import (ssc_ssh_read_catalog,
                                    ssc_sscid_from_otherid,
                                    ssc_ssh_subset_groups,
                                    gtsm3_era5_cds_ssh_read_catalog,
                                    gtsm3_era5_cds_ssh_retrieve_data,
                                    _remove_accents,
                                    )
import logging

source_list = ["uhslc", "psmsl-gnssir", "ssc", "ioc", "rwsddl", 
               "cmems", "cmems-nrt", # requires CMEMS credentials
               "gtsm3-era5-cds", # requires CDS credentials
               ] 
if os.path.exists(r"p:\metocean-data\licensed\GESLA3"):
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
    
    # check if station names/ids can be converted to S64 as done in
    # dfm_tools.observations._make_hydrotools_consistent()
    # fixed for uhslc in https://github.com/Deltares/dfm_tools/issues/1172
    # ssc still has accents but does not have data so .astype("S64") will not be called
    if source not in ["ssc"]:
        ssc_catalog_gpd["station_name"].astype("S64")
        ssc_catalog_gpd["station_id"].astype("S64")
        ssc_catalog_gpd["station_name_unique"].astype("S64")


@pytest.mark.requiressecrets
@pytest.mark.unittest
@pytest.mark.parametrize("source", source_list)
def test_ssh_retrieve_data(source, tmp_path):
    # ssc does not contain data, only station locations, so early return
    if source=="ssc":
        return
    
    time_min, time_max = '2020-01-01','2020-02-01'
    
    ssc_catalog_gpd = dfmt.ssh_catalog_subset(source=source)
    if source=="rwsddl":
        # order of rows in rwsddl locations dataframe is python-version-dependent
        # make sure we always test on the same hist station (no realtime data)
        bool_hoekvhld = ssc_catalog_gpd["Code"].isin(["hoekvanholland"])
        ssc_catalog_gpd = ssc_catalog_gpd.loc[bool_hoekvhld]
    
    if source == "ioc":
        # stat_index=0 fails since 2025-08-15 since AMTSI was added
        # this station has only data from that date onwards
        stat_index = 2
    else:
        stat_index = 0
    ssc_catalog_gpd_sel = ssc_catalog_gpd.iloc[[stat_index]]
    
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
    bool_hoekvhld = ssc_catalog_gpd["Code"].isin(["hoekvanholland"])
    ssc_catalog_gpd_sel = ssc_catalog_gpd.loc[bool_hoekvhld]
    assert len(ssc_catalog_gpd_sel) == 1
    
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
    bool_stations = locations.index.isin(['hoekvanholland', 'ijmuiden.buitenhaven','scheveningen'])
    bool_procestype = locations['ProcesType'].isin(['meting'])
    bool_grootheid = locations['Grootheid.Code'].isin(['WATHTE'])
    bool_groepering = locations['Groepering.Code'].isin([''])
    selected = locations.loc[bool_procestype & bool_grootheid & bool_hoedanigheid & bool_groepering & bool_stations]
    selected_withtimemax = dfmt.observations.rwsddl_ssh_get_time_max(selected)
    assert len(selected) == 3
    assert len(selected_withtimemax) == 3
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
def test_ssc_add_linked_stations():
    ssc_catalog_gpd = ssc_ssh_read_catalog(linked_stations=True)
    abas_row = ssc_catalog_gpd.loc["SSC-abas"]
    dist_dict = abas_row["dist_dict"][0]
    assert set(dist_dict.keys()) == set({'IOC: abas', 'UHSLC: 347'})
    assert np.isclose(abas_row["dist_min"], 0.00047381430963381466)
    assert np.isclose(abas_row["dist_max"], 0.0074595241134952145)


@pytest.mark.unittest
def test_ssh_catalog_toxynfile(tmp_path):
    ssc_catalog_gpd = dfmt.ssh_catalog_subset(source="ssc")
    file_xyn = tmp_path / 'test_ssc_obs.xyn'
    dfmt.ssh_catalog_toxynfile(ssc_catalog_gpd, file_xyn)
    assert os.path.isfile(file_xyn)


@pytest.mark.unittest
def test_ssh_catalog_tokmlfile(tmp_path):
    ssc_catalog_gpd = dfmt.ssh_catalog_subset(source="ssc")
    file_kml = tmp_path / 'test_ssc.kml'
    dfmt.ssh_catalog_tokmlfile(ssc_catalog_gpd, file_kml)
    assert os.path.isfile(file_kml)


@pytest.mark.unittest
def test_gtsm3_era5_cds_ssh_retrieve_data_invalidfreq_nonetimes():
    df = gtsm3_era5_cds_ssh_read_catalog()
    
    with pytest.raises(ValueError) as e:
        gtsm3_era5_cds_ssh_retrieve_data(
            row=df.iloc[0],
            time_min=None,
            time_max=None,
            time_freq='10min',
            )
    assert "time frequency for retrieving gtsm3-era5-cds data should" in str(e.value)


@pytest.mark.unittest
def test_remove_accents(caplog):
    # ø is a non-ascii character replaced by o in the function, to avoid dropping
    output_str = _remove_accents("Måløy")
    assert output_str == "Maloy"
    
    # ð is a non-ascii character that is not accounted for in the function, so it is dropped
    with caplog.at_level(logging.WARNING):
        _remove_accents("Måløyð")
    assert "_remove_accents() dropped characters: 'Måløyð' became 'Maloy'" in caplog.text
