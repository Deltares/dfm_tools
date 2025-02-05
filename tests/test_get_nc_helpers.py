# -*- coding: utf-8 -*-
"""
Created on Mon Feb  3 11:42:28 2025

@author: veenstra
"""
import pytest
import os
import dfm_tools as dfmt


@pytest.mark.unittest
def test_get_ncvarproperties():
    uds = dfmt.data.fm_curvedbend_map()
    vars_pd = dfmt.get_ncvarproperties(uds)
    expected_columns = ['shape', 'dimensions', 'dtype', 'units',
           'standard_name', 'long_name', 'mesh', 'location', 'grid_mapping',
           'bounds', 'formula_terms']
    
    assert len(vars_pd) == 43
    assert len(vars_pd.columns) == 26
    for colname in expected_columns:
        assert colname in vars_pd.columns
    assert vars_pd.index.str.startswith("mesh2d_").sum() == 40
    assert vars_pd.index.str.startswith("time").sum() == 2
    assert "projected_coordinate_system" in vars_pd.index


@pytest.mark.unittest
def test_rename_waqvars():
    # originally tested with r'p:\archivedprojects\11203850-coastserv\06-Model\waq_model\simulations\run0_20200319\DFM_OUTPUT_kzn_waq\kzn_waq_0000_map.nc'
    # but converted to a portable testcase by setting up a dummy dataset
    file_nc = dfmt.data.fm_curvedbend_map(return_filepath=True)
    uds = dfmt.open_partitioned_dataset(file_nc)
    
    # add dummy waq variable
    waq_attrs = {'mesh': 'mesh2d',
     'location': 'face',
     'cell_methods': 'mesh2d_nFaces: mean',
     'long_name': 'Chlfa',
     'units': '(mg/m3)',
     'grid_mapping': 'wgs84',
     'description': 'Chlfa - Chlorophyll-a concentration (mg/m3) in flow element'}
    uds['mesh2d_water_quality_output_61'] = uds['mesh2d_flowelem_bl'].assign_attrs(waq_attrs)
    uds['mesh2d_water_quality_output_62'] = uds['mesh2d_flowelem_bl'].assign_attrs(waq_attrs)
    
    uds = dfmt.rename_waqvars(uds)
    assert 'mesh2d_Chlfa' in uds.data_vars
    # in case of duplicates, only the first one is renamed
    assert 'mesh2d_water_quality_output_62' in uds.data_vars


@pytest.mark.unittest
def test_rename_fouvars_regular():
    file_nc_fou = dfmt.data.fm_westernscheldt_fou(return_filepath=True)
    uds = dfmt.open_partitioned_dataset(file_nc_fou)
    uds_renamed = dfmt.rename_fouvars(uds)
    
    assert 'mesh2d_ux_mean_20140101000000_20140101000600' in uds_renamed.data_vars
    assert 'mesh2d_uy_mean_20140101000000_20140101000600' in uds_renamed.data_vars
    
    
@pytest.mark.requireslocaldata
@pytest.mark.unittest
def test_rename_fouvars_tidal():
    file_nc_fou = r'p:\1230882-emodnet_hrsm\GTSMv3.0EMODnet\EMOD_MichaelTUM_yearcomponents\GTSMv4.1_yeartide_2014_2.20.06\output\gtsm_model_0000_fou.nc'
    uds = dfmt.open_partitioned_dataset(file_nc_fou)
    uds_renamed = dfmt.rename_fouvars(uds,drop_tidal_times=False)
    uds_renamed_clean = dfmt.rename_fouvars(uds)
    
    assert 'mesh2d_s1_mean_20131231000000_20150101000000' in uds_renamed.data_vars
    assert 'mesh2d_s1_min_20131231000000_20150101000000' in uds_renamed.data_vars
    assert 'mesh2d_s1_mindepth_20131231000000_20150101000000' in uds_renamed.data_vars
    assert 'mesh2d_s1_ampM2_20131231000000_20141227180000' in uds_renamed.data_vars
    assert 'mesh2d_s1_mean_20131231000000_20150101000000' in uds_renamed_clean.data_vars
    assert 'mesh2d_s1_min_20131231000000_20150101000000' in uds_renamed_clean.data_vars
    assert 'mesh2d_s1_mindepth_20131231000000_20150101000000' in uds_renamed_clean.data_vars
    assert 'mesh2d_s1_ampM2' in uds_renamed_clean.data_vars
