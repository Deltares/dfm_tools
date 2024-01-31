#!/usr/bin/env python

import pytest
import os
import dfm_tools as dfmt
import numpy as np
import hydrolib.core.dflowfm as hcdfm
import pandas as pd
from dfm_tools.get_nc import (calc_dist_pythagoras,
                              calc_dist_haversine,
                              intersect_edges_withsort,
                              reconstruct_zw_zcc_fromz,
                              reconstruct_zw_zcc_fromzsigma,
                              reconstruct_zw_zcc_fromsigma
                              )


@pytest.mark.parametrize("file_nc, expected_size", [pytest.param(dfmt.data.fm_grevelingen_map(return_filepath=True), (44796, 4, 2), id='from partitioned map Grevelingen'),
                                                    pytest.param(dfmt.data.fm_curvedbend_map(return_filepath=True), (550, 4, 2), id='from map curvedbend'),
                                                    pytest.param(dfmt.data.fm_grevelingen_net(return_filepath=True), (44804,4,2), id='fromnet Grevelingen')])
@pytest.mark.systemtest
def test_facenodecoordinates_shape(file_nc, expected_size):
    
    uds = dfmt.open_partitioned_dataset(file_nc)
    facenodecoordinates = uds.grid.face_node_coordinates
    
    assert facenodecoordinates.shape == expected_size


@pytest.mark.unittest
def test_getmapdata():
    """
    Checks whether ghost cells are properly taken care of (by asserting shape). And also whether varname can be found from attributes in case of Chlfa.
    
    """
    file_nc = dfmt.data.fm_grevelingen_map(return_filepath=True)
    expected_size = (44796, 36)
    varname = 'mesh2d_sa1'
    
    data_xr_map = dfmt.open_partitioned_dataset(file_nc)
    data_xr_map = dfmt.rename_waqvars(data_xr_map)
    data_varsel = data_xr_map[varname].isel(time=2)
    
    assert data_varsel.shape == expected_size


@pytest.mark.requireslocaldata
@pytest.mark.unittest
def test_getmapdata_waq():
    """
    Checks whether ghost cells are properly taken care of (by asserting shape). And also whether varname can be found from attributes in case of Chlfa.
    
    """
    file_nc = os.path.join(r'p:\archivedprojects\11203850-coastserv\06-Model\waq_model\simulations\run0_20200319\DFM_OUTPUT_kzn_waq', 'kzn_waq_0*_map.nc')
    expected_size = (17385, 39)
    varname = 'mesh2d_Chlfa'
    
    data_xr_map = dfmt.open_partitioned_dataset(file_nc)
    data_xr_map = dfmt.rename_waqvars(data_xr_map)
    data_varsel = data_xr_map[varname].isel(time=2)
    
    assert data_varsel.shape == expected_size


@pytest.mark.requireslocaldata
@pytest.mark.unittest
def test_rename_waqvars():
    file_nc = os.path.join(r'p:\archivedprojects\11203850-coastserv\06-Model\waq_model\simulations\run0_20200319\DFM_OUTPUT_kzn_waq', 'kzn_waq_0000_map.nc')
    uds = dfmt.open_partitioned_dataset(file_nc)
    uds = dfmt.rename_waqvars(uds)
    
    assert 'mesh2d_Chlfa' in uds.data_vars


@pytest.mark.requireslocaldata
@pytest.mark.unittest
def test_rename_fouvars_regular():
    file_nc_fou = r'p:\archivedprojects\11203379-005-mwra-updated-bem\03_model\02_final\A72_ntsu0_kzlb2\DFM_OUTPUT_MB_02_fou\MB_02_0000_fou.nc'
    uds = dfmt.open_partitioned_dataset(file_nc_fou)
    uds_renamed = dfmt.rename_fouvars(uds)
    
    assert 'mesh2d_tem_mean_20160201000000_20160301000000' in uds_renamed.data_vars
    assert 'mesh2d_uy_mean_20160101000000_20170101000000' in uds_renamed.data_vars
    
    
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


@pytest.mark.unittest
def test_calc_dist_pythagoras():
    """
    all crossings for cb3 and testline (which has linebend in cell en with line crosses same cell twice)
    """
    edge_index = np.array([0,0,0,1,1,1,2,2])

    edges = np.array([[[2084.67741935, 3353.02419355],
                       [2255.79637097, 3307.15725806]],
                      [[2255.79637097, 3307.15725806],
                       [2222.27822581, 3206.60282258]],
                      [[2222.27822581, 3206.60282258],
                       [2128.78024194, 3266.58266129]]])

    intersections = np.array([[[2084.67741935, 3353.02419355],
                               [2144.15041424, 3337.08297842]],
                              [[2144.15041424, 3337.08297842],
                               [2202.53662217, 3321.43306702]],
                              [[2202.53662217, 3321.43306702],
                               [2255.79637097, 3307.15725806]],
                              [[2255.79637097, 3307.15725806],
                               [2246.9810802,  3280.71138574]],
                              [[2246.9810802,  3280.71138574],
                               [2239.02015401, 3256.82860719]],
                              [[2239.02015401, 3256.82860719],
                               [2222.27822581, 3206.60282258]],
                              [[2222.27822581, 3206.60282258],
                               [2173.05750857, 3238.17837704]],
                              [[2173.05750857, 3238.17837704],
                               [2128.78024194, 3266.58266129]]])
    
    edge_len = calc_dist_pythagoras(edges[:,0,0], edges[:,1,0], edges[:,0,1], edges[:,1,1])
    edge_len_cum = np.cumsum(edge_len)
    edge_len_cum0 = np.concatenate([[0],edge_len_cum[:-1]])
    crs_dist_starts = calc_dist_pythagoras(edges[edge_index,0,0], intersections[:,0,0], edges[edge_index,0,1], intersections[:,0,1]) + edge_len_cum0[edge_index]
    crs_dist_stops  = calc_dist_pythagoras(edges[edge_index,0,0], intersections[:,1,0], edges[edge_index,0,1], intersections[:,1,1]) + edge_len_cum0[edge_index]
    
    crs_dist_starts_check = np.array([  0.        ,  61.57239204, 122.01963352, 177.15945184,
                                      205.03584892, 230.21050794, 283.15313349, 341.63128877])
    crs_dist_stops_check = np.array([ 61.57239204, 122.01963352, 177.15945184, 205.03584892,
                                     230.21050794, 283.15313349, 341.63128877, 394.23622869])
    
    assert np.allclose(crs_dist_starts, crs_dist_starts_check)
    assert np.allclose(crs_dist_stops, crs_dist_stops_check)


@pytest.mark.unittest
def test_calc_dist_haversine():
    """
    first 15 crossings for DSCM
    """
    edge_index = np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
 
    edges = np.array([[[ 8.92659074, 56.91538014],
                       [ 8.58447136, 58.66874192]]])
  
    intersections = np.array([[[ 8.8856893 , 57.125     ],
                               [ 8.88406329, 57.13333333]],
                              [[ 8.88406329, 57.13333333],
                               [ 8.88243727, 57.14166667]],
                              [[ 8.88243727, 57.14166667],
                               [ 8.88081125, 57.15      ]],
                              [[ 8.88081125, 57.15      ],
                               [ 8.87918524, 57.15833333]],
                              [[ 8.87918524, 57.15833333],
                               [ 8.87755922, 57.16666667]],
                              [[ 8.87755922, 57.16666667],
                               [ 8.87593321, 57.175     ]],
                              [[ 8.87593321, 57.175     ],
                               [ 8.875     , 57.17978268]],
                              [[ 8.875     , 57.17978268],
                               [ 8.87430719, 57.18333333]],
                              [[ 8.87430719, 57.18333333],
                               [ 8.87268117, 57.19166667]],
                              [[ 8.87268117, 57.19166667],
                               [ 8.87105516, 57.2       ]],
                              [[ 8.87105516, 57.2       ],
                               [ 8.86942914, 57.20833333]],
                              [[ 8.86942914, 57.20833333],
                               [ 8.86780312, 57.21666667]],
                              [[ 8.86780312, 57.21666667],
                               [ 8.86617711, 57.225     ]],
                              [[ 8.86617711, 57.225     ],
                               [ 8.86455109, 57.23333333]],
                              [[ 8.86455109, 57.23333333],
                               [ 8.86292507, 57.24166667]]])
    
    edge_len = calc_dist_haversine(edges[:,0,0], edges[:,1,0], edges[:,0,1], edges[:,1,1])
    edge_len_cum = np.cumsum(edge_len)
    edge_len_cum0 = np.concatenate([[0],edge_len_cum[:-1]])
    crs_dist_starts = calc_dist_haversine(edges[edge_index,0,0], intersections[:,0,0], edges[edge_index,0,1], intersections[:,0,1]) + edge_len_cum0[edge_index]
    crs_dist_stops  = calc_dist_haversine(edges[edge_index,0,0], intersections[:,1,0], edges[edge_index,0,1], intersections[:,1,1]) + edge_len_cum0[edge_index]
    
    crs_dist_starts_check = np.array([23439.77082715, 24371.57628696, 25303.38057682, 26235.18142118,
                                      27166.97986164, 28098.77713142, 29030.57089118, 29565.34666042,
                                      29962.36237409, 30894.15262176, 31825.93935877, 32757.7238182 ,
                                      33689.50704171, 34621.28675393, 35553.06418784])
    crs_dist_stops_check = np.array([24371.57628696, 25303.38057682, 26235.18142118, 27166.97986164,
                                     28098.77713142, 29030.57089118, 29565.34666042, 29962.36237409,
                                     30894.15262176, 31825.93935877, 32757.7238182 , 33689.50704171,
                                     34621.28675393, 35553.06418784, 36484.84038514])
    
    assert np.allclose(crs_dist_starts, crs_dist_starts_check)
    assert np.allclose(crs_dist_stops, crs_dist_stops_check)


@pytest.mark.unittest
def test_intersect_edges():
    """
    ordering of xu.ugrid2d.intersect_edges return arrays is wrong, but we test it here since this test will fail once sorting is fixed in xugrid or numba.celltree. If so, depracate dfmt.intersect_edges_withsort()
    """
    
    file_nc = dfmt.data.fm_curvedbend_map(return_filepath=True) #sigmalayer
                    
    line_array = np.array([[2084.67741935, 3353.02419355], #with linebend in cell en with line crossing same cell twice
                           [2255.79637097, 3307.15725806],
                           [2222.27822581, 3206.60282258],
                           [2128.78024194, 3266.58266129]])
    
    uds = dfmt.open_partitioned_dataset(file_nc)
    
    edges = np.stack([line_array[:-1],line_array[1:]],axis=1)
    edge_index, face_index, intersections = uds.grid.intersect_edges(edges)
    
    assert (edge_index == np.array([0, 0, 0, 1, 1, 1, 2, 2])).all()
    assert (face_index == np.array([ 91, 146, 147, 202, 147, 201, 146, 201])).all()


@pytest.mark.unittest
def test_intersect_edges_withsort():
    """
    ordering of xu.ugrid2d.intersect_edges return arrays is wrong, so dfmt.intersect_edges_withsort() combines it with sorting. The line array clearly shows different ordering of the resulting face_index array
    """
    
    file_nc = dfmt.data.fm_curvedbend_map(return_filepath=True) #sigmalayer
    
    line_array = np.array([[2084.67741935, 3353.02419355], #with linebend in cell en with line crossing same cell twice
                           [2255.79637097, 3307.15725806],
                           [2222.27822581, 3206.60282258],
                           [2128.78024194, 3266.58266129]])
    
    uds = dfmt.open_partitioned_dataset(file_nc)
    
    edges = np.stack([line_array[:-1],line_array[1:]],axis=1)
    edge_index, face_index, intersections = intersect_edges_withsort(uds,edges)
    
    assert (edge_index == np.array([0, 0, 0, 1, 1, 1, 2, 2])).all()
    assert (face_index == np.array([ 91, 146, 147, 147, 202, 201, 201, 146])).all()


@pytest.mark.unittest
def test_zlayermodel_correct_layers():
    """
    we assert max/min only, since zlayers can also be valid when the top/bottom layer contains nan values
    """
    
    file_nc = dfmt.data.fm_grevelingen_map(return_filepath=True) #zlayer
    data_frommap_merged = dfmt.open_partitioned_dataset(file_nc)
    
    timestep = 3
    
    data_frommap_timesel = data_frommap_merged.isel(time=timestep) #select data for all layers
    data_frommap_merged_fullgrid = reconstruct_zw_zcc_fromz(data_frommap_timesel)
    
    vals_wl = data_frommap_merged_fullgrid['mesh2d_s1']
    vals_bl = data_frommap_merged_fullgrid['mesh2d_flowelem_bl'].to_numpy()
    
    vals_zw_max = data_frommap_merged_fullgrid['mesh2d_flowelem_zw'].max(dim='nmesh2d_interface')
    vals_zw_min = data_frommap_merged_fullgrid['mesh2d_flowelem_zw'].min(dim='nmesh2d_interface')
    assert np.allclose(vals_zw_max, vals_wl)
    assert np.allclose(vals_zw_min, vals_bl)
    
    # check z-centers in one specific cell
    zcc = data_frommap_merged_fullgrid['mesh2d_flowelem_zcc']
    zcc_onecell = zcc.isel(nmesh2d_face=5000).load().fillna(-999)
    zcc_onecell_expected = np.array([-999]*32 + [ -4.1403    , -3.375     , -2.125, -0.84786265])
    assert np.allclose(zcc_onecell, zcc_onecell_expected)
    
    # check z-interfaces in one specific cell
    zw = data_frommap_merged_fullgrid['mesh2d_flowelem_zw']
    zw_onecell = zw.isel(nmesh2d_face=5000).load().fillna(-999)
    zw_onecell_expected = np.array([-999]*32 + [-4.2806, -4.0, -2.75, -1.5, -0.195725292])
    assert np.allclose(zw_onecell, zw_onecell_expected)
    
    # check if all non-dry centers are below waterlevel and above bed
    vals_zcc_max = data_frommap_merged_fullgrid['mesh2d_flowelem_zcc'].max(dim='nmesh2d_layer').to_numpy()
    vals_zcc_min = data_frommap_merged_fullgrid['mesh2d_flowelem_zcc'].min(dim='nmesh2d_layer').to_numpy()
    bool_dry = vals_wl == vals_bl
    assert ((vals_zcc_max < vals_wl) | bool_dry).all()
    assert ((vals_zcc_min > vals_bl) | bool_dry).all()


@pytest.mark.unittest
def test_zlayermodel_correct_layers_nanabovewl():
    """
    we assert max/min only, since zlayers can also be valid when the top/bottom layer contains nan values
    """
    
    file_nc = dfmt.data.fm_grevelingen_map(return_filepath=True) #zlayer
    uds = dfmt.open_partitioned_dataset(file_nc)
    
    # set wl/bl below the second interface (first is at -0.25, second is at -1.5)
    uds['mesh2d_s1'] = uds['mesh2d_s1'].clip(max=-1.8)  # assuming this variable is available, which is not guaranteed
    uds['mesh2d_flowelem_bl'] = uds['mesh2d_flowelem_bl'].clip(max=-1.8)  # assuming this variable is available, which is not guaranteed

    timestep = 3
    uds_timesel = uds.isel(time=timestep) #select data for all layers
    uds_timesel_fullgrid = reconstruct_zw_zcc_fromz(uds_timesel)
    
    # check z-centers in one specific cell
    zcc = uds_timesel_fullgrid['mesh2d_flowelem_zcc']
    zcc_onecell = zcc.isel(nmesh2d_face=5000).load().fillna(-999)
    zcc_onecell_expected = np.array([-999]*32 + [ -4.1403    , -3.375     , -2.275, -999.])
    assert np.allclose(zcc_onecell, zcc_onecell_expected)
    
    # check z-interfaces in one specific cell
    zw = uds_timesel_fullgrid['mesh2d_flowelem_zw']
    zw_onecell = zw.isel(nmesh2d_face=5000).load().fillna(-999)
    zw_onecell_expected = np.array([-999]*32 + [-4.2806, -4.0, -2.75, -1.8, -999.])
    assert np.allclose(zw_onecell, zw_onecell_expected)


@pytest.mark.requireslocaldata
@pytest.mark.unittest
def test_zsigmalayermodel_correct_layers():
    """
    we assert top/max/min only, since zlayers can also be valid when the bottom layer contains nan values
    """
    
    file_nc = r'p:\dflowfm\maintenance\JIRA\05000-05999\05477\c103_ws_3d_fourier\DFM_OUTPUT_westerscheldt01_0subst\westerscheldt01_0subst_map.nc' #zsigma model without fullgrid output but with new ocean_sigma_z_coordinate variable
    data_frommap_merged = dfmt.open_partitioned_dataset(file_nc)
    
    timestep = 1
    
    data_frommap_timesel = data_frommap_merged.isel(time=timestep) #select data for all layers
    data_frommap_merged_fullgrid = reconstruct_zw_zcc_fromzsigma(data_frommap_timesel)
    
    vals_wl = data_frommap_merged_fullgrid['mesh2d_s1'].to_numpy()
    vals_bl = data_frommap_merged_fullgrid['mesh2d_flowelem_bl'].to_numpy()

    # check z-centers in one specific cell
    zcc = data_frommap_merged_fullgrid['mesh2d_flowelem_zcc']
    zcc_onecell = zcc.isel(mesh2d_nFaces=5000).load().fillna(-999)
    zcc_onecell_expected = np.array([-999, -999, -999, -4.86792313, -0.250909869])
    assert np.allclose(zcc_onecell, zcc_onecell_expected)
    
    # check z-interfaces in one specific cell
    zw = data_frommap_merged_fullgrid['mesh2d_flowelem_zw']
    zw_onecell = zw.isel(mesh2d_nFaces=5000).load().fillna(-999)
    zw_onecell_expected = np.array([-999., -999., -999.,   -7.17642976, -2.5594165, 2.05759676])
    assert np.allclose(zw_onecell, zw_onecell_expected)
    
    vals_zw_max = data_frommap_merged_fullgrid['mesh2d_flowelem_zw'].max(dim='mesh2d_nInterfaces').to_numpy()
    vals_zw_min = data_frommap_merged_fullgrid['mesh2d_flowelem_zw'].min(dim='mesh2d_nInterfaces').to_numpy()
    vals_zw_top = data_frommap_merged_fullgrid['mesh2d_flowelem_zw'].isel(mesh2d_nInterfaces=-1).to_numpy()
    assert np.allclose(vals_zw_max, vals_wl)
    assert np.allclose(vals_zw_min, vals_bl)
    assert np.allclose(vals_zw_max, vals_zw_top)

    # check if all non-dry cell centers are below waterlevel and above bed
    vals_zcc_max = data_frommap_merged_fullgrid['mesh2d_flowelem_zcc'].max(dim='mesh2d_nLayers').to_numpy()
    vals_zcc_min = data_frommap_merged_fullgrid['mesh2d_flowelem_zcc'].min(dim='mesh2d_nLayers').to_numpy()
    vals_zcc_top = data_frommap_merged_fullgrid['mesh2d_flowelem_zcc'].isel(mesh2d_nLayers=-1).to_numpy()
    bool_dry = vals_wl == vals_bl
    assert ((vals_zcc_max < vals_wl) | bool_dry).all()
    assert ((vals_zcc_min > vals_bl) | bool_dry).all()
    assert (np.isclose(vals_zcc_max, vals_zcc_top) | bool_dry).all()


@pytest.mark.unittest
def test_sigmalayermodel_correct_layers():
    """
    we assert top/bot/max/min, since max/top and min/bottom sigmalayers are always aligned and do not contain nans
    """
    file_nc = dfmt.data.fm_curvedbend_map(return_filepath=True) #sigmalayer
    data_frommap_merged = dfmt.open_partitioned_dataset(file_nc)
    
    timestep = 3
    
    data_frommap_timesel = data_frommap_merged.isel(time=timestep) #select data for all layers
    data_frommap_merged_fullgrid = reconstruct_zw_zcc_fromsigma(data_frommap_timesel)
    
    vals_wl = data_frommap_merged_fullgrid['mesh2d_s1'].to_numpy()
    vals_bl = data_frommap_merged_fullgrid['mesh2d_flowelem_bl'].to_numpy()

    vals_zw_max = data_frommap_merged_fullgrid['mesh2d_flowelem_zw'].max(dim='mesh2d_nInterfaces').to_numpy()
    vals_zw_min = data_frommap_merged_fullgrid['mesh2d_flowelem_zw'].min(dim='mesh2d_nInterfaces').to_numpy()
    vals_zw_top = data_frommap_merged_fullgrid['mesh2d_flowelem_zw'].isel(mesh2d_nInterfaces=-1).to_numpy()
    vals_zw_bot = data_frommap_merged_fullgrid['mesh2d_flowelem_zw'].isel(mesh2d_nInterfaces=0).to_numpy()
    assert (np.abs(vals_zw_max-vals_wl)<1e-6).all()
    assert (np.abs(vals_zw_min-vals_bl)<1e-6).all()
    assert (np.abs(vals_zw_max-vals_zw_top)<1e-6).all()
    assert (np.abs(vals_zw_min-vals_zw_bot)<1e-6).all()

    vals_zcc_max = data_frommap_merged_fullgrid['mesh2d_flowelem_zcc'].max(dim='mesh2d_nLayers').to_numpy()
    vals_zcc_min = data_frommap_merged_fullgrid['mesh2d_flowelem_zcc'].min(dim='mesh2d_nLayers').to_numpy()
    vals_zcc_top = data_frommap_merged_fullgrid['mesh2d_flowelem_zcc'].isel(mesh2d_nLayers=-1).to_numpy()
    vals_zcc_bot = data_frommap_merged_fullgrid['mesh2d_flowelem_zcc'].isel(mesh2d_nLayers=0).to_numpy()
    assert (vals_zcc_max < vals_wl).all() # < works since there are no bl>0 in this model
    assert (vals_zcc_min > vals_bl).all()
    assert (np.abs(vals_zcc_max-vals_zcc_top)<1e-6).all()
    assert (np.abs(vals_zcc_min-vals_zcc_bot)<1e-6).all()


@pytest.mark.requireslocaldata
def test_timmodel_to_dataframe():
    
    file_tim = r'p:\archivedprojects\11206811-002-d-hydro-grevelingen\simulaties\Jaarsom2017_dfm_006_zlayer\boundary_conditions\hist\jaarsom_2017\sources_sinks\FlakkeeseSpuisluis.tim'
    
    data_tim = hcdfm.TimModel(file_tim)
    
    refdate = '2016-01-01'
    tim_pd = dfmt.TimModel_to_DataFrame(data_tim, parse_column_labels=True, refdate=refdate)
    
    assert tim_pd.index[0] == pd.Timestamp('2016-01-01 00:00:00')
    assert len(tim_pd) == 39603
    assert tim_pd.columns[-1] == 'Phaeocystis_P (g/m3)'

