# -*- coding: utf-8 -*-
"""
Created on Wed Dec 13 12:17:41 2023

@author: veenstra
"""

import pytest
import numpy as np
import dfm_tools as dfmt
import xarray as xr
from dfm_tools.get_nc import (calc_dist_pythagoras,
                              calc_dist_haversine,
                              intersect_edges_withsort,
                              reconstruct_zw_zcc_fromz,
                              reconstruct_zw_zcc_fromzsigma,
                              reconstruct_zw_zcc_fromsigma
                              )


@pytest.mark.unittest
def test_zlayermodel_correct_layers():
    """
    we assert max/min only, since zlayers can also be valid when the top/bottom layer contains nan values
    """
    
    file_nc = dfmt.data.fm_grevelingen_map(return_filepath=True) #zlayer
    uds = dfmt.open_partitioned_dataset(file_nc)
    
    timestep = 3
    
    uds_timesel = uds.isel(time=timestep) #select data for all layers
    uds_fullgrid = reconstruct_zw_zcc_fromz(uds_timesel)
    
    vals_wl = uds_fullgrid['mesh2d_s1']
    vals_bl = uds_fullgrid['mesh2d_flowelem_bl'].to_numpy()
    
    vals_zw_max = uds_fullgrid['mesh2d_flowelem_zw'].max(dim='nmesh2d_interface')
    vals_zw_min = uds_fullgrid['mesh2d_flowelem_zw'].min(dim='nmesh2d_interface')
    assert np.allclose(vals_zw_max, vals_wl)
    assert np.allclose(vals_zw_min, vals_bl)
    
    # check z-centers in one specific cell
    zcc = uds_fullgrid['mesh2d_flowelem_zcc']
    zcc_onecell = zcc.isel(nmesh2d_face=5000).load().fillna(-999)
    zcc_onecell_expected = np.array([-999]*32 + [ -4.1403    , -3.375     , -2.125, -0.84786265])
    assert np.allclose(zcc_onecell, zcc_onecell_expected)
    
    # check z-interfaces in one specific cell
    zw = uds_fullgrid['mesh2d_flowelem_zw']
    zw_onecell = zw.isel(nmesh2d_face=5000).load().fillna(-999)
    zw_onecell_expected = np.array([-999]*32 + [-4.2806, -4.0, -2.75, -1.5, -0.195725292])
    assert np.allclose(zw_onecell, zw_onecell_expected)
    
    # check if all non-dry centers are below waterlevel and above bed
    vals_zcc_max = uds_fullgrid['mesh2d_flowelem_zcc'].max(dim='nmesh2d_layer').to_numpy()
    vals_zcc_min = uds_fullgrid['mesh2d_flowelem_zcc'].min(dim='nmesh2d_layer').to_numpy()
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


@pytest.mark.unittest
def test_zsigmalayermodel_correct_layers():
    """
    we assert top/max/min only, since zlayers can also be valid when the bottom layer contains nan values
    """
    
    # zsigma model without fullgrid output but with new ocean_sigma_z_coordinate variable
    file_nc = dfmt.data.fm_westernscheldt_map(return_filepath=True)
    uds = dfmt.open_partitioned_dataset(file_nc)
    
    timestep = 1
    
    uds_timesel = uds.isel(time=timestep) #select data for all layers
    uds_fullgrid = reconstruct_zw_zcc_fromzsigma(uds_timesel)
    
    vals_wl = uds_fullgrid['mesh2d_s1'].to_numpy()
    vals_bl = uds_fullgrid['mesh2d_flowelem_bl'].to_numpy()

    # check z-centers in one specific cell
    zcc = uds_fullgrid['mesh2d_flowelem_zcc']
    zcc_onecell = zcc.isel(mesh2d_nFaces=5000).load().fillna(-999)
    zcc_onecell_expected = np.array([-999, -999, -999, -4.86792313, -0.250909869])
    assert np.allclose(zcc_onecell, zcc_onecell_expected)
    
    # check z-interfaces in one specific cell
    zw = uds_fullgrid['mesh2d_flowelem_zw']
    zw_onecell = zw.isel(mesh2d_nFaces=5000).load().fillna(-999)
    zw_onecell_expected = np.array([-999., -999., -999.,   -7.17642976, -2.5594165, 2.05759676])
    assert np.allclose(zw_onecell, zw_onecell_expected)
    
    vals_zw_max = uds_fullgrid['mesh2d_flowelem_zw'].max(dim='mesh2d_nInterfaces').to_numpy()
    vals_zw_min = uds_fullgrid['mesh2d_flowelem_zw'].min(dim='mesh2d_nInterfaces').to_numpy()
    vals_zw_top = uds_fullgrid['mesh2d_flowelem_zw'].isel(mesh2d_nInterfaces=-1).to_numpy()
    assert np.allclose(vals_zw_max, vals_wl)
    assert np.allclose(vals_zw_min, vals_bl)
    assert np.allclose(vals_zw_max, vals_zw_top)

    # check if all non-dry cell centers are below waterlevel and above bed
    vals_zcc_max = uds_fullgrid['mesh2d_flowelem_zcc'].max(dim='mesh2d_nLayers').to_numpy()
    vals_zcc_min = uds_fullgrid['mesh2d_flowelem_zcc'].min(dim='mesh2d_nLayers').to_numpy()
    vals_zcc_top = uds_fullgrid['mesh2d_flowelem_zcc'].isel(mesh2d_nLayers=-1).to_numpy()
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
    uds_fullgrid = reconstruct_zw_zcc_fromsigma(data_frommap_timesel)
    
    vals_wl = uds_fullgrid['mesh2d_s1'].to_numpy()
    vals_bl = uds_fullgrid['mesh2d_flowelem_bl'].to_numpy()

    vals_zw_max = uds_fullgrid['mesh2d_flowelem_zw'].max(dim='mesh2d_nInterfaces').to_numpy()
    vals_zw_min = uds_fullgrid['mesh2d_flowelem_zw'].min(dim='mesh2d_nInterfaces').to_numpy()
    vals_zw_top = uds_fullgrid['mesh2d_flowelem_zw'].isel(mesh2d_nInterfaces=-1).to_numpy()
    vals_zw_bot = uds_fullgrid['mesh2d_flowelem_zw'].isel(mesh2d_nInterfaces=0).to_numpy()
    assert (np.abs(vals_zw_max-vals_wl)<1e-6).all()
    assert (np.abs(vals_zw_min-vals_bl)<1e-6).all()
    assert (np.abs(vals_zw_max-vals_zw_top)<1e-6).all()
    assert (np.abs(vals_zw_min-vals_zw_bot)<1e-6).all()

    vals_zcc_max = uds_fullgrid['mesh2d_flowelem_zcc'].max(dim='mesh2d_nLayers').to_numpy()
    vals_zcc_min = uds_fullgrid['mesh2d_flowelem_zcc'].min(dim='mesh2d_nLayers').to_numpy()
    vals_zcc_top = uds_fullgrid['mesh2d_flowelem_zcc'].isel(mesh2d_nLayers=-1).to_numpy()
    vals_zcc_bot = uds_fullgrid['mesh2d_flowelem_zcc'].isel(mesh2d_nLayers=0).to_numpy()
    assert (vals_zcc_max < vals_wl).all() # < works since there are no bl>0 in this model
    assert (vals_zcc_min > vals_bl).all()
    assert (np.abs(vals_zcc_max-vals_zcc_top)<1e-6).all()
    assert (np.abs(vals_zcc_min-vals_zcc_bot)<1e-6).all()


@pytest.mark.unittest
def test_polyline_mapslice():
    uds = dfmt.data.fm_curvedbend_map()
    timestep = 72
    line_array = np.array([[ 104.15421399, 2042.7077107 ],
                            [2913.47878063, 2102.48057382]])
    
    uds_crs = dfmt.polyline_mapslice(uds.isel(time=timestep), line_array)
    assert len(uds_crs.grid.node_x) == 720
    assert uds_crs.grid.face_node_connectivity.shape == (180, 4)
    assert np.isclose(uds_crs.grid.node_x.min(), 484.02691589067194)
    assert np.isclose(uds_crs.grid.node_x.max(), 1737.0864257281437)
    assert np.isclose(uds_crs.grid.node_y.min(), -5)
    assert np.isclose(uds_crs.grid.node_y.max(), 0.9261683648147339)


@pytest.mark.unittest
def test_get_dataset_atdepths_mapfile():
    file_nc = dfmt.data.fm_curvedbend_map(return_filepath=True)
    uds = dfmt.open_partitioned_dataset(file_nc)

    # z0
    depths = [-1,-4,0,-6]
    uds_atdepths = dfmt.get_Dataset_atdepths(data_xr=uds, depths=depths, reference='z0')
    tem_values = uds_atdepths.mesh2d_tem1.isel(time=-1, mesh2d_nFaces=10).to_numpy()
    exp_values = np.array([np.nan, 15., 15., 15.])
    assert uds_atdepths.mesh2d_tem1.shape == (73, 550, 4)
    assert np.allclose(tem_values, exp_values, equal_nan=True)
    assert np.isclose(uds_atdepths.mesh2d_tem1.sum(), 1511040)
    
    # waterlevel
    depths = [-1,-4,0,-6]
    uds_atdepths = dfmt.get_Dataset_atdepths(data_xr=uds, depths=depths, reference='waterlevel')
    tem_values = uds_atdepths.mesh2d_tem1.isel(time=-1, mesh2d_nFaces=10).to_numpy()
    exp_values = np.array([np.nan, 15., 15., np.nan])
    assert uds_atdepths.mesh2d_tem1.shape == (73, 550, 4)
    assert np.allclose(tem_values, exp_values, equal_nan=True)
    assert np.isclose(uds_atdepths.mesh2d_tem1.sum(), 1204500)
    
    # bedlevel
    depths = [-1,-4,0,-6]
    uds_atdepths = dfmt.get_Dataset_atdepths(data_xr=uds, depths=depths, reference='bedlevel')
    tem_values = uds_atdepths.mesh2d_tem1.isel(time=-1, mesh2d_nFaces=10).to_numpy()
    exp_values = np.array([np.nan, np.nan, np.nan, 15.])
    assert uds_atdepths.mesh2d_tem1.shape == (73, 550, 4)
    assert np.allclose(tem_values, exp_values, equal_nan=True)
    assert np.isclose(uds_atdepths.mesh2d_tem1.sum(), 602250)
    
    # single depth
    depths = -4
    uds_atdepths = dfmt.get_Dataset_atdepths(data_xr=uds, depths=depths, reference='z0')
    tem_values = uds_atdepths.mesh2d_tem1.isel(time=-1, mesh2d_nFaces=10).to_numpy()
    exp_values = 15.
    assert uds_atdepths.mesh2d_tem1.shape == (73, 550)
    assert np.isclose(tem_values, exp_values)
    assert np.isclose(uds_atdepths.mesh2d_tem1.sum(), 602250)
    
    # nonexistent reference
    with pytest.raises(KeyError) as e:
        dfmt.get_Dataset_atdepths(data_xr=uds, depths=depths, reference='nonexistent')
    assert 'unknown reference "nonexistent"' in str(e.value)


@pytest.mark.unittest
def test_get_dataset_atdepths_hisfile():
    
    file_nc = dfmt.data.fm_grevelingen_his(return_filepath=True)
    ds = xr.open_dataset(file_nc)#, preprocess=dfmt.preprocess_hisnc)

    depths = [-1,-4,0,-6]
    uds_atdepths = dfmt.get_Dataset_atdepths(data_xr=ds, depths=depths, reference='z0')
    uds_selzt = uds_atdepths.isel(stations=2).isel(time=slice(40,100))

    tem_values = uds_selzt.temperature.isel(time=-1).to_numpy()
    exp_values = np.array([5.81065253, 5.43289777, 5.36916911, 5.36916911])

    assert uds_selzt.temperature.shape == (60,4)
    assert np.allclose(tem_values, exp_values)
    assert np.isclose(uds_selzt.temperature.sum(), 1295.56826688)


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


@pytest.mark.parametrize("sort", [pytest.param(x, id=f"sort={x}") for x in [False, True]])
@pytest.mark.unittest
def test_intersect_edges(sort):
    """
    ordering of xu.ugrid2d.intersect_edges return arrays is wrong.
    
    `sort=False` will fail once sorting is fixed in xugrid or numba.celltree.
    If so, depracate dfmt.intersect_edges_withsort().
    
    `sort=True` tests dfmt.intersect_edges_withsort(), it includes sorting. The
    line array clearly shows different ordering of the resulting face_index.
    
    Once intersect_edges is fixed, move this tests to test_xugrid_helpers.py
    """
    
    uds = dfmt.data.fm_curvedbend_map()
    
    line_array = np.array([[2084.67741935, 3353.02419355], #with linebend in cell en with line crossing same cell twice
                           [2255.79637097, 3307.15725806],
                           [2222.27822581, 3206.60282258],
                           [2128.78024194, 3266.58266129]])
    
    edges = np.stack([line_array[:-1],line_array[1:]],axis=1)
    if sort:
        edge_index, face_index, intersections = intersect_edges_withsort(uds, edges)
        expected_face_index = [ 91, 146, 147, 147, 202, 201, 201, 146]
    else:
        edge_index, face_index, intersections = uds.grid.intersect_edges(edges)
        expected_face_index = [ 91, 146, 147, 202, 147, 201, 146, 201]
    
    assert (edge_index == np.array([0, 0, 0, 1, 1, 1, 2, 2])).all()
    assert (face_index == np.array(expected_face_index)).all()
    assert intersections.shape == (8, 2, 2)


@pytest.mark.unittest
def test_rasterize_ugrid():
    uds = dfmt.data.fm_curvedbend_map()
    
    # rasterize uds with resolution
    ds = dfmt.rasterize_ugrid(uds, resolution=100)
    ds_s1_sel = ds.mesh2d_s1.isel(time=-1, x=slice(None,4), y=slice(None,4))
    expected_values = np.array([[0.99710127, 0.99663291, 0.99585074, 0.99469613],
           [0.99431713, 0.99332472, 0.99181335, 0.98968208],
           [0.99102747, 0.98969936, 0.98777662, 0.98510291],
           [0.98731461, 0.98584209, 0.98372741, 0.98082678]])
    assert np.allclose(ds_s1_sel, expected_values)

    # rasterize uda with resolution
    da = dfmt.rasterize_ugrid(uds.mesh2d_s1, resolution=100)
    da_s1_sel = da.isel(time=-1, x=slice(None,4), y=slice(None,4))
    assert np.allclose(da_s1_sel, expected_values)
    
    # rasterize uds with ds_like
    xmin, ymin, xmax, ymax = 0, 0, 4000, 4000
    d = 120
    regx = np.arange(xmin + 0.5 * d, xmax, d)
    regy = np.arange(ymin + 0.5 * d, ymax, d)
    ds_like = xr.DataArray(np.empty((len(regy), len(regx))), {"y": regy, "x": regx}, ["y", "x"])
    ds = dfmt.rasterize_ugrid(uds, ds_like=ds_like)
    ds_s1_sel = ds.mesh2d_s1.isel(time=-1, x=slice(None,4), y=slice(None,3))
    expected_values = np.array([[0.99710127, 0.99663291, 0.99585074, 0.99469613],
           [0.99431713, 0.99332472, 0.99181335, 0.98968208],
           [0.99102747, 0.98969936, 0.98777662, 0.98510291]])
    assert np.allclose(ds_s1_sel, expected_values)

    # assert the TypeError
    with pytest.raises(TypeError) as e:
        dfmt.rasterize_ugrid(1)
    assert "rasterize_ugrid expected" in str(e.value)


@pytest.mark.unittest
def test_plot_ztdata():
    file_nc = dfmt.data.fm_grevelingen_his(return_filepath=True)
    ds = xr.open_dataset(file_nc)#, preprocess=dfmt.preprocess_hisnc)
    ds_sel = ds.isel(stations=-1)
    
    with pytest.raises(ValueError) as e:
        dfmt.plot_ztdata(data_xr_sel=ds, varname='salinity')
    assert "unexpected number of dimensions in requested" in str(e.value)
    dfmt.plot_ztdata(data_xr_sel=ds_sel, varname='salinity')
    dfmt.plot_ztdata(data_xr_sel=ds_sel, varname='salinity', only_contour=True)

