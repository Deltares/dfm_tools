# -*- coding: utf-8 -*-
"""
Created on Mon Jul 17 13:08:34 2023

@author: veenstra
"""

import os
import pytest
import dfm_tools as dfmt
import numpy as np
import pandas as pd
import datetime as dt
import xarray as xr
import shapely
import geopandas as gpd
from dfm_tools.interpolate_grid2bnd import (tidemodel_componentlist,
                                            components_translate_upper,
                                            interp_regularnc_to_plipointsDataset,
                                            check_time_extent,
                                            ext_add_boundary_object_per_polyline,
                                            ds_apply_conventions,
                                            ds_apply_conversion_dict,
                                            open_prepare_dataset,
                                            _ds_sel_time_outside,
                                            )
from dfm_tools.hydrolib_helpers import get_ncbnd_construct
import hydrolib.core.dflowfm as hcdfm
from dfm_tools.errors import OutOfRangeError
import warnings


@pytest.mark.unittest
def test_check_time_extent(cmems_dataset_4times):
    ds = cmems_dataset_4times
    
    # prior to ds timerange
    try:
        check_time_extent(ds, tstart="2011-01-01", tstop="2011-01-01")
    except OutOfRangeError as e:
        assert "tstart" in str(e)
    
    # after ds timerange
    try:
        check_time_extent(ds, tstart="2021-01-01", tstop="2021-01-01")
    except OutOfRangeError as e:
        assert "tstop" in str(e)
    
    # in ds timerange
    check_time_extent(ds, tstart="2020-01-01", tstop="2020-01-01")
    
    # on ds timerange
    check_time_extent(ds, tstart="2019-12-31 12:00:00", tstop="2020-01-03 12:00:00")


@pytest.mark.unittest
def test_conversion_dict():
    """
    since notations of equations are sometimes updated, convenient to have this test
    that checks the conversion values that were once used.
    """
    dict_keys_waq = ['tracerbndOXY', 'tracerbndNO3', 'tracerbndPO4', 'tracerbndSi', 
                     'tracerbndPON1', 'tracerbndPOP1', 'tracerbndPOC1', 
                     'tracerbndDON', 'tracerbndDOP', 'tracerbndDOC', 'tracerbndOpal']
    
    conversion_expected = {'tracerbndOXY': 0.032,
     'tracerbndNO3': 0.014,
     'tracerbndPO4': 0.03097,
     'tracerbndSi': 0.02808,
     'tracerbndPON1': 0.004226415094339623,
     'tracerbndPOP1': 0.0005843396226415094,
     'tracerbndPOC1': 0.024,
     'tracerbndDON': 0.013693584905660377,
     'tracerbndDOP': 0.0005843396226415094,
     'tracerbndDOC': 0.11678671698113208,
     'tracerbndOpal': 0.0018252}
    
    conversion_dict = dfmt.get_conversion_dict()
    for quan in dict_keys_waq:
        assert np.abs(conversion_dict[quan]['conversion'] - conversion_expected[quan]) < 1e-9


@pytest.mark.unittest
def test_tidemodel_componentlist():
    comp_list = tidemodel_componentlist(tidemodel='FES2014', convention=False)
    comp_list_convention = tidemodel_componentlist(tidemodel='FES2014', convention=True)
    
    assert len(comp_list) == 34
    assert len(comp_list_convention) == 34
    assert comp_list[1] == 'eps2'
    assert comp_list_convention[1] == 'EPSILON2'


@pytest.mark.unittest
def test_components_translate_upper():
    comp_list = components_translate_upper(['m2','eps2','e2'])
    assert comp_list == ['M2','EPSILON2','EPSILON2']


@pytest.mark.unittest
@pytest.mark.requireslocaldata
def test_plipointsDataset_fews_accepted():
    """
    check if FEWS netcdf bnd export is correctly processed to a 
    hcdfm ForcingModel object including the necessary conversions
    """
    file_nc_fews = r'p:\dflowfm\maintenance\JIRA\06000-06999\06187\C01\salinity_DCSM-FM_OB_all.nc'
    data_interp = xr.open_dataset(file_nc_fews)
    
    #convert plipointsDataset to hydrolib ForcingModel
    forcingmodel_object = dfmt.plipointsDataset_to_ForcingModel(plipointsDataset=data_interp)
    
    forcing0 = forcingmodel_object.forcing[0]
    assert isinstance(forcingmodel_object, hcdfm.ForcingModel)
    assert isinstance(forcing0, hcdfm.T3D)
    assert forcing0.quantityunitpair[1].unit == 'ppt'
    
    # test whether so was renamed to salinitybnd
    assert forcing0.quantityunitpair[1].quantity == 'salinitybnd'


@pytest.mark.systemtest
def test_plipointsDataset_to_ForcingModel_drop_allnan_points():
    #construct polyfile gdf
    point_x = [-71.5, -71.5, -71.5, -71.5,]
    point_y = [12.5, 12.6, 12.7, 12.8]
    polyobject_pd = pd.DataFrame(dict(x=point_x, y=point_y))
    poly_object = dfmt.DataFrame_to_PolyObject(polyobject_pd, name="abc_bnd")
    polyfile_object = hcdfm.PolyFile()
    polyfile_object.objects.append(poly_object)
    gdf_points = dfmt.PolyFile_to_geodataframe_points(polyfile_object)
    
    # actual cmems data
    no3_values = [[[       np.nan,        np.nan,        np.nan,        np.nan,
             5.2622618e-05],
            [       np.nan,        np.nan,        np.nan,        np.nan,
             5.1448391e-05],
            [9.0356334e-05, 1.8063841e-04,        np.nan, 7.0045113e-05,
             4.9546972e-05],
            [4.0657567e-05, 4.6050434e-05, 4.7762936e-05, 4.2734890e-05,
             4.0186114e-05],
            [2.6492798e-05, 2.8967228e-05, 3.0298394e-05, 3.1432221e-05,
             3.2138258e-05],
            [2.4448847e-05, 2.6705173e-05, 2.6929094e-05, 2.5548252e-05,
             2.5827681e-05]]]
    
    ds = xr.Dataset()
    ds['longitude'] = xr.DataArray([-72.  , -71.75, -71.5 , -71.25, -71.  ], dims='longitude')
    ds['latitude'] = xr.DataArray([12.  , 12.25, 12.5 , 12.75, 13.  , 13.25], dims='latitude')
    ds['time'] = xr.DataArray([639204.],dims='time').assign_attrs({'units': 'hours since 1950-01-01'})
    ds['tracerbndNO3'] =  xr.DataArray(no3_values, dims=('time','latitude','longitude')).assign_attrs({'units': 'mmol m-3'})
    ds = xr.decode_cf(ds)
    
    # interpolate to points and convert to forcingmodel
    data_interp = interp_regularnc_to_plipointsDataset(data_xr_reg=ds, gdf_points=gdf_points)
    forcingmodel_object = dfmt.plipointsDataset_to_ForcingModel(plipointsDataset=data_interp)
    
    # check if the resulting forcingmodel does not have 4 but 2 points
    # the first two points were skipped because they were all nan
    assert len(forcingmodel_object.forcing) == 2
    assert forcingmodel_object.forcing[0].name == 'abc_bnd_0003'
    assert forcingmodel_object.forcing[1].name == 'abc_bnd_0004'


def test_ds_sel_time_outside(cmems_dataset_4times):
    ds = cmems_dataset_4times
    
    # exact outer bounds
    tstart = "2019-12-31 12:00"
    tstop = "2020-01-03 12:00"
    ds_sel = _ds_sel_time_outside(ds, tstart, tstop)
    assert ds_sel.time.values[0] <= pd.Timestamp(tstart)
    assert ds_sel.time.values[1] > pd.Timestamp(tstart)
    assert ds_sel.time.values[-2] < pd.Timestamp(tstop)
    assert ds_sel.time.values[-1] >= pd.Timestamp(tstop)
    assert len(ds_sel.time) == 4
    
    # exact inner bounds
    tstart = "2020-01-01 12:00"
    tstop = "2020-01-02 12:00"
    ds_sel = _ds_sel_time_outside(ds, tstart, tstop)
    assert ds_sel.time.values[0] <= pd.Timestamp(tstart)
    assert ds_sel.time.values[1] > pd.Timestamp(tstart)
    assert ds_sel.time.values[-2] < pd.Timestamp(tstop)
    assert ds_sel.time.values[-1] >= pd.Timestamp(tstop)
    assert len(ds_sel.time) == 2

    # inexact inner bounds
    tstart = "2020-01-01"
    tstop = "2020-01-03"
    ds_sel = _ds_sel_time_outside(ds, tstart, tstop)
    assert ds_sel.time.values[0] <= pd.Timestamp(tstart)
    assert ds_sel.time.values[1] > pd.Timestamp(tstart)
    assert ds_sel.time.values[-2] < pd.Timestamp(tstop)
    assert ds_sel.time.values[-1] >= pd.Timestamp(tstop)
    assert len(ds_sel.time) == 4
    
    # tstart out of bounds
    tstart = "2019-12-30 12:00"
    tstop = "2020-01-03 12:00"
    with pytest.raises(OutOfRangeError) as e:
        _ds_sel_time_outside(ds, tstart, tstop)
    assert "requested tstart 2019-12-30 12:00:00 outside" in str(e.value)
    
    # tstart out of bounds
    tstart = "2019-12-31 12:00"
    tstop = "2030-01-03 12:00"
    with pytest.raises(OutOfRangeError) as e:
        _ds_sel_time_outside(ds, tstart, tstop)
    assert "requested tstop 2030-01-03 12:00:00 outside" in str(e.value)


@pytest.mark.systemtest
def test_open_prepare_dataset_correctdepths(tmp_path, cmems_dataset_4times):
    """
    to validate open_prepare_dataset behaviour for depths, in the past the depth values got lost and replaced by depth idx
    """
    
    ds_moretime = cmems_dataset_4times
    file_nc = tmp_path / 'temp_cmems_dummydata.nc'
    ds_moretime.to_netcdf(file_nc)
    
    ds_moretime_import = open_prepare_dataset(dir_pattern=file_nc, quantity='salinitybnd', tstart='2020-01-01 12:00', tstop='2020-01-02 12:00')
    assert len(ds_moretime_import.time) == 2


@pytest.mark.unittest
def test_ds_apply_conventions(cmems_dataset_4times):
    # generate datset with depths defined positive down
    ds_moretime = cmems_dataset_4times
    ds_moretime['depth'] = -1 * ds_moretime['depth']
    ds_moretime['depth'].attrs['positive'] = 'down'
    ds_converted = ds_apply_conventions(data_xr=ds_moretime)
    
    ncbnd_construct = get_ncbnd_construct()
    varn_depth = ncbnd_construct['varn_depth']
    dimn_depth = ncbnd_construct['dimn_depth']
    assert varn_depth in ds_converted.variables
    assert dimn_depth in ds_converted.dims
    
    depth_expected = np.array([-0.494025, -1.541375, -2.645669, -3.819495, -5.078224])
    depth_actual = ds_converted[varn_depth].to_numpy()
    assert np.allclose(depth_actual, depth_expected)
    

@pytest.mark.unittest
def test_ds_apply_conversion_dict_rename(cmems_dataset_4times):
    conversion_dict = dfmt.get_conversion_dict()
    ds_moretime = cmems_dataset_4times
    ds_converted = ds_apply_conversion_dict(data_xr=ds_moretime, conversion_dict=conversion_dict, quantity='salinitybnd')
    assert 'so' in ds_moretime.data_vars
    assert 'salinitybnd' in ds_converted.data_vars
    assert np.allclose(ds_converted['salinitybnd'], ds_moretime['so'], equal_nan=True)


@pytest.mark.unittest
def test_ds_apply_conversion_dict_rename_and_factor(cmems_dataset_4times):
    conversion_dict = dfmt.get_conversion_dict()
    ds_moretime = cmems_dataset_4times
    ds_moretime = ds_moretime.rename_vars({'so':'o2'})
    ds_moretime['o2'] = ds_moretime['o2'].assign_attrs({'units':'dummy'})
    ds_converted = ds_apply_conversion_dict(data_xr=ds_moretime, conversion_dict=conversion_dict, quantity='tracerbndOXY')
    assert 'o2' in ds_moretime.data_vars
    assert 'tracerbndOXY' in ds_converted.data_vars
    assert np.allclose(ds_converted['tracerbndOXY'], ds_moretime['o2']*0.032, equal_nan=True)
    assert ds_converted['tracerbndOXY'].attrs['units'] == 'g/m3'


@pytest.mark.unittest
def test_open_prepare_dataset_slightly_different_latlons(tmp_path, cmems_dataset_4times):
    """
    to check whether an error is raised when trying to combine datasets with slightly 
    different coordinates: https://github.com/Deltares/dfm_tools/issues/574
    
    """
    ds1 = cmems_dataset_4times.isel(time=slice(None,2))
    ds2 = cmems_dataset_4times.isel(time=slice(2,None))
    
    # deliberately alter longitude coordinate slightly
    ds_lon = ds1.longitude.to_numpy().copy()
    ds_lon[1] += 1e-8
    ds2['longitude'] = xr.DataArray(ds_lon,dims='longitude')
    
    file_nc1 = tmp_path / 'temp_cmems_2day_p1.nc'
    file_nc2 = tmp_path / 'temp_cmems_2day_p2.nc'
    ds1.to_netcdf(file_nc1)
    ds2.to_netcdf(file_nc2)
    file_nc = tmp_path / 'temp_cmems_2day_*.nc'
    
    with pytest.raises(ValueError) as e:
        open_prepare_dataset(file_nc, quantity='salinitybnd', tstart='2020-01-01 12:00:00', tstop='2020-01-02 12:00:00')
    
    # ValueError: cannot align objects with join='exact' where index/labels/sizes are not equal along these coordinates (dimensions): 'longitude' ('longitude',)
    assert "cannot align objects with join='exact' " in str(e.value)


@pytest.mark.unittest
def test_interp_regularnc_to_plipointsDataset(cmems_dataset_notime):
    """
    Linear interpolation to a new dimension in dfmt.interp_regularnc_to_plipoints() 
    resulted in unexpected nan values since scipy 1.10.0.
    More info in https://github.com/pydata/xarray/issues/7701 and related issues.
    Since this interpn is way faster than interp1d on a grid or all separate xy points, 
    we sticked to interpn but filled nan values with a nearest interpolation. 
    This only happens if there are non-nan values in the surrounding cells.
    In this test we check whether this combination of interpolation results in the same values as interp1d.
    This is mainly relevant when interpolating to existing lat/lon values, so that is the case we test here.
    This method gives as much valid values as possible given the input dataset, but does not fill nans where it should not do that.
    """
    
    ncbnd_construct = get_ncbnd_construct()
    dimn_point = ncbnd_construct['dimn_point']
    dimn_depth = ncbnd_construct['dimn_depth']
    varn_depth = ncbnd_construct['varn_depth']
    varn_pointname = ncbnd_construct['varn_pointname']
    
    ds = cmems_dataset_notime.copy()
    ds = ds.rename_dims({'depth':dimn_depth})
    ds = ds.rename_vars({'depth':varn_depth})
    so_np = ds['so'].to_numpy()
    lons = ds['longitude'].to_numpy()
    lats = ds['latitude'].to_numpy()
    
    for ipoint in range(3):
        x_xr = xr.DataArray([lons[ipoint]],dims=(dimn_point))
        y_xr = xr.DataArray([lats[ipoint]],dims=(dimn_point))
        
        # interp1d # these are actually irrelevant now, are not used for testing
        interp_with_floats = ds.interp(longitude=x_xr[0], latitude=y_xr[0], method='linear').so #selecting one value from the da drops the new plipoints dimension
        interp_with_da_existing = ds.interp(longitude=x_xr.values, latitude=y_xr.values, method='linear').so.isel(longitude=0,latitude=0) #using the DataArray values keeps lat/lon dimenions, gives the same interp result
        
        geom = shapely.points(x_xr, y_xr)
        gdf = gpd.GeoDataFrame(data={varn_pointname:[f'name_{ipoint+1:04d}']}, geometry=geom)
        interp_with_da_newdim = dfmt.interp_regularnc_to_plipointsDataset(ds, gdf, load=True)
        interp_da_actual = interp_with_da_newdim.so.isel({dimn_point:0})
        
        #define expected values since in some cases like point 0 this is not the same as interp1d returns
        interp_da_expected = xr.DataArray(so_np[:,ipoint,ipoint],dims=(dimn_depth))
        
        print(ipoint)
        print(interp_with_floats.to_numpy())
        print(interp_with_da_existing.to_numpy())
        print(interp_da_actual.to_numpy())
        print(interp_da_expected.to_numpy())
        
        assert (interp_with_floats.isnull()==interp_with_da_existing.isnull()).all()
        assert (interp_da_expected.isnull()==interp_da_actual.isnull()).all()
    
    """
    prints with scipy 1.11.3:
    0
    [35.819576 35.802788 35.786003       nan       nan]
    [35.819576 35.802788 35.786003       nan       nan]
    [35.819576 35.802788 35.786003 35.776848 35.781425]
    [35.819576 35.802788 35.786003 35.776848 35.781425]
    1
    [35.824154 35.810417 35.796684       nan       nan]
    [35.824154 35.810417 35.796684       nan       nan]
    [35.824154 35.810417 35.796684       nan       nan]
    [35.824154 35.810417 35.796684       nan       nan]
    2
    [35.82873 35.81805      nan      nan      nan]
    [35.82873 35.81805      nan      nan      nan]
    [35.82873 35.81805      nan      nan      nan]
    [35.82873 35.81805      nan      nan      nan]
    """


@pytest.mark.unittest
def test_interp_regularnc_to_plipointsDataset_checkvardimnames(cmems_dataset_notime):
    ncbnd_construct = get_ncbnd_construct()
    dimn_point = ncbnd_construct['dimn_point']
    dimn_depth = ncbnd_construct['dimn_depth']
    varn_depth = ncbnd_construct['varn_depth']
    varn_pointx = ncbnd_construct['varn_pointx']
    varn_pointy = ncbnd_construct['varn_pointy']
    varn_pointname = ncbnd_construct['varn_pointname']
    
    ds = cmems_dataset_notime.copy()
    ds = ds.rename_dims({'depth':dimn_depth})
    ds = ds.rename_vars({'depth':varn_depth})
    lons = ds['longitude'].to_numpy()
    lats = ds['latitude'].to_numpy()
    
    x_xr = xr.DataArray([lons[0]],dims=(dimn_point))
    y_xr = xr.DataArray([lats[0]],dims=(dimn_point))
    
    geom = shapely.points(x_xr, y_xr)
    gdf = gpd.GeoDataFrame(data={varn_pointname:['name_0001']}, geometry=geom)
    interp_with_da_newdim = dfmt.interp_regularnc_to_plipointsDataset(ds, gdf, load=True)
    
    # check if only expected dims/vars are present
    varn_inda = set(list(interp_with_da_newdim.variables))
    varn_expected = set(['so', varn_depth, varn_pointx, varn_pointy, varn_pointname])
    dimn_inda = set(list(interp_with_da_newdim.dims))
    dimn_expected = set([dimn_point, dimn_depth])
    assert varn_inda == varn_expected
    assert dimn_inda == dimn_expected


@pytest.mark.systemtest
@pytest.mark.requireslocaldata
def test_interpolate_tide_to_plipoints(data_dcsm_gdf):
    nanvalue = -999
    gdf_points = data_dcsm_gdf
    
    tidemodel_list = ['tpxo80_opendap', 'FES2014', 'FES2012', 'EOT20', 'GTSMv4.1']#, 'GTSMv4.1_opendap']
    for tidemodel in tidemodel_list:
        print(tidemodel)
        dtstart = dt.datetime.now()
        component_list_tidemodel = tidemodel_componentlist(tidemodel, convention=True)
        
        if tidemodel=='tpxo80_opendap': # 4.7 sec (all components: 5.8 sec)
            amp_expected = np.array([1.09643936, 1.08739412, 1.08555067])
            phs_expected = np.array([81.92059326171875, 82.513671875, 82.69258117675781])
        elif tidemodel=='FES2014': # 10.5 sec (all components: 97.6 sec)
            amp_expected = np.array([1.1147955656051636, 1.1004363298416138, 1.0878639221191406])
            phs_expected = np.array([81.8884048461914, 82.19799041748047, 82.39784240722656])
        elif tidemodel=='FES2012': # 9.5 sec (all components: 92.4 sec)
            amp_expected = np.array([nanvalue, 1.0839321613311768, 1.0718189477920532])
            phs_expected = np.array([nanvalue, 82.33390808105469, 82.60922241210938])
        elif tidemodel=='EOT20': # 9.9 sec (all components: 61.9 sec)
            amp_expected = np.array([nanvalue, 1.10231938,  1.08968516])
            phs_expected = np.array([nanvalue, 82.54157149, 82.76990983])
        elif tidemodel in ['GTSMv4.1','GTSMv4.1_opendap']: # 13.8 sec vs 39.2 sec (all components: 82.4 sec vs 498.4 sec)
            amp_expected = np.array([1.13028932, 1.10648024, 1.09396541])
            phs_expected = np.array([81.21875763, 81.41669464, 81.66479492])
        
        for component_list in [['M2','S2']]: # [None]: # 
            data_interp = dfmt.interpolate_tide_to_plipoints(tidemodel=tidemodel, gdf_points=gdf_points, component_list=component_list)
            
            compnames_now = data_interp['compno'].to_numpy().tolist()
            if component_list is None:
                compnames_expected = component_list_tidemodel
            else:
                compnames_expected = component_list
            
            data_interp_M2 = data_interp.sel(compno='M2').fillna(nanvalue)
            
            amp_now = data_interp_M2['amplitude'].to_numpy()
            phs_now = data_interp_M2['phase'].to_numpy()
            
            assert compnames_now == compnames_expected
            assert (np.abs(amp_expected-amp_now)<1e-6).all()
            assert (np.abs(phs_expected-phs_now)<1e-6).all()
        
        print(f'>> tide interpolation from {tidemodel} took: ',end='')
        print(f'{(dt.datetime.now()-dtstart).total_seconds():.2f} sec')

    
@pytest.mark.unittest
def test_interpolate_tide_to_forcingmodel(data_dcsm_gdf):
    """
    This tests adds to test_interpolate_tide_to_plipoints, since it also interpolates to ForcingModel
    Furthermore, it runs on Github since it does not depend on local data
    """
    
    tidemodel = 'tpxo80_opendap'
    component_list = ['M2']
    gdf_points = data_dcsm_gdf

    data_interp = dfmt.interpolate_tide_to_plipoints(tidemodel=tidemodel, gdf_points=gdf_points, 
                                                     component_list=component_list, load=True)
    ForcingModel_object = dfmt.plipointsDataset_to_ForcingModel(plipointsDataset=data_interp)
    
    forcing0 = ForcingModel_object.forcing[0]
    assert isinstance(ForcingModel_object, hcdfm.ForcingModel)
    assert isinstance(forcing0, hcdfm.Astronomic)
    
    assert forcing0.quantityunitpair[0].unit == '-'
    assert forcing0.quantityunitpair[0].quantity == 'astronomic component'
    assert forcing0.quantityunitpair[1].unit == 'm'
    assert forcing0.quantityunitpair[1].quantity == 'waterlevelbnd amplitude'
    assert forcing0.quantityunitpair[2].unit == 'degrees'
    assert forcing0.quantityunitpair[2].quantity == 'waterlevelbnd phase'


@pytest.mark.unittest
def test_interp_uds_to_plipoints():
    """
    very basic test for function, 
    should be made more strict with learnings from postprocess_interpolate_uds_toplipoints.py
    """
    
    ncbnd_construct = get_ncbnd_construct()
    varn_pointname = ncbnd_construct['varn_pointname']
    
    file_nc = dfmt.data.fm_grevelingen_map(return_filepath=True)
    uds = dfmt.open_partitioned_dataset(file_nc)
    
    # adding xy=[1,1] deliberately to test if there are nans included
    gdf = gpd.GeoDataFrame(geometry=gpd.points_from_xy([51500,55000,1],[418800,421500,1]))
    gdf[varn_pointname] = ['pt_0001','pt_0002','pt_0003']
    
    # fig, ax = plt.subplots()
    # uds.mesh2d_flowelem_bl.ugrid.plot(ax=ax)
    # gdf.plot(ax=ax)
    
    # interpolate depths (zsigma to z):
    ds_atdepths = dfmt.get_Dataset_atdepths(data_xr=uds, depths=[-5,-1])
    ds_atdepths = ds_atdepths.rename({'depth_from_z0':'depth'})
    
    #interpolate to plipoints
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        ds_plipoints = dfmt.interp_uds_to_plipoints(uds=ds_atdepths, gdf=gdf) #workaround for plipoints out of the model domain
    assert varn_pointname in ds_plipoints.coords
    assert 'depth' in ds_plipoints.coords # TODO: maybe should be z like ncbnd_construct['varn_depth'], although get_Dataset_atdepths() makes depth variable
    
    retrieved = ds_plipoints.mesh2d_sa1.isel(time=-1).to_numpy()
    expected = np.array([[29.00111981, 28.99379263],
                         [28.95170139, 28.63716083],
                         [np.nan, np.nan]])
    
    assert np.allclose(retrieved, expected, equal_nan=True)


@pytest.mark.unittest
def test_ext_add_boundary_object_per_polyline_onepolyline(tmp_path):
    file_pli = os.path.join(tmp_path,'test_model.pli')
    with open(file_pli,'w') as f:
        f.write("""name1
        2    2
        1.0    2.0
        3.0    4.0
        """)
    
    # new ext
    ext_new = hcdfm.ExtModel()
    forcingmodel_object = hcdfm.ForcingModel()
    boundary_object = hcdfm.Boundary(quantity='waterlevelbnd',
                                     locationfile=file_pli,
                                     forcingfile=forcingmodel_object)
    ext_add_boundary_object_per_polyline(ext_new=ext_new, boundary_object=boundary_object)
    assert len(ext_new.boundary) == 1
    polyfile_names = [os.path.basename(x.locationfile.filepath) for x in ext_new.boundary]
    assert polyfile_names == ['test_model.pli']


@pytest.mark.unittest
def test_ext_add_boundary_object_per_polyline_twopolylines(tmp_path):
    file_pli = os.path.join(tmp_path,'test_model.pli')
    with open(file_pli,'w') as f:
        f.write("""name1
        2    2
        1.0    2.0
        3.0    4.0
        name2
        2    2
        1.0    2.0
        3.0    4.0
        """)
    
    # new ext
    ext_new = hcdfm.ExtModel()
    forcingmodel_object = hcdfm.ForcingModel()
    boundary_object = hcdfm.Boundary(quantity='waterlevelbnd',
                                     locationfile=file_pli,
                                     forcingfile=forcingmodel_object)
    ext_add_boundary_object_per_polyline(ext_new=ext_new, boundary_object=boundary_object)
    assert len(ext_new.boundary) == 2
    polyfile_names = [os.path.basename(x.locationfile.filepath) for x in ext_new.boundary]
    assert polyfile_names == ['name1.pli', 'name2.pli']


@pytest.mark.unittest
def test_ext_add_boundary_object_per_polyline_wrong_name(tmp_path):
    file_pli = os.path.join(tmp_path,'test_model.pli')
    with open(file_pli,'w') as f:
        f.write("""test_model
        2    2
        1.0    2.0
        3.0    4.0
        name2
        2    2
        1.0    2.0
        3.0    4.0
        """)
    
    # new ext
    ext_new = hcdfm.ExtModel()
    forcingmodel_object = hcdfm.ForcingModel()
    boundary_object = hcdfm.Boundary(quantity='waterlevelbnd',
                                     locationfile=file_pli,
                                     forcingfile=forcingmodel_object)
    with pytest.raises(ValueError) as e:
        ext_add_boundary_object_per_polyline(ext_new=ext_new, boundary_object=boundary_object)
    assert "The names of one of the polylines in the polyfile is the same as the polyfilename" in str(e.value)


@pytest.mark.unittest
def test_open_dataset_extra_deprecated():
    with pytest.raises(DeprecationWarning):
        dfmt.open_dataset_extra()


@pytest.mark.unittest
def test_interp_hisnc_to_plipoints(tmp_path):
    ncbnd_construct = get_ncbnd_construct()
    dimn_point = ncbnd_construct['dimn_point']
    
    file_pli = os.path.join(tmp_path,'grev_internal.pli')
    with open(file_pli,'w') as f:
        f.write("""grev_internal
        4    2
        51655.45  420059.29
        55133.26  421410.41
        59848.92  422508.19
        61676.24  417399.27
        """)

    ds_his = dfmt.data.fm_grevelingen_his()
    da = ds_his[['waterlevel']]

    ds_interp = dfmt.interp_hisnc_to_plipoints(data_xr_his=da, file_pli=file_pli, kdtree_k=4)
    assert ds_interp.sizes[dimn_point] == 4
    assert ds_interp.sizes["time"] == 145

    last_points = ds_interp.waterlevel.isel(time=-1).to_numpy()
    last_points_expected = np.array(
        [-0.20098736, -0.19833495, -0.19681664, -0.19467336]
        )
    assert np.allclose(last_points, last_points_expected)
    pointnames_np = ds_interp.station_id.to_numpy()
    pointnames = [str(x) for x in pointnames_np]
    pointnames_expected = [
        'grev_internal_0001',
        'grev_internal_0002',
        'grev_internal_0003',
        'grev_internal_0004',
        ]
    assert pointnames == pointnames_expected
    
    assert set(ds_interp.coords) == set(['lat', 'lon', 'station_id', 'time'])
