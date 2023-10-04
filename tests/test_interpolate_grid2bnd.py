# -*- coding: utf-8 -*-
"""
Created on Mon Jul 17 13:08:34 2023

@author: veenstra
"""

import pytest
import dfm_tools as dfmt
import numpy as np
import datetime as dt
import xarray as xr


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
    comp_list = dfmt.tidemodel_componentlist(tidemodel='FES2014', convention=False)
    comp_list_convention = dfmt.tidemodel_componentlist(tidemodel='FES2014', convention=True)
    
    assert len(comp_list) == 34
    assert len(comp_list_convention) == 34
    assert comp_list[1] == 'eps2'
    assert comp_list_convention[1] == 'EPSILON2'


@pytest.mark.unittest
def test_components_translate_upper():
    comp_list = dfmt.components_translate_upper(['m2','eps2','e2'])
    assert comp_list == ['M2','EPSILON2','EPSILON2']


@pytest.mark.unittest
def test_xarray_interp_to_newdim():
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
    
    ds = xr.Dataset()
    so_np = np.array([[[35.819576, 35.82568 , 35.82873 ],
                       [35.819576, 35.824154, 35.831783],
                       [35.822628, 35.824154, 35.82873 ]],
                      
                      [[35.802788, 35.80584 , 35.815   ],
                       [35.815   , 35.810417, 35.821102],
                       [35.824154, 35.813473, 35.81805 ]],
                      
                      [[35.786003, 35.789055, np.nan],
                       [35.807365, 35.796684, np.nan],
                       [35.824154, 35.80584 , np.nan]],
                      
                      [[35.776848, np.nan,    np.nan],
                       [35.792107, np.nan,    np.nan],
                       [35.822628, np.nan,    np.nan]],
                                              
                      [[35.781425, np.nan,    np.nan],
                       [35.792107, np.nan,    np.nan],
                       [35.789055, np.nan,    np.nan]]])
    ds['so'] = xr.DataArray(so_np,dims=('depth','latitude','longitude'))
    lons = [-9.6, -9.5, -9.4]
    lats = [42.9, 43.0, 43.1]
    ds['longitude'] = xr.DataArray(lons, dims=('longitude'))
    ds['latitude'] = xr.DataArray(lats, dims=('latitude'))
    
    for ipoint in range(3):
        x_xr = xr.DataArray([lons[ipoint]],dims=('plipoints'))
        y_xr = xr.DataArray([lats[ipoint]],dims=('plipoints'))
        
        # interp1d # these are actually irrelevant now, are not used for testing
        interp_with_floats = ds.interp(longitude=x_xr[0], latitude=y_xr[0], method='linear').so #selecting one value from the da drops the new plipoints dimension
        interp_with_da_existing = ds.interp(longitude=x_xr.values, latitude=y_xr.values, method='linear').so.isel(longitude=0,latitude=0) #using the DataArray values keeps lat/lon dimenions, gives the same interp result
        
        # combination of linear and nearest interpn
        interp_with_da_newdim_lin = ds.interp(longitude=x_xr, latitude=y_xr, method='linear').so.isel(plipoints=0) #using the DataArray introduces a plipoints dimension, which gives different interp result
        interp_with_da_newdim_near = ds.interp(longitude=x_xr, latitude=y_xr, method='nearest').so.isel(plipoints=0) #using the DataArray introduces a plipoints dimension, which gives different interp result
        interp_with_da_newdim = interp_with_da_newdim_lin.combine_first(interp_with_da_newdim_near)
        
        #define expected values since in some cases like point 0 this is not the same as interp1d returns
        interp_da_expected = xr.DataArray(so_np[:,ipoint,ipoint],dims=('depth'))
        
        print(ipoint)
        print(interp_with_floats.to_numpy())
        print(interp_with_da_existing.to_numpy())
        print(interp_with_da_newdim.to_numpy())
        print(interp_da_expected.to_numpy())
        
        assert (interp_with_floats.isnull()==interp_with_da_existing.isnull()).all()
        assert (interp_da_expected.isnull()==interp_with_da_newdim.isnull()).all()

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


@pytest.mark.systemtest
@pytest.mark.requireslocaldata
def test_interpolate_tide_to_plipoints():
    
    nPoints = 3# None #amount of Points to process per PolyObject in the plifile (use int for testing, use None for all Points)
    file_pli = r'p:\archivedprojects\11208054-004-dcsm-fm\models\model_input\bnd_cond\pli\DCSM-FM_OB_all_20181108.pli'
    nanvalue = -999
    
    tidemodel_list = ['tpxo80_opendap', 'FES2014', 'FES2012', 'EOT20', 'GTSMv4.1']#, 'GTSMv4.1_opendap']
    for tidemodel in tidemodel_list:
        print(tidemodel)
        dtstart = dt.datetime.now()
        component_list_tidemodel = dfmt.tidemodel_componentlist(tidemodel, convention=True)
        
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
        
        for component_list in [['M2','S2','M4']]: # [None]: # 
            data_interp = dfmt.interpolate_tide_to_plipoints(tidemodel=tidemodel, file_pli=file_pli, component_list=component_list, nPoints=nPoints)
            
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
            
            # file_bc_out = f'tide_{tidemodel}.bc'
            # ForcingModel_object = dfmt.plipointsDataset_to_ForcingModel(plipointsDataset=data_interp)
            # ForcingModel_object.save(filepath=file_bc_out)
    
        print(f'>> tide interpolation from {tidemodel} took: ',end='')
        print(f'{(dt.datetime.now()-dtstart).total_seconds():.2f} sec')

