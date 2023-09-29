# -*- coding: utf-8 -*-
"""
Created on Mon Jul 17 13:08:34 2023

@author: veenstra
"""

import pytest
import dfm_tools as dfmt
import numpy as np
import datetime as dt


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
    comp_list = dfmt.tidemodel_componentlist(tidemodel='FES2014', convention=True)
    assert len(comp_list) == 32
    assert comp_list[1] == 'EPSILON2'


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

