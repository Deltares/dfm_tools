# -*- coding: utf-8 -*-
"""
Created on Mon Jul 17 13:08:34 2023

@author: veenstra
"""

import pytest
import dfm_tools as dfmt
import numpy as np


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


# TODO: test preprocess_interpolate_nc_to_bc.py and compare output to references (make extra systemtest if needed)
# TODO: add gtsm (if fast)
# TODO: make faster by not loading entire multifile dataset in case of only 3 components
# TODO: for 0to360 model like fes, also add point just before or after 0 meridian
@pytest.mark.systemtest
@pytest.mark.requireslocaldata # TODO: this is not necessary in case of tpxo since data is on opendap, but plifile is local
def test_interpolate_tide_to_plipoints():
    
    nPoints = 3# None #amount of Points to process per PolyObject in the plifile (use int for testing, use None for all Points)
    file_pli = r'p:\archivedprojects\11208054-004-dcsm-fm\models\model_input\bnd_cond\pli\DCSM-FM_OB_all_20181108.pli'
    
    tidemodel_list = ['tpxo80','FES2014'] #tpxo80, FES2014, GTSM4.1preliminary
    for tidemodel in tidemodel_list:
        print(tidemodel)
        if tidemodel=='tpxo80':
            amp_expected = np.array([1.09643936, 1.08739412, 1.08555067])
            phs_expected = np.array([81.92059326171875, 82.513671875, 82.69258117675781])
            component_list_full = ['M2','S2','N2','K2','K1','O1','P1','Q1','MF','MM','M4','MS4','MN4']
        elif tidemodel=='FES2014':
            amp_expected = np.array([1.1147955656051636, 1.1004363298416138, 1.0878639221191406])
            phs_expected = np.array([81.8884048461914, 82.19799041748047, 82.39784240722656])
            component_list_full = [] #TODO: add for FES2014
        
        for component_list in [['M2','S2','M4']]:#, None]:
            data_interp = dfmt.interpolate_tide_to_plipoints(tidemodel=tidemodel, file_pli=file_pli, component_list=component_list, nPoints=nPoints)
            
            compnames_now = data_interp['compno'].to_numpy().tolist()
            if component_list is None:
                compnames_expected = component_list_full
            else:
                compnames_expected = component_list
            
            data_interp_M2 = data_interp.sel(compno='M2')
            
            amp_now = data_interp_M2['amplitude'].to_numpy().tolist()
            
            phs_now = data_interp_M2['phase'].to_numpy().tolist()
            
            assert compnames_now == compnames_expected
            assert (np.abs(amp_expected-amp_now)<1e-6).all()
            assert (np.abs(phs_expected-phs_now)<1e-6).all()


