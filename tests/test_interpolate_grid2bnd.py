# -*- coding: utf-8 -*-
"""
Created on Mon Jul 17 13:08:34 2023

@author: veenstra
"""

import pytest
import dfm_tools as dfmt
import numpy as np


@pytest.mark.systemtest
@pytest.mark.requireslocaldata #TODO: this is not necessary in case of tpxo, but plifile is local
def test_interpolate_tide_to_plipoints_tpxo():
    
    nPoints = 3# None #amount of Points to process per PolyObject in the plifile (use int for testing, use None for all Points)
    file_pli = r'p:\archivedprojects\11208054-004-dcsm-fm\models\model_input\bnd_cond\pli\DCSM-FM_OB_all_20181108.pli'
    
    tidemodel = 'tpxo80' #FES2014, FES2012, EOT20, GTSM4.1preliminary
    for component_list in [None, ['M2','S2','M4']]:
        data_interp = dfmt.interpolate_tide_to_plipoints(tidemodel=tidemodel, file_pli=file_pli, component_list=component_list, nPoints=nPoints)
        
        compnames_now = data_interp['compnames'].to_numpy().tolist()
        if component_list is None:
            compnames_expected = ['M2','S2','N2','K2','K1','O1','P1','Q1','MF','MM','M4','MS4','MN4']
        else:
            compnames_expected = component_list
        
        data_interp_M2 = data_interp.sel(compno='M2')
        
        amp_expected = np.array([1.09643936, 1.08739412, 1.08555067])
        amp_now = data_interp_M2['amplitude'].to_numpy().tolist()
        
        phs_expected = np.array([81.92059326171875, 82.513671875, 82.69258117675781])
        phs_now = data_interp_M2['phase'].to_numpy().tolist()
        
        assert compnames_now == compnames_expected
        assert (np.abs(amp_expected-amp_now)<1e-6).all()
        assert (np.abs(phs_expected-phs_now)<1e-6).all()
        