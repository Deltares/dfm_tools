# -*- coding: utf-8 -*-
"""
Created on Thu Jun 22 09:41:49 2023

@author: veenstra
"""

import pytest
import pandas as pd
import numpy as np
from dfm_tools.download import round_timestamp_to_outer_noon


@pytest.mark.unittest
def test_round_timestamp_to_outer_noon():
    date_min_mod_list = pd.date_range('2020-01-01','2020-01-02',freq='30min')
    date_max_mod_list = pd.date_range('2020-01-03','2020-01-04',freq='30min')
    
    for date_min_mod, date_max_mod in zip(date_min_mod_list, date_max_mod_list):
        
        date_min, date_max = round_timestamp_to_outer_noon(date_min_mod, date_max_mod)
        
        assert date_min.hour == 12
        assert date_min.minute == 0
        assert date_max.hour == 12
        assert date_max.minute == 0
        assert date_min <= pd.Timestamp(date_min_mod)
        assert date_max >= pd.Timestamp(date_max_mod)
        assert np.abs((date_min - pd.Timestamp(date_min_mod)).total_seconds()) < 24*3600
        assert np.abs((date_max - pd.Timestamp(date_max_mod)).total_seconds()) < 24*3600