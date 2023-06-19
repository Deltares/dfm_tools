# -*- coding: utf-8 -*-
"""
Created on Mon Jun 19 22:03:43 2023

@author: veenstra
"""

import xarray as xr

def test_xarray_pandas_resample():
    """
    this testcase fails with pandas 2.0.2
    this testcase succeeds with pandas 1.5.3
    """
    ds = xr.tutorial.load_dataset("air_temperature")
    ds.resample(time='D')

