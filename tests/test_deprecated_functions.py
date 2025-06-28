# -*- coding: utf-8 -*-
"""
Created on Sat Jun 28 09:28:52 2025

@author: veenstra
"""

import dfm_tools as dfmt
import pytest

def test_deprecated_functions():
    for func in dir(dfmt.deprecated_functions):
        if func.startswith("__"):
            continue
        with pytest.raises(DeprecationWarning):
            exec(f"dfmt.{func}()")
