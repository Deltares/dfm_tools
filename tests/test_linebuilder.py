# -*- coding: utf-8 -*-
"""
Created on Thu Jul  3 09:35:25 2025

@author: veenstra
"""
import dfm_tools as dfmt
import pytest
import matplotlib.pyplot as plt
from matplotlib.backend_bases import MouseEvent, MouseButton
import numpy as np

@pytest.mark.unittest
def test_linebuilder():
    """
    https://matplotlib.org/stable/api/backend_bases_api.html#matplotlib.backend_bases.LocationEvent
    """
    # adding point: x=-0.037480, y=0.030774
    # adding point: x=-0.008206, y=0.034345
    # adding point: x=-0.004435, y=0.000119
    # removing last point if present
    # adding point: x=0.009980, y=0.029881
    # adding point: x=0.021069, y=0.019464

    fig, ax = plt.subplots()
    linebuilder = dfmt.LineBuilder(ax=ax, block=False)
    
    for event in [
        MouseEvent(name="outside", canvas=fig.canvas, x=1, y=1, key="control", button=MouseButton.LEFT, dblclick=False),
        MouseEvent(name="loc0", canvas=fig.canvas, x=100, y=100, key="control", button=MouseButton.LEFT, dblclick=False),
        MouseEvent(name="loc1", canvas=fig.canvas, x=120, y=120, key="control", button=MouseButton.LEFT, dblclick=False),
        MouseEvent(name="loc2", canvas=fig.canvas, x=110, y=110, key="control", button=MouseButton.LEFT, dblclick=False),
        MouseEvent(name="mistake", canvas=fig.canvas, x=190, y=190, key="control", button=MouseButton.LEFT, dblclick=False),
        MouseEvent(name="undo", canvas=fig.canvas, x=190, y=190, key="control", button=MouseButton.RIGHT, dblclick=False),
        MouseEvent(name="loclast", canvas=fig.canvas, x=150, y=150, key="control", button=MouseButton.LEFT, dblclick=True)
    ]:
        linebuilder(event)
    
    expected_values = np.array(
        [[0.04032258, 0.12770563],
         [0.08064516, 0.18181818],
         [0.06048387, 0.1547619 ]]
        )
    # assert for the correct expected values
    assert np.allclose(linebuilder.line_array, expected_values)
