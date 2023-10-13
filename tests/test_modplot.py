# -*- coding: utf-8 -*-
"""
Created on Fri Oct 13 10:27:55 2023

@author: veenstra
"""

import pytest
import numpy as np
import matplotlib.pyplot as plt
plt.close('all')
import dfm_tools as dfmt


@pytest.mark.unittest
def test_velovect():
    """
    this test will fail with matplotlib<3.4.0
    """
    x = np.linspace(-4,4,120)
    y = np.linspace(-3,3,100)
    X,Y = np.meshgrid(x,y)
    U = -1 - X**2 + Y
    V = 1 + X - Y**2
    speed = np.sqrt(U*U + V*V)
    grains = 15
    
    # dfmt.velovect requires matplotlib>=3.4.0
    fig,ax = plt.subplots()
    dfmt.velovect(ax, X, Y, U, V, color=speed, cmap='winter', arrowstyle='fancy', 
                  linewidth=speed/5, integration_direction='forward',
                  density=5, grains=grains)


@pytest.mark.unittest
def test_broken_streamlines():
    """
    this test will fail with matplotlib<3.6.0
    """
    x = np.linspace(-4,4,120)
    y = np.linspace(-3,3,100)
    X,Y = np.meshgrid(x,y)
    U = -1 - X**2 + Y
    V = 1 + X - Y**2
    speed = np.sqrt(U*U + V*V)
    
    # broken_streamlines requires matplotlib>=3.6.0
    fig,ax = plt.subplots()
    ax.streamplot(X, Y, U, V, color=speed, cmap='winter', arrowstyle='fancy', linewidth=speed/5, broken_streamlines=False)
    