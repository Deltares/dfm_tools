# -*- coding: utf-8 -*-
"""
Created on Sun Jul  9 22:35:28 2023

@author: veenstra
"""

import pytest
import geopandas as gpd
import matplotlib.pyplot as plt
from shapely import Polygon, LineString
import dfm_tools as dfmt


@pytest.mark.unittest
def test_get_coastlines_gdb_global():
    """
    opens global coarse coastlines at lower 30 degrees, so includes L6/Antarctic
    """
    min_area = 1000 # km2
    bbox = (-180, -90, 180, -60) # xmin, ymin, xmax, ymax
    res = 'c'
    gdf = dfmt.get_coastlines_gdb(res=res, min_area=min_area, bbox=bbox)
    
    assert len(gdf) == 17


@pytest.mark.unittest
def test_get_coastlines_gdb_allres():
    """
    opens coastline gdb of Ireland and UK, only islands > 1000km2
    """
    min_area = 1000 # km2
    bbox = (-10, 55, 0, 60) # xmin, ymin, xmax, ymax
    for res,ngeoms in zip('fhilc',[5,5,5,5,4]):
        gdf = dfmt.get_coastlines_gdb(res=res, min_area=min_area, bbox=bbox)
        assert len(gdf) == ngeoms
    assert isinstance(gdf,gpd.GeoDataFrame)
    assert isinstance(gdf.geometry[0],Polygon)
    assert 'area' in gdf.columns
    assert 'geometry' in gdf.columns


@pytest.mark.unittest
def test_get_coastlines_gdb_lakes_islandsinlakes():
    """
    opens coastline gdb of lake Kivu (so lake and island in lake, L2 and L3)
    """
    
    bbox = (28.9, -2.4, 29.3, -1.6) # xmin, ymin, xmax, ymax
    gdf = dfmt.get_coastlines_gdb(res='h', min_area=0, bbox=bbox)
    
    assert len(gdf) == 9


@pytest.mark.unittest
def test_plot_coastlines():
    
    fig, ax = plt.subplots()
    ax.set_xlim(-10,0)
    ax.set_ylim(50,60)
    dfmt.plot_coastlines(res='l', min_area=500)


@pytest.mark.unittest
def test_get_borders_gdb_global():
    """
    opens global coarse borders
    """
    bbox = (-180, -90, 180, 90) # xmin, ymin, xmax, ymax
    res = 'c'
    gdf = dfmt.get_borders_gdb(res=res, bbox=bbox)
    
    assert len(gdf) == 451


@pytest.mark.unittest
def test_get_borders_gdb_allres():
    """
    opens borders gdb of Ireland and UK
    """
    bbox = (-10, 40, 5, 60) # xmin, ymin, xmax, ymax
    for res in 'fhilc':
        gdf = dfmt.get_borders_gdb(res=res, bbox=bbox)
        assert len(gdf) == 9
    assert isinstance(gdf,gpd.GeoDataFrame)
    assert isinstance(gdf.geometry[0],LineString)
    assert 'geometry' in gdf.columns


@pytest.mark.unittest
def test_plot_borders():
    
    fig, ax = plt.subplots()
    ax.set_xlim(-10,5)
    ax.set_ylim(40,60)
    dfmt.plot_borders(res='l')
