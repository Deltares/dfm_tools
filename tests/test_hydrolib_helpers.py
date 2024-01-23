# -*- coding: utf-8 -*-
"""
Created on Fri Oct 20 16:09:20 2023

@author: veenstra
"""

import pytest
import dfm_tools as dfmt
import hydrolib.core.dflowfm as hcdfm


@pytest.mark.unittest
def test_geodataframe_with_Polygon_to_PolyFile():
    """
    converting a geodataframe with Polygon geometries to hcdfm.PolyFile
    """
    lon_min, lon_max, lat_min, lat_max = -6, 2, 48.5, 51.2
    bbox = (lon_min, lat_min, lon_max, lat_max)
    
    # gdf containing shapely.Polygon geometries
    coastlines_gdf = dfmt.get_coastlines_gdb(bbox=bbox, res='h')
    
    polyfile = dfmt.geodataframe_to_PolyFile(coastlines_gdf)
    assert isinstance(polyfile, hcdfm.PolyFile)


@pytest.mark.unittest
def test_geodataframe_with_LineString_to_PolyFile(tmp_path):
    """
    converting a geodataframe with Linestring geometries to hcdfm.PolyFile
    """
    # write polygon file
    file_pol = tmp_path / 'temp_coastlines.pol'
    lon_min, lon_max, lat_min, lat_max = -68.45, -68.1, 12, 12.35
    bbox = (lon_min, lat_min, lon_max, lat_max)
    coastlines_gdf = dfmt.get_coastlines_gdb(bbox=bbox, res='h')
    polyfile = dfmt.geodataframe_to_PolyFile(coastlines_gdf)
    polyfile.save(file_pol)
    
    # read polygon file
    polyfile_object = hcdfm.PolyFile(file_pol)

    # gdf containing shapely.LineString geometries
    gdf_polyfile = dfmt.PolyFile_to_geodataframe_linestrings(polyfile_object,crs=None)
    
    polyfile = dfmt.geodataframe_to_PolyFile(gdf_polyfile)
    assert isinstance(polyfile, hcdfm.PolyFile)
