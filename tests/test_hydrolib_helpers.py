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
def test_geodataframe_with_LineString_to_PolyFile():
    """
    converting a geodataframe with Linestring geometries to hcdfm.PolyFile
    """
    
    file_pli = r'p:\archivedprojects\11208054-004-dcsm-fm\models\model_input\bnd_cond\pli\DCSM-FM_OB_all_20181108.pli'
    polyfile_object = hcdfm.PolyFile(file_pli)

    # gdf containing shapely.LineString geometries
    gdf_polyfile = dfmt.PolyFile_to_geodataframe_linestrings(polyfile_object,crs=None)
    
    polyfile = dfmt.geodataframe_to_PolyFile(gdf_polyfile)
    
    assert isinstance(polyfile, hcdfm.PolyFile)
