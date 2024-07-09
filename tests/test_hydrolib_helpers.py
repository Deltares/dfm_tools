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


@pytest.fixture(scope='session')
def bnd_gdf():
    dxy = 0.02
    crs = 4326
    lon_min, lon_max, lat_min, lat_max = -68.31, -68.27, 12.10, 12.21
    mk_object = dfmt.make_basegrid(lon_min, lon_max, lat_min, lat_max, dx=dxy, dy=dxy, crs=crs)
    bnd_gdf = dfmt.generate_bndpli_cutland(mk=mk_object, res='h', buffer=0.01)
    return bnd_gdf


@pytest.mark.unittest
def test_geodataframe_to_PolyFile_name_none(bnd_gdf):
    polyfile_obj = dfmt.geodataframe_to_PolyFile(bnd_gdf)
    names = [x.metadata.name for x in polyfile_obj.objects]
    assert names == ['L1', 'L2']


@pytest.mark.unittest
def test_geodataframe_to_PolyFile_name_some(bnd_gdf):
    polyfile_obj = dfmt.geodataframe_to_PolyFile(bnd_gdf, name="test_model")
    names = [x.metadata.name for x in polyfile_obj.objects]
    assert names == ['test_model1', 'test_model2']


@pytest.mark.unittest
def test_geodataframe_to_PolyFile_duplicated_names(bnd_gdf):
    # deliberately giving all polylines the same name
    bnd_gdf['name'] = 'duplicate_bnd'
    with pytest.raises(ValueError) as e:
        dfmt.geodataframe_to_PolyFile(bnd_gdf)
    assert 'duplicate polyline names found in polyfile' in str(e.value)

    # deliberately giving all polylines the same empty name
    bnd_gdf['name'] = ''
    with pytest.raises(ValueError) as e:
        dfmt.geodataframe_to_PolyFile(bnd_gdf)
    assert "duplicate polyline names found in polyfile" in str(e.value)


@pytest.mark.unittest
def test_geodataframe_to_PolyFile_incorrect_name(bnd_gdf):
    with pytest.raises(ValueError) as e:
        dfmt.geodataframe_to_PolyFile(bnd_gdf, name='1')
    assert 'name should start with a letter' in str(e.value)
    
    with pytest.raises(ValueError) as e:
        dfmt.geodataframe_to_PolyFile(bnd_gdf, name='-')
    assert 'name should start with a letter' in str(e.value)
    
    with pytest.raises(ValueError) as e:
        dfmt.geodataframe_to_PolyFile(bnd_gdf, name='')
    assert 'name is not allowed to be an empty string' in str(e.value)
