# -*- coding: utf-8 -*-
"""
Created on Fri Oct 20 16:09:20 2023

@author: veenstra
"""

import os
import pytest
import numpy as np
import pandas as pd
import geopandas as gpd
import dfm_tools as dfmt
import hydrolib.core.dflowfm as hcdfm
from dfm_tools.hydrolib_helpers import get_ncbnd_construct


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
def test_geodataframe_to_PolyFile_name_default(bnd_gdf):
    polyfile_obj = dfmt.geodataframe_to_PolyFile(bnd_gdf)
    names = [x.metadata.name for x in polyfile_obj.objects]
    assert names == ['L1', 'L2']


@pytest.mark.unittest
def test_geodataframe_to_PolyFile_name_some(bnd_gdf):
    polyfile_obj = dfmt.geodataframe_to_PolyFile(bnd_gdf, name="test_model")
    names = [x.metadata.name for x in polyfile_obj.objects]
    assert names == ['test_model1', 'test_model2']


@pytest.mark.unittest
def test_geodataframe_to_PolyFile_name_invalidtype(bnd_gdf):
    with pytest.raises(TypeError) as e:
        dfmt.geodataframe_to_PolyFile(bnd_gdf, name=None)
    assert 'name should be a string' in str(e.value)


@pytest.mark.unittest
def test_geodataframe_to_PolyFile_name_incorrect(bnd_gdf):
    with pytest.raises(ValueError) as e:
        dfmt.geodataframe_to_PolyFile(bnd_gdf, name='1')
    assert 'names in polyfile do not all start with a letter' in str(e.value)
    
    with pytest.raises(ValueError) as e:
        dfmt.geodataframe_to_PolyFile(bnd_gdf, name='-')
    assert 'names in polyfile do not all start with a letter' in str(e.value)
    
    with pytest.raises(ValueError) as e:
        dfmt.geodataframe_to_PolyFile(bnd_gdf, name='')
    assert 'names in polyfile do not all start with a letter' in str(e.value)


@pytest.mark.unittest
def test_geodataframe_to_PolyFile_namecolumn_some(bnd_gdf):
    bnd_gdf['name'] = ['test_model1','test_model2']
    polyfile_obj = dfmt.geodataframe_to_PolyFile(bnd_gdf)
    names = [x.metadata.name for x in polyfile_obj.objects]
    assert names == ['test_model1', 'test_model2']


@pytest.mark.unittest
def test_geodataframe_to_PolyFile_namecolumn_name_both(bnd_gdf):
    # name argument is ignored if name column is provided
    # not per se desired, but also not completely wrong
    bnd_gdf['name'] = ['test_model1','test_model2']
    polyfile_obj = dfmt.geodataframe_to_PolyFile(bnd_gdf, name='dummy')
    names = [x.metadata.name for x in polyfile_obj.objects]
    assert names == ['test_model1', 'test_model2']


@pytest.mark.unittest
def test_geodataframe_to_PolyFile_namecolumn_duplicated_names(bnd_gdf):
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
def test_geodataframe_to_PolyFile_namecolumn_numeric_start(bnd_gdf):
    bnd_gdf['name'] = ['1','2']
    with pytest.raises(ValueError) as e:
        dfmt.geodataframe_to_PolyFile(bnd_gdf)
    assert 'names in polyfile do not all start with a letter' in str(e.value)


@pytest.mark.unittest
def test_read_polyfile_as_gdf_points_linestrings(tmp_path):
    ncbnd_construct = get_ncbnd_construct()
    varn_pointname = ncbnd_construct['varn_pointname']
    
    # write a dummy polyfile with first three points of reference file
    # p:\archivedprojects\11208054-004-dcsm-fm\models\model_input\bnd_cond\pli\DCSM-FM_OB_all_20181108.pli
    poly_str = """DCSM-FM_OB_all_20181108
    3    2
    -9.25  43.0
    -9.50  43.0
    -9.75  43.0
    """
    file_pli = os.path.join(tmp_path, "dummy.pli")
    with open(file_pli, "w") as f:
        f.write(poly_str)
    
    # read polyfile with hydrolib-core
    polyfile_object = hcdfm.PolyFile(file_pli)
    
    # convert to geodataframe with points
    gdf_points = dfmt.PolyFile_to_geodataframe_points(polyfile_object)
    
    reference_names = ['DCSM-FM_OB_all_20181108_0001', 'DCSM-FM_OB_all_20181108_0002', 'DCSM-FM_OB_all_20181108_0003']
    reference_x = np.array([-9.25, -9.5 , -9.75])
    reference_y = np.array([43., 43., 43.])
    
    assert isinstance(gdf_points, gpd.GeoDataFrame)
    assert len(gdf_points) == 3
    assert (gdf_points.geometry.type == "Point").all()
    assert np.allclose(gdf_points.geometry.x.to_numpy(), reference_x)
    assert np.allclose(gdf_points.geometry.y.to_numpy(), reference_y)
    assert varn_pointname in gdf_points.columns
    assert gdf_points[varn_pointname].tolist() == reference_names

    # convert to geodataframe with linestrings
    gdf_lines = dfmt.PolyFile_to_geodataframe_linestrings(polyfile_object)
    line0_geom = gdf_lines.geometry.iloc[0]
    
    reference_names = ['DCSM-FM_OB_all_20181108']
    reference_x = np.array([-9.25, -9.5 , -9.75])
    reference_y = np.array([43., 43., 43.])
    
    assert isinstance(gdf_lines, gpd.GeoDataFrame)
    assert len(gdf_lines) == 1
    assert (gdf_lines.geometry.type == "LineString").all()
    assert np.allclose(line0_geom.xy[0], reference_x)
    assert np.allclose(line0_geom.xy[1], reference_y)
    assert 'name' in gdf_lines.columns
    assert gdf_lines['name'].tolist() == reference_names


@pytest.mark.unittest
def test_timmodel_to_dataframe(tmp_path):
    # original data from p:\archivedprojects\11206811-002-d-hydro-grevelingen\simulaties\Jaarsom2017_dfm_006_zlayer\boundary_conditions\hist\jaarsom_2017\sources_sinks\FlakkeeseSpuisluis.tim
    file_tim = os.path.join(tmp_path, "FlakkeeseSpuisluis_dummy.tim")
    str_tim = """* COLUMN_N=25 
     * COLUMN1=Time (min) since the reference date (2016-01-01)
     * COLUMN2=Flow discharge (m3/s), positive in
     * COLUMN3=Salinity (psu)
     * COLUMN4=Temperature (deg)
    * Time(min)   Flow         Salinity       Temperature
    0.00000e+00   0.0000e+00   2.883150e+01   8.4000e+00
    6.56639e+05   0.0000e+00   2.826540e+01   1.0295e+01
    6.56640e+05  -2.1468e+02   2.825750e+01   1.0300e+01
    6.56650e+05  -2.1384e+02   2.808100e+01   1.0400e+01
    6.56660e+05  -2.1384e+02   2.814450e+01   1.0350e+01
    6.56670e+05  -2.1128e+02   2.789250e+01   1.0300e+01
    6.56680e+05  -2.0955e+02   2.788500e+01   1.0350e+01
    6.56690e+05  -2.0781e+02   2.791250e+01   1.0250e+01
    """
    with open(file_tim, 'w') as f:
        f.write(str_tim)
    
    data_tim = hcdfm.TimModel(file_tim)
    
    refdate = '2016-01-01'
    tim_pd = dfmt.TimModel_to_DataFrame(data_tim, parse_column_labels=True, refdate=refdate)
    
    assert tim_pd.index[0] == pd.Timestamp(refdate)
    assert len(tim_pd) == 8
    assert tim_pd.columns[-1] == 'Temperature (deg)'
    
    data_tim_out = dfmt.DataFrame_to_TimModel(tim_pd, refdate=refdate)
    assert len(data_tim_out.comments) == 4
    assert len(data_tim_out.timeseries) == 8
