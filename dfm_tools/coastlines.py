# -*- coding: utf-8 -*-
"""
Created on Mon Apr 24 03:46:44 2023

@author: veenstra
"""

import os
import geopandas
import pandas as pd


def get_coastlines_gdb(res:str='h', bbox:tuple = (-180, -90, 180, 90), min_area:float = 0, crs:str = None, include_fields:list = ['area']) -> geopandas.geoseries.GeoSeries:
    """
    get coastlines from GSHHS database: https://www.ngdc.noaa.gov/mgg/shorelines/data/gshhg/latest/readme.txt
    geopandas docs https://geopandas.org/en/stable/docs/reference/api/geopandas.read_file.html
    
    Parameters
    ----------
    res : str, optional
        f(ull), h(igh), i(ntermediate), l(ow), c(oarse) resolution. The default is 'h'.
    bbox : tuple, optional
        (minx, miny, maxx, maxy), also includes shapes that are partly in the bbox. The default is (-180, -90, 180, 90).
    min_area : float, optional
        in km2, min_area>0 speeds up process. The default is 0.
    crs : str, optional
        
    include_fields : list, optional
        which shapefile fields to include, None gives all. The default is ['area'].

    Returns
    -------
    coastlines_gdb : TYPE
        DESCRIPTION.

    """
    #TODO: maybe download+unzip+load+cache from https://www.ngdc.noaa.gov/mgg/shorelines/
    #TODO: add auto deciding of res/min_area based on bbox size (in WGS84 coords)
    #TODO: maybe also add for rivers/borders?
    
    if res not in 'fhilc':
        raise KeyError(f'invalid res="{res}", resolution options are f(ull), h(igh), i(ntermediate), l(ow), c(oarse)')
    if crs is not None: #convert bbox from input crs to WGS84
        bbox_points = geopandas.points_from_xy(x=[bbox[0],bbox[2]], y=[bbox[1],bbox[3]], crs=crs)
        bbox_points = bbox_points.to_crs('EPSG:4326') #convert to WGS84
        bbox = (bbox_points.x[0], bbox_points.y[0], bbox_points.x[1], bbox_points.y[1])
        
    #L(evel)1 (coastlines except antarctica) and L6 (antarctic grounding line)
    dir_GSHHS_shp = r'p:\1230882-emodnet_hrsm\data\landboundary_GSHHS\gshhg-shp-2.3.7\GSHHS_shp'
    file_shp_L1 = os.path.join(dir_GSHHS_shp,res,f'GSHHS_{res}_L1.shp') #coastlines
    file_shp_L6 = os.path.join(dir_GSHHS_shp,res,f'GSHHS_{res}_L6.shp') #Antarctic grounding-line polygons
    file_shp_L2 = os.path.join(dir_GSHHS_shp,res,f'GSHHS_{res}_L2.shp') #lakes
    file_shp_L3 = os.path.join(dir_GSHHS_shp,res,f'GSHHS_{res}_L3.shp') #islands-in-lakes
    
    coastlines_gdb_L1 = geopandas.read_file(file_shp_L1, include_fields=include_fields, where=f"area>{min_area}", bbox=bbox)
    coastlines_gdb_L6 = geopandas.read_file(file_shp_L6, include_fields=include_fields, where=f"area>{min_area}", bbox=bbox)
    coastlines_gdb_list = [coastlines_gdb_L1,coastlines_gdb_L6]
    if len(coastlines_gdb_L1)<2: #if max one L1 polygon is selected, automatically add lakes and islands-in-lakes
        coastlines_gdb_L2 = geopandas.read_file(file_shp_L2, include_fields=include_fields, where=f"area>{min_area}", bbox=bbox)
        coastlines_gdb_L3 = geopandas.read_file(file_shp_L3, include_fields=include_fields, where=f"area>{min_area}", bbox=bbox)
        coastlines_gdb_list = [coastlines_gdb_L1,coastlines_gdb_L6,coastlines_gdb_L2,coastlines_gdb_L3]
    coastlines_gdb = pd.concat(coastlines_gdb_list)
    
    if crs:
        coastlines_gdb = coastlines_gdb.to_crs(crs)
    
    return coastlines_gdb


def plot_coastlines(ax, res:str='h', min_area:float = 0, crs=None, **kwargs):
    """
    get coastlines with get_coastlines_gdb and bbox depending on axlims, plot on ax and set axlims back to original values
    """
    #TODO: if ax is GeoAxis, get crs from ax
    
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    bbox = (xlim[0], ylim[0], xlim[1], ylim[1])
    
    if 'edgecolor' not in kwargs:
        kwargs['edgecolor'] = 'k'
    if 'facecolor' not in kwargs:
        kwargs['facecolor'] = 'none'
    if 'linewidth' not in kwargs:
        kwargs['linewidth'] = 0.7
    
    coastlines_gdb = get_coastlines_gdb(bbox=bbox, res=res, min_area=min_area, crs=crs)
    if coastlines_gdb.empty:
        return
    
    coastlines_gdb.plot(ax=ax, **kwargs)
    
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)

