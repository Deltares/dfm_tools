# -*- coding: utf-8 -*-
"""
Created on Mon Apr 24 14:29:19 2023

@author: veenstra
"""

import hydrolib.core.dflowfm as hcdfm
import geopandas
import matplotlib.pyplot as plt
plt.close('all')
import dfm_tools as dfmt
from shapely.geometry import Point, Polygon
import pandas as pd

dir_testinput = r'c:\DATA\dfm_tools_testdata'

def pointlike_to_geodataframe(polyline_object, crs='EPSG:4326', add_pointnames=True):
    #conversion to dataframe
    polyobject_pd = dfmt.pointlike_to_DataFrame(polyline_object)
    plipoints_list = []
    for i in range(0,len(polyobject_pd)):
        plipoints_list.append(Point(polyobject_pd['x'][i], polyobject_pd['y'][i]))
    
    #make gdf of points (1 point per row)
    gdf = geopandas.GeoDataFrame({'geometry': plipoints_list}, crs=crs)
    if add_pointnames:
        gdf['pointnames'] = pd.Series(polyobject_pd.index).apply(lambda x: f'{polyline_object.metadata.name}_{x+1:04d}')
    
    return gdf


def PolyFile_to_geodataframe(polyfile_object, crs='EPSG:4326'):
    plilines_list = []
    plinames_list = []
    #fig,ax = plt.subplots()
    for iPO, polyline_object in enumerate(polyfile_object.objects):
        print(f'processing PolyObject {iPO+1} of {len(polyfile_object.objects)}: name={polyline_object.metadata.name}')
        gdf = pointlike_to_geodataframe(polyline_object, crs=crs, add_pointnames=True)
        
        #make gdf of points (1 point per row)
        plilines_list.append(Polygon(gdf.geometry))
        plinames_list.append(polyline_object.metadata.name)
    gdf_polyfile = geopandas.GeoDataFrame({'name': plinames_list, 'geometry': plilines_list}, crs=crs)
    return gdf_polyfile


def geodataframerow_to_geodataframe(gdf_row, add_pointnames=True):
    x,y = gdf_row.geometry.exterior.coords.xy
    plipoints_list = []
    for xone,yone in zip(x.tolist(),y.tolist()):
        plipoints_list.append(Point(xone,yone))
    gdf = geopandas.GeoDataFrame({'geometry': plipoints_list}, crs=crs)
    if add_pointnames:
        name = gdf_row['name']
        gdf['pointnames'] = pd.Series(gdf.index).apply(lambda x: f'{name}_{x+1:04d}')
    
    return gdf


file_pli = r'p:\11208054-004-dcsm-fm\models\model_input\bnd_cond\pli\DCSM-FM_OB_all_20181108.pli'
                 # os.path.join(dir_testinput,'GSHHS_high_min1000km2.ldb')
                 # os.path.join(dir_testinput,'DFM_grevelingen_3D\\Grevelingen-FM_BL_fxw.pliz')
                 # os.path.join(dir_testinput,r'Maeslant.tek')
                 # os.path.join(dir_testinput,r'nima-1013-lo-wl_validname.tek')
                 
#load pol/tek/pli/ldb file
polyfile_object = hcdfm.PolyFile(file_pli)

crs = 'epsg:28992'
gdf_polyline = pointlike_to_geodataframe(polyfile_object.objects[0],crs=crs)
gdf_polyfile = PolyFile_to_geodataframe(polyfile_object,crs=crs)

#gdf_polyfile = gdf_polyfile.to_crs(4326) # convert to decimal degrees


gdf_polyline_frompolyfile = geodataframerow_to_geodataframe(gdf_row=gdf_polyfile.iloc[0])






