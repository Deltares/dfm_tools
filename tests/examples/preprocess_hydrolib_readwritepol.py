# -*- coding: utf-8 -*-
"""
Created on Wed Aug 24 13:25:41 2022

@author: veenstra
"""

import os
import datetime as dt
import matplotlib.pyplot as plt
plt.close('all')
import dfm_tools as dfmt
import hydrolib.core.dflowfm as hcdfm

dir_testinput = r'c:\DATA\dfm_tools_testdata'
dir_output = '.'

dtstart = dt.datetime.now()

file_pli_list = [r'p:\archivedprojects\11208054-004-dcsm-fm\models\model_input\bnd_cond\pli\DCSM-FM_OB_all_20181108.pli',
                 os.path.join(dir_testinput,'GSHHS_high_min1000km2.ldb'),
                 os.path.join(dir_testinput,'DFM_grevelingen_3D\\Grevelingen-FM_BL_fxw.pliz'),
                 os.path.join(dir_testinput,r'Maeslant.tek'),
                 os.path.join(dir_testinput,r'nima-1013-lo-wl_validname.tek'),
                 ]

for file_pli in file_pli_list:
    #load pol/tek/pli/ldb file
    polyfile_object = hcdfm.PolyFile(file_pli)
    print(f'processing PolyFile: {os.path.basename(file_pli)}')
    
    #geopandas
    if '.tek' in file_pli:
        polyobject_pd = dfmt.tekalobject_to_DataFrame(polyfile_object.objects[0])
        fig,ax = plt.subplots()
        polyobject_pd.plot(ax=ax)
        ax.set_title(f'1 PolyObjects, name of last is {polyobject_pd.index.name}')
    
    else:
        #df_polyline = dfmt.pointlike_to_DataFrame(polyfile_object.objects[0])
        #gdf_polyline = dfmt.pointlike_to_geodataframe_points(polyfile_object.objects[0],crs=crs)
        gdf_polyfile = dfmt.PolyFile_to_geodataframe_linestrings(polyfile_object,crs=None) #TODO: z-column (and next ones) are missing
        fig,ax = plt.subplots()
        gdf_polyfile.plot(ax=ax)
        pf_name = gdf_polyfile.iloc[-1]['name']
        ax.set_title(f'{len(gdf_polyfile)} PolyObjects, name of last is {pf_name}')
        
        #get extents of all objects in polyfile
        pol_bounds = gdf_polyfile.geometry.bounds
        print(pol_bounds['minx'].min(),pol_bounds['maxx'].max(),pol_bounds['miny'].min(),pol_bounds['maxy'].max())

time_passed = (dt.datetime.now()-dtstart).total_seconds()
print(f'>>time passed: {time_passed:.2f} sec')
