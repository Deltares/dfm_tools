# -*- coding: utf-8 -*-
"""
Created on Wed Oct 27 22:04:02 2021

@author: veenstra

"""

import os
import numpy as np
import matplotlib.pyplot as plt
plt.close('all')
import geopandas as gpd #conda install --channel conda-forge geopandas (breaks dfm_tools environment because of Qt issue)
from shapely.geometry import Polygon
import dfm_tools as dfmt
import datetime as dt

dtstart_script = dt.datetime.now()

dir_testinput = r'c:\DATA\dfm_tools_testdata'
dir_output = '.'


varlist = ['mesh2d_sa1','mesh2d_s1']#'Chlfa']#,'mesh2d_s1']
dir_shp = dir_output
if not os.path.exists(dir_shp):
    os.makedirs(dir_shp)
#file_nc = os.path.join(r'p:\archivedprojects\11203850-coastserv\06-Model\waq_model\simulations\run0_20200319\DFM_OUTPUT_kzn_waq', 'kzn_waq_0*_map.nc') #TODO: "ValueError: All NaN slice encountered" (is normally a runtimewarning, but not here)
file_nc = os.path.join(dir_testinput,'DFM_3D_z_Grevelingen','computations','run01','DFM_OUTPUT_Grevelingen-FM','Grevelingen-FM_0*_map.nc')

data_xr_map = dfmt.open_partitioned_dataset(file_nc)
vars_pd = dfmt.get_ncvarproperties(data_xr_map)

for timestep in [3]:#[0,10,20,30]:
    data_map_timesel = data_xr_map.isel(time=timestep)
    
    #data_sel = dfmt.get_mapdata_atdepths(data_xr_map=data_map_timesel, depths=0, reference='waterlevel') #top layer: 0m from waterlevel
    data_sel = dfmt.get_mapdata_atdepths(data_xr_map=data_map_timesel, depths=-5, reference='z0') #4m from model reference
    #data_sel = dfmt.get_mapdata_atdepths(data_xr_map=data_map_timesel, depths=0, reference='bedlevel') #bottomlayer: 0m from bedlevel
    
    print('creating geodataframe with cells from ugrid_verts')
    #converting ugrid_verts to list of Polygons
    ugrid_all_verts = dfmt.get_ugrid_verts(data_map_timesel)
    pol_shp_list = [Polygon(pol_data[~np.isnan(pol_data).all(axis=1)]) for pol_data in ugrid_all_verts]
    newdata = gpd.GeoDataFrame({'geometry': pol_shp_list},crs="EPSG:28992") #way more time efficient than doing it the loop
    
    for iV, varname in enumerate(varlist):
        varname_found = dfmt.get_varnamefromattrs(data_map_timesel,varname)
        
        data_sel_var = data_sel[varname_found]
        newdata[varname] = data_sel_var.to_numpy()
        
        fig, ax = plt.subplots()
        data_sel_var.ugrid.plot(cmap='viridis',edgecolor='face')
        fig.tight_layout()
    
    timestamp = data_map_timesel.time.dt.strftime('%Y%m%d').data
    file_shp = os.path.join(dir_shp,f'shp_{timestamp}')
    newdata.to_file(file_shp)

print(f'script runtime: {(dt.datetime.now()-dtstart_script).total_seconds():.2f} sec')
