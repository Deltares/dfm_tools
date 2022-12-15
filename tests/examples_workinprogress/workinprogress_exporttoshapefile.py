# -*- coding: utf-8 -*-
"""
Created on Wed Oct 27 22:04:02 2021

@author: veenstra

"""

import os
import numpy as np
import matplotlib.pyplot as plt
plt.close('all')
import geopandas as gpd
from shapely.geometry import Polygon
import dfm_tools as dfmt
import datetime as dt

dtstart_script = dt.datetime.now()

dir_testinput = r'c:\DATA\dfm_tools_testdata'
dir_output = '.'
if not os.path.exists(dir_output):
    os.makedirs(dir_output)

#file_nc = os.path.join(r'p:\archivedprojects\11203850-coastserv\06-Model\waq_model\simulations\run0_20200319\DFM_OUTPUT_kzn_waq', 'kzn_waq_0*_map.nc')
file_nc = os.path.join(dir_testinput,'DFM_3D_z_Grevelingen','computations','run01','DFM_OUTPUT_Grevelingen-FM','Grevelingen-FM_0*_map.nc')
basename = os.path.basename(file_nc).replace('.','').replace('_0*_','_0000_')

if 'Grevelingen' in file_nc:
    crs = "EPSG:28992"
else:
    crs = "EPSG:4326"

varlist = ['mesh2d_sa1','Chlfa']#,'mesh2d_s1']

data_xr_map = dfmt.open_partitioned_dataset(file_nc)
vars_pd = dfmt.get_ncvarproperties(data_xr_map)

for timestep in [2,3]:#[0,10,20,30]:
    data_map_timesel = data_xr_map.isel(time=timestep)
    
    #data_sel = dfmt.get_Dataset_atdepths(data_xr=data_map_timesel, depths=0, reference='waterlevel') #top layer: 0m from waterlevel
    data_sel = dfmt.get_Dataset_atdepths(data_xr=data_map_timesel, depths=-4, reference='z0') #4m from model reference
    #data_sel = dfmt.get_Dataset_atdepths(data_xr=data_map_timesel, depths=0, reference='bedlevel') #bottomlayer: 0m from bedlevel
    
    #creating geodataframe with cells from ugrid_verts
    ugrid_all_verts = dfmt.get_ugrid_verts(data_map_timesel)
    pol_shp_list = [Polygon(verts_one[~np.isnan(verts_one).all(axis=1)]) for verts_one in ugrid_all_verts]
    newdata = gpd.GeoDataFrame({'geometry': pol_shp_list},crs=crs)
    
    for iV, varname in enumerate(varlist):
        try: #check if varname is present as key or in attributes
            varname_found = dfmt.get_varnamefromattrs(data_map_timesel,varname)
        except:
            print(f'varname {varname} not found in dataset')
            continue
        
        data_sel_var = data_sel[varname_found]
        newdata[varname] = data_sel_var.to_numpy() #can only have faces dimension (no time/layer)
        
        fig, ax = plt.subplots()
        data_sel_var.ugrid.plot(cmap='viridis',edgecolor='face')
        fig.tight_layout()
    
    timestamp = data_map_timesel.time.dt.strftime('%Y%m%d').data
    file_shp = os.path.join(dir_output,f'shp_{basename}_{timestamp}') #TODO: add depth+reference to filename
    newdata.to_file(file_shp)

print(f'script runtime: {(dt.datetime.now()-dtstart_script).total_seconds():.2f} sec')
