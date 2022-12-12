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

dir_testinput = r'c:\DATA\dfm_tools_testdata'
dir_output = '.'


varlist = ['mesh2d_sa1']#'Chlfa']#,'mesh2d_s1']
dir_shp = dir_output
if not os.path.exists(dir_shp):
    os.makedirs(dir_shp)
#file_nc = os.path.join(r'p:\archivedprojects\11203850-coastserv\06-Model\waq_model\simulations\run0_20200319\DFM_OUTPUT_kzn_waq', 'kzn_waq_0*_map.nc') #TODO: "ValueError: All NaN slice encountered" (is normally a runtimewarning, but not here)
file_nc = os.path.join(dir_testinput,'DFM_3D_z_Grevelingen','computations','run01','DFM_OUTPUT_Grevelingen-FM','Grevelingen-FM_0*_map.nc')
file_nc_nostar = file_nc.replace('0*','0000') #TODO: introduce support for * in dfm_tools definitions

data_xr_map = dfmt.open_partitioned_dataset(file_nc)
vars_pd = dfmt.get_ncvarproperties(file_nc=file_nc_nostar)

ugrid_all_verts = dfmt.get_ugrid_verts(data_xr_map)

pol_shp_list = [] #TODO: maybe there is a more direct way to convert xugrid to shapefile?
for pol_data in ugrid_all_verts: #[range(5000),:,:]
    pol_data_nonan = pol_data[~np.isnan(pol_data).all(axis=1)]
    pol_shp = Polygon(pol_data_nonan)
    pol_shp_list.append(pol_shp)

print('creating geodataframe with cells')
newdata = gpd.GeoDataFrame({'geometry': pol_shp_list},crs="EPSG:28992") #way more time efficient than doing it the loop

for iV, varname in enumerate(varlist):
    newdata[varname] = None

for timestep in [3]:#[0,10,20,30]:
    for iV, varname in enumerate(varlist):
        varname_found = dfmt.get_varnamefromattrs(file_nc_nostar,varname)
        data_map_timesel = data_xr_map.isel(time=timestep)
        
        #data_sel = dfmt.get_mapdata_atdepths(data_xr_map=data_map_timesel, depths=0, reference='waterlevel') #top layer: 0m from waterlevel
        #data_sel = dfmt.get_mapdata_atdepths(data_xr_map=data_map_timesel, depths=-5, reference='z0') #4m from model reference
        data_sel = dfmt.get_mapdata_atdepths(data_xr_map=data_map_timesel, depths=0, reference='bedlevel') #bottomlayer: 0m from bedlevel
        
        data_sel_var = data_sel[varname_found]
        
        newdata[varname] = data_sel_var.to_numpy()
        timestamp = data_sel_var.time.dt.strftime('%Y%m%d').data
        file_shp = os.path.join(dir_shp,f'shp_{varname}_{timestamp}')
        newdata.to_file(file_shp)
        
        fig, ax = plt.subplots()
        data_sel_var.ugrid.plot(cmap='viridis',edgecolor='face')
        fig.tight_layout()
        fig.savefig(os.path.join(file_shp,f'{varname}_{timestamp}'))

