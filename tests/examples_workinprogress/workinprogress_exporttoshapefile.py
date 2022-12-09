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


varlist = ['Chlfa']#,'mesh2d_s1']
dir_shp = dir_output
if not os.path.exists(dir_shp):
    os.makedirs(dir_shp)
file_nc = os.path.join(r'p:\archivedprojects\11203850-coastserv\06-Model\waq_model\simulations\run0_20200319\DFM_OUTPUT_kzn_waq', 'kzn_waq_0*_map.nc')
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
newdata = gpd.GeoDataFrame({'geometry': pol_shp_list},crs="EPSG:4326") #way more time efficient than doing it the loop

for iV, varname in enumerate(varlist):
    newdata[varname] = None

for timestep in [6]:#[0,10,20,30]:
    for iV, varname in enumerate(varlist):
        varname_found = dfmt.get_varnamefromattrs(file_nc_nostar,varname)
        data_var = data_xr_map[varname_found].isel(time=timestep)
        if 'nmesh2d_layer' in data_var.dims: #TODO: dimname could also be layno
            data_var = data_var.ffill(dim='nmesh2d_layer').bfill(dim='nmesh2d_layer') #replace nans in toplayer with first non-nan value in layers below (and for bottomlayer with first non-nan value in layers above)
            data_var_sel = data_var.isel(nmesh2d_layer=-1) #since nans are now filled, this provides the toplayer of the model
        else:
            data_var_sel = data_var
        
        newdata[varname] = data_var_sel.to_numpy()
        timestamp = data_var_sel.time.dt.strftime('%Y%m%d').data
        file_shp = os.path.join(dir_shp,f'shp_{varname}_{timestamp}')
        newdata.to_file(file_shp)
        
        fig, ax = plt.subplots(figsize=(6,7))
        data_var_sel.ugrid.plot(cmap='viridis')
        fig.tight_layout()
        fig.savefig(os.path.join(file_shp,f'{varname}_{timestamp}'))

