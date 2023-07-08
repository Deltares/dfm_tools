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

dir_output = '.'
if not os.path.exists(dir_output):
    os.makedirs(dir_output)
export_kml = True

#file_nc = os.path.join(r'p:\archivedprojects\11203850-coastserv\06-Model\waq_model\simulations\run0_20200319\DFM_OUTPUT_kzn_waq', 'kzn_waq_0*_map.nc')
file_nc = dfmt.data.fm_grevelingen_map(return_filepath=True)
basename = os.path.basename(file_nc).replace('.','').replace('_0*_','_0000_')

if 'Grevelingen' in file_nc:
    crs = "EPSG:28992"
else:
    crs = "EPSG:4326"

varlist = ['mesh2d_sa1','mesh2d_Chlfa']#,'mesh2d_s1']

data_xr_map = dfmt.open_partitioned_dataset(file_nc)
data_xr_map = dfmt.rename_waqvars(data_xr_map)
vars_pd = dfmt.get_ncvarproperties(data_xr_map)

for iT, timestep in enumerate([2,3]):#[0,10,20,30]:
    data_map_timesel = data_xr_map.isel(time=timestep)
    
    #data_sel = dfmt.get_Dataset_atdepths(data_xr=data_map_timesel, depths=0, reference='waterlevel') #top layer: 0m from waterlevel
    data_sel = dfmt.get_Dataset_atdepths(data_xr=data_map_timesel, depths=-4, reference='z0') #4m from model reference
    #data_sel = dfmt.get_Dataset_atdepths(data_xr=data_map_timesel, depths=0, reference='bedlevel') #bottomlayer: 0m from bedlevel
    
    #creating geodataframe with cells from ugrid_verts
    ugrid_all_verts = data_map_timesel.grid.face_node_coordinates
    pol_shp_list = [Polygon(verts_one[~np.isnan(verts_one).all(axis=1)]) for verts_one in ugrid_all_verts]
    newdata = gpd.GeoDataFrame({'geometry': pol_shp_list},crs=crs)

    if iT==0 and export_kml: #export without variable
        #https://stackoverflow.com/questions/36222857/convert-geodataframe-polygons-to-kml-file
        import fiona
        fiona.supported_drivers['KML'] = 'rw'
        file_kml = os.path.join(dir_output,f'{basename}.kml') #TODO: add depth+reference to filename
        if os.path.exists(file_kml):#to avoid "DriverError: unsupported driver: 'LIBKML'"
            os.remove(file_kml)
        newdata.to_file(file_kml, driver='KML')
    
    for iV, varname in enumerate(varlist):
        if not hasattr(data_map_timesel,varname):
            print(f'varname {varname} not found in dataset')
            continue
        
        data_sel_var = data_sel[varname]
        newdata[varname] = data_sel_var.to_numpy() #can only have faces dimension (no time/layer)
        
        fig, ax = plt.subplots()
        data_sel_var.ugrid.plot(cmap='viridis')
        fig.tight_layout()
    
    timestamp = data_map_timesel.time.dt.strftime('%Y%m%d').data
    file_shp = os.path.join(dir_output,f'shp_{basename}_{timestamp}') #TODO: add depth+reference to filename
    newdata.to_file(file_shp) #TODO: solve "UserWarning: Column names longer than 10 characters will be truncated when saved to ESRI Shapefile."
    
print(f'script runtime: {(dt.datetime.now()-dtstart_script).total_seconds():.2f} sec')
