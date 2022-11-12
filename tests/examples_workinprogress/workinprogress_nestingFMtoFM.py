# -*- coding: utf-8 -*-
"""
Created on Sat Nov 12 20:11:29 2022

@author: veenstra
"""

import os
import xarray as xr
try: #0.3.1 release
    #from hydrolib.core.io.net.models import NetworkModel, Network
    from hydrolib.core.io.polyfile.models import PolyFile
except: #main branch and next release
    #from hydrolib.core.io.dflowfm.net.models import NetworkModel, Network
    from hydrolib.core.io.dflowfm.polyfile.models import PolyFile
from pathlib import Path
import pandas as pd
from dfm_tools.get_nc import get_netdata
from dfm_tools.hydrolib_helpers import pointlike_to_DataFrame
from dfm_tools.xarray_helpers import preprocess_hisnc
import numpy as np
import matplotlib.pyplot as plt
plt.close('all')
from scipy.spatial import KDTree
import contextily as ctx

#KDTree, finds minimal eucledian distance between points (maybe haversine would be better). Alternatively sklearn.neighbors.BallTree: tree = BallTree(data_lonlat_pd)
#TODO: add coordinate conversion of pli coordinates
#TODO: add max distance for nestpoints (eg sqrt of max cell size of large grid? How to determine to use 3/4/more points)
#TODO: add conversion of sigma/z-sigma model (time-varying z coords) to T3D with fixed z reference

dir_output = '.'


#NESTING PART 1
#get net and make KDTree of cellcenters
#file_net = r'p:\1230882-emodnet_hrsm\global_tide_surge_model\trunk\gtsm4.1\step11_global_1p25eu_net.nc' #TODO: cannot get netdata from this network, do in dfm_tools, hydrolib, xugrid, meshkernel?
file_net = r'p:\1230882-emodnet_hrsm\GTSMv5.0\runs\reference_GTSMv4.1_wiCA\output\gtsm_model_0000_map.nc'
crs_net = 'EPSG:4326'
data_net = get_netdata(file_net)
face_coords_pd = pd.DataFrame(dict(x=data_net.mesh2d_node_x,y=data_net.mesh2d_node_y))
tree_nest1 = KDTree(face_coords_pd)

#get and plot pli coordinates
file_pli_list = [Path(r'p:\i1000668-tcoms\03_newModel\01_input\02_bnd\pli\IndianOcean.pli'),
                 Path(r'p:\i1000668-tcoms\03_newModel\01_input\02_bnd\pli\PacificOcean.pli'),
                 Path(r'p:\i1000668-tcoms\03_newModel\01_input\02_bnd\pli\SeaOfJapan.pli'),
                 Path(r'p:\i1000668-tcoms\03_newModel\01_input\02_bnd\pli\TorresStrait.pli'),]

data_pol_pd_allfiles = pd.DataFrame()
for file_pli in file_pli_list:
    polyfile_object = PolyFile(file_pli)
    data_pol_pd_list = [pointlike_to_DataFrame(polyobj) for polyobj in polyfile_object.objects]
    data_pol_pd = pd.concat(data_pol_pd_list)
    data_pol_pd_allfiles = pd.concat([data_pol_pd_allfiles,data_pol_pd])

fig,ax = plt.subplots()
ax.plot(data_pol_pd_allfiles['x'],data_pol_pd_allfiles['y'],label='polyline')
ax.legend()
ctx.add_basemap(ax=ax,attribution=False,crs=crs_net)

#get and plot unique cell center coordinates
plicoords_distance1, plicoords_cellidx = tree_nest1.query(data_pol_pd_allfiles,k=4) #TODO: do on spherical globe (like gtsm obs snapping procedure)
cellidx_uniq = np.unique(plicoords_cellidx)
cellcoords = face_coords_pd.iloc[cellidx_uniq]
ax.plot(cellcoords['x'],cellcoords['y'],'xr',label='cellcenters')
cellcoords['name'] = 'nestpoint_' + cellcoords.index.astype(str)

#write nesting obspoints to file
file_obs = os.path.join(dir_output,'netpoint_obs.xyn')
cellcoords.to_csv(file_obs,sep='\t',index=False,header=False, float_format='%11.6f') #TODO: add hydrolib function once it exists


#NESTING PART 2
#get his
file_his = r'p:\1230882-emodnet_hrsm\GTSMv5.0\runs\reference_GTSMv4.1_wiCA\output\gtsm_model_0000_his.nc'
data_xr_his = xr.open_mfdataset(file_his,preprocess=preprocess_hisnc)

hisstations_pd = data_xr_his.stations.to_dataframe()
hisstations_pd_xy = hisstations_pd[['station_x_coordinate','station_y_coordinate']]
tree_nest2 = KDTree(hisstations_pd_xy)

plicoords_distance2, plicoords_nestpointidx = tree_nest2.query(data_pol_pd_allfiles, k=4)


#TODO FIX THIS PART
for iP, pli_row in data_pol_pd_allfiles.iterrows(): #TODO: maybe more efficient with xarray (multidim station retrieval)
    distance_oneset = plicoords_distance2[iP,:]
    weights_oneset = (1/distance_oneset) / (1/distance_oneset).sum()
    stations_oneset = hisstations_pd.index[plicoords_nestpointidx[iP,:]].tolist()
    data_xr_his_3points = data_xr_his.sel(stations=stations_oneset)
    data_xr_his_weightedmean = (data_xr_his_3points*weights_oneset).sum(dim='stations')
    
    #plot
    fig,ax = plt.subplots()
    data_xr_var_3points.plot.line(x='time')
    data_xr_var_3points.plot.line(x='time')






















