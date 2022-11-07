# -*- coding: utf-8 -*-
"""
Created on Sun Oct 23 03:12:57 2022

@author: veenstra
"""

import os
import xarray as xr
import matplotlib.pyplot as plt
plt.close('all')
from scipy.spatial import KDTree
import pandas as pd
import numpy as np

dir_testinput = r'c:\DATA\dfm_tools_testdata'
dir_output = '.'

file_nc = os.path.join(dir_testinput,'DFM_3D_z_Grevelingen\\computations\\run01\\DFM_OUTPUT_Grevelingen-FM\\Grevelingen-FM_0000_his.nc')

station = ['GTSO-01','GTSO-02','GTSO-03','GTSO-04','GTSO-05','GTSO-06','GTSO-07',
           'GTSO-08','GTSO-09','GTSO-10','GTSO-11','GTSO-12','GTSO-13','GTSO-14',
           'GTSO-15','GTSO-16','GTSO-17','GTSO-18','GTSO-19','GTSO-20',
           'Bommenede','Grevelingen hevel West','Brouwerssluis binnen','Brouwerssluis binnen-hand']

data_xr = xr.open_dataset(file_nc)
data_xr_var = data_xr['waterlevel']


xyvals = np.array([[ 47703.73517405, 418471.4021889 ],
                   [ 51269.95437417, 419249.48637802],
                   [ 55160.37531975, 421389.21789809]])

#says lon/lat but is x/y (RD coordinates)
path_lonlat_pd = pd.DataFrame(xyvals,columns=['lon','lat'])

data_lon_flat = data_xr_var['station_x_coordinate'].to_numpy()
data_lat_flat = data_xr_var['station_y_coordinate'].to_numpy()
data_lonlat_pd = pd.DataFrame({'lon':data_lon_flat,'lat':data_lat_flat})

#KDTree, finds minimal eucledian distance between points (maybe haversine would be better)
tree = KDTree(data_lonlat_pd) #alternatively sklearn.neighbors.BallTree: tree = BallTree(data_lonlat_pd)
distance, data_lonlat_idx = tree.query(path_lonlat_pd, k=3) #TODO: maybe add outofbounds treshold for distance
#data_lonlat_pd.iloc[data_lonlat_idx]

for iS, stat_row in path_lonlat_pd.iterrows():
    distance_oneset = distance[iS,:]
    weights_oneset = (1/distance_oneset) / (1/distance_oneset).sum()
    data_xr_var_3points = data_xr_var.isel(stations=data_lonlat_idx[iS,:])
    data_xr_var_weightedmean = (data_xr_var_3points*weights_oneset).sum(dim='stations')
    
    #plot
    fig,ax = plt.subplots()
    data_xr_var_3points.plot.line(x='time')
    data_xr_var_3points.plot.line(x='time')


