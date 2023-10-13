# -*- coding: utf-8 -*-
"""
Created on Sat Nov 12 20:11:29 2022

@author: veenstra
"""

import os
import xarray as xr
from pathlib import Path
import pandas as pd
from scipy.spatial import KDTree
import numpy as np
import matplotlib.pyplot as plt
plt.close('all')
import contextily as ctx
import dfm_tools as dfmt
import hydrolib.core.dflowfm as hcdfm


#TODO: add coordinate conversion of pli coordinates
#TODO: add max distance for nestpoints (eg sqrt of max cell size of large grid? How to determine to use 3/4/more points)
#TODO: add conversion of sigma/z-sigma model (time-varying z coords) to T3D with fixed z reference?

dir_output = '.'


#NESTING PART 1
#get net and make KDTree of cellcenters
#file_net = r'p:\1230882-emodnet_hrsm\global_tide_surge_model\trunk\gtsm4.1\step11_global_1p25eu_net.nc' #TODO: cannot get netdata from this network, do in dfm_tools, hydrolib, xugrid, meshkernel?
file_net = r'p:\1230882-emodnet_hrsm\GTSMv5.0\runs\reference_GTSMv4.1_wiCA\output\gtsm_model_0000_map.nc'
crs_net = 'EPSG:4326'
#data_xr_mapmerged = dfmt.open_partitioned_dataset(file_net)#.replace('_0000_','_0*_')) #TODO: ValueError: Dimensions {'nNetElem'} do not exist. Expected one or more of Frozen({}) (xugrid flexible gridnames support?)
data_xr = xr.open_dataset(file_net)
face_coords_pd = pd.DataFrame(dict(x=data_xr.FlowElemContour_x.mean(dim='nFlowElemContourPts'),y=data_xr.FlowElemContour_y.mean(dim='nFlowElemContourPts')))
tree_nest1 = KDTree(face_coords_pd)

#get and plot pli coordinates
file_pli_list = [Path(r'p:\i1000668-tcoms\03_newModel\01_input\02_bnd\pli\IndianOcean.pli'),
                 Path(r'p:\i1000668-tcoms\03_newModel\01_input\02_bnd\pli\PacificOcean.pli'),
                 Path(r'p:\i1000668-tcoms\03_newModel\01_input\02_bnd\pli\SeaOfJapan.pli'),
                 Path(r'p:\i1000668-tcoms\03_newModel\01_input\02_bnd\pli\TorresStrait.pli'),]

fig,ax = plt.subplots()
for file_pli in file_pli_list:
    polyfile_object = hcdfm.PolyFile(file_pli)
    data_pol_pd_list = [dfmt.pointlike_to_DataFrame(polyobj) for polyobj in polyfile_object.objects]
    data_pol_pd = pd.concat(data_pol_pd_list)
    
    #get and plot unique cell center coordinates
    plicoords_distance1, plicoords_cellidx = tree_nest1.query(data_pol_pd,k=4) #TODO: do on spherical globe (like gtsm obs snapping procedure)
    cellidx_uniq = np.unique(plicoords_cellidx)
    cellcoords = face_coords_pd.iloc[cellidx_uniq]
    cellcoords = cellcoords.reset_index() #revert from face numbers to 0-based, prevents SettingWithCopyWarning
    ax.plot(cellcoords['x'],cellcoords['y'],'x',label=f'{file_pli.name}_cellcenters')
    maxnumlen = cellcoords.index.astype(str).str.len().max()
    cellcoords['name'] = pd.Series(cellcoords.index).apply(lambda x: f'nestpoint_{x+1:0{maxnumlen}d}')
    
    #write nesting obspoints to file
    basename = file_pli.name.replace('.pli','')
    file_obs = os.path.join(dir_output,f'{basename}_obs.xyn')
    cellcoords.to_csv(file_obs,sep='\t',index=False,header=False, float_format='%11.6f') #TODO: add hydrolib function once it exists

    ax.plot(data_pol_pd['x'],data_pol_pd['y'],label=file_pli.name)

ax.legend()
ctx.add_basemap(ax=ax,attribution=False,crs=crs_net)


#NESTING PART 2
for file_pli in file_pli_list:
    kdtree_k = 4
    file_his = r'p:\1230882-emodnet_hrsm\GTSMv5.0\runs\reference_GTSMv4.1_wiCA\output\gtsm_model_0000_his.nc'
    data_xr_his = xr.open_mfdataset(file_his,preprocess=dfmt.preprocess_hisnc)
    data_xr_his_selvars = data_xr_his[['waterlevel']]#,'velocity_magnitude']]
    
    data_interp = dfmt.interp_hisnc_to_plipoints(data_xr_his=data_xr_his_selvars,file_pli=file_pli,kdtree_k=kdtree_k)
    fig,ax = plt.subplots()
    data_interp.waterlevel.drop_vars('plipoints').T.plot(ax=ax) #TODO: this does not work properly "TypeError: Dimensions of C (5473, 294) are incompatible with X (294) and/or Y (5474); see help(pcolormesh)" >> removing plipoints (names) solves the issue but not super desireable.
    
    # rename_dict = {'waterlevel':'waterlevelbnd',
    #                 'velocity_magnitude':'velocitybnd'}
    # data_interp = data_interp.rename(rename_dict)
    
    ForcingModel_object = dfmt.plipointsDataset_to_ForcingModel(plipointsDataset=data_interp)
    file_bc_out = file_pli.name.replace('.pli','.bc')
    
    print('saving bc file')
    ForcingModel_object.save(filepath=file_bc_out) #TODO REPORT: writing itself is fast, but takes quite a while to start writing (probably because of conversion)
    print('done')
















