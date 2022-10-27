# -*- coding: utf-8 -*-
"""
Created on Wed Oct 26 13:59:19 2022

@author: veenstra
"""

import os
import numpy as np
import time as tm
import xarray as xr
import matplotlib.pyplot as plt
plt.close('all')
#import contextily as ctx

from dfm_tools.get_nc import get_ncmodeldata
from dfm_tools.get_nc_helpers import get_ncvarproperties#, get_varnamefrom_keyslongstandardname
from dfm_tools.regulargrid import scatter_to_regulargrid

model = 'DCSM'
method = 'linear' #for scipy.interpolate.griddata: https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.griddata.html

time_start = tm.time()
if model=='DCSM':
    file_nc = r'p:\11200463-mixitin-mixotrophs\24.NS_D3DFM\Model_runs\NS_54\DFM_OUTPUT_DCSM-FM_4nm_waq\DCSM-FM_4nm_waq_0000_map.nc'
    multipart = True
    lat_min, lat_max, lat_res = 51.0, 57.0, 0.04
    lon_min, lon_max, lon_res = -2.0,  6.0, 0.04
    maskland_dist = 0.1 #smaller value gives all nans
    times_idx = [0,2]#'all' #np.arange(40,44,2)
    layers_idx = [49]
    varname_list = ['mesh2d_greenC_1']#,'mesh2d_diatC_1','mesh2d_zooC_1','mesh2d_cmC_1','mesh2d_diatChl_1','mesh2d_greenChl_1','mesh2d_cmChl_1']
elif model=='TTTZ':
    file_nc = r'p:\11202428-hisea\03-Model\Greece_model\waq_model\run2_DYNAMO_20200727\DFM_OUTPUT_tttz_waq\tttz_waq_0000_map.nc'
    multipart = True
    lat_min, lat_max, lat_res = 36.6, 38.0, 0.02
    lon_min, lon_max, lon_res = 22.8, 24.3, 0.02
    maskland_dist = 0.01
    times_idx = 'all' #one timestep in TTTZ model
    layers_idx = [45,46]
    varname_list = ['mesh2d_s1','mesh2d_sa1']#, 'mesh2d_ucx', 'mesh2d_ucy', 'mesh2d_tem1', 'mesh2d_sa1', 'mesh2d_water_quality_output_17', 'mesh2d_OXY']
dir_output = './data_output'
outname = os.path.basename(file_nc).replace('.nc','_regular.nc')

"""
file_nc : string, is the path to the 0'th partition of the DFlowFM map file output
nx : integer, is the number of points in the x-direction (longitude) you wish to interpolate to. The points are evenly spaced.
ny : integer, is the number of points in the y-direction (latitude) you wish to interpolate to. The points are evenly spaced.
times_idx : numpy array or 'all', an array of times you want to do the interpolation for
layers_idx : numpy array, 'all', or integer. The number of layers you wish to include. The script detect if there are layers or not. 
"""

if not os.path.exists(dir_output):
    os.makedirs(dir_output)
file_nc_out = os.path.join(dir_output, outname)

reg_x_vec = np.arange(lon_min,lon_max,lon_res)
reg_y_vec = np.arange(lat_min,lat_max,lat_res)

vars_pd = get_ncvarproperties(file_nc=file_nc)

data_frommap_x = get_ncmodeldata(file_nc=file_nc, varname='mesh2d_face_x', multipart=multipart)
data_frommap_y = get_ncmodeldata(file_nc=file_nc, varname='mesh2d_face_y', multipart=multipart)
#get_varnamefrom_keyslongstandardname(file_nc,'meshd2_layer_z')
data_frommap_z = get_ncmodeldata(file_nc=file_nc, varname='mesh2d_layer_z', layer=layers_idx, multipart=multipart) #TODO: depth is now not really used
if len(data_frommap_z.shape)>1:
    raise Exception('converting to regulargrid is currently only possible for sigmalayers')
#data_frommap_sigma = get_ncmodeldata(file_nc=file_nc, varname='mesh2d_layer_sigma', layer=layers_idx, multipart=multipart)
#data_frommap_sigma_z = get_ncmodeldata(file_nc=file_nc, varname='mesh2d_layer_sigma_z', layer=layers_idx, multipart=multipart)

data_xr_source = xr.open_dataset(file_nc)
data_vars_list = list(data_xr_source.data_vars)

#setup Dataset with constant variables
data_xr_out = xr.Dataset()

timevar = data_xr_source.time
if times_idx != 'all':
    timevar = timevar.isel(time=times_idx)
data_xr_out['time'] = timevar

lonvar_attrs = {'axis':'X',
                'reference':'geographical coordinates, WGS84 projection',
                'units':'degrees_east',
                '_CoordinateAxisType':'Lon',
                'long_name':'longitude',
                'valid_max':'180',
                'valid_min':'-180'}
lonvar = xr.DataArray(reg_x_vec, dims=('lon'), attrs=lonvar_attrs)
data_xr_out['longitude'] = lonvar
data_xr_out = data_xr_out.set_coords('longitude')

latvar_attrs = {'axis':'Y',
                'reference':'geographical coordinates, WGS84 projection',
                'units':'degrees_north',
                '_CoordinateAxisType':'Lat',
                'long_name':'latitude',
                'valid_max':'90',
                'valid_min':'-90'}
latvar = xr.DataArray(reg_y_vec, dims=('lat'), attrs=latvar_attrs)
data_xr_out['latitude'] = latvar
data_xr_out = data_xr_out.set_coords('latitude')

layervar_attrs = {'axis':'Z',
                  'long_name':'layer'}
layervar = xr.DataArray(layers_idx, dims=('layer'), attrs=layervar_attrs)
data_xr_out['layerno'] = layervar

depthvar_attrs = {'axis':'Z', #TODO: get these from the original dataset after moving to xarray
                  'units':'m',
                  #'_CoordinateZisPositive':'down',
                  #'_CoordinateAxisType':'Height',
                  'long_name':'Depth'}
depthvar = xr.DataArray(data_frommap_z, dims=('layer'), attrs=depthvar_attrs)
data_xr_out['mesh2d_layer_z'] = depthvar


#loop over other datavariables with tima and faces dimensions
for varname in data_vars_list:
    if varname_list is not None: #optionally only use selected variables
        if varname not in varname_list:
            continue
    var_xr = data_xr_source[varname]
    if not 'time' in var_xr.dims: #process only variables with time varying values
        continue
    if not 'mesh2d_nFaces' in var_xr.dims: #process only variables with values on faces
        continue
    
    print(varname)
    
    if 'mesh2d_nLayers' in var_xr.dims:
        data_frommap_var = get_ncmodeldata(file_nc=file_nc, varname=varname, timestep=times_idx, layer=layers_idx, return_xarray=True, multipart=multipart)
        ntimes = data_frommap_var.shape[var_xr.dims.index('time')]
        nlayers = data_frommap_var.shape[var_xr.dims.index('mesh2d_nLayers')]
        field_array = np.empty((len(reg_y_vec),len(reg_x_vec),ntimes,nlayers))
        for t in range(ntimes):
            for l in range(nlayers):
                field_array[:,:,t,l] = scatter_to_regulargrid(xcoords=data_frommap_x, ycoords=data_frommap_y, reg_x_vec=reg_x_vec, reg_y_vec=reg_y_vec,
                                                              values=data_frommap_var.isel(time=t,mesh2d_nLayers=l), method=method, maskland_dist=maskland_dist)[2]
        fieldvar_attrs = data_xr_source[varname].attrs
        fieldvar = xr.DataArray(field_array, dims=('lat', 'lon', 'time', 'layer'), attrs=fieldvar_attrs, name=varname)
        data_xr_out[varname] = fieldvar
    else:
        data_frommap_var = get_ncmodeldata(file_nc=file_nc, varname=varname, timestep=times_idx, return_xarray=True, multipart=multipart)
        ntimes = data_frommap_var.shape[var_xr.dims.index('time')]
        field_array = np.empty((len(reg_y_vec),len(reg_x_vec),ntimes))
        for t in range(ntimes):
            field_array[:,:,t] = scatter_to_regulargrid(xcoords=data_frommap_x, ycoords=data_frommap_y, reg_x_vec=reg_x_vec, reg_y_vec=reg_y_vec,
                                                        values=data_frommap_var.isel(time=t), method=method, maskland_dist=maskland_dist)[2]
        
        fieldvar_attrs = data_xr_source[varname].attrs
        fieldvar = xr.DataArray(field_array, dims=('lat', 'lon', 'time'), attrs=fieldvar_attrs, name=varname)
        data_xr_out[varname] = fieldvar

    print('done with variable %s' % varname)
    

time_elapsed = tm.time() - time_start
print('Duration: %f s' %time_elapsed)

for varname in varname_list:
    fig, ax = plt.subplots()#figsize=(12, 6))
    if 'layer' in data_xr_out[varname].dims:
        data_xr_out[varname].isel(time=0,layer=-1).plot(ax=ax,x='longitude',y='latitude',alpha=0.8)
    else:
        data_xr_out[varname].isel(time=0).plot(ax=ax,x='longitude',y='latitude',alpha=0.8)
    #source = ctx.providers.Esri.WorldImagery # ctx.providers.Stamen.Terrain (default), ctx.providers.CartoDB.Voyager, ctx.providers.NASAGIBS.ViirsEarthAtNight2012, ctx.providers.Stamen.Watercolor
    #ctx.add_basemap(ax, attribution=False, crs='epsg:4326', source=source)
        
#data_xr_out.to_netcdf(file_nc_out)
#data_xr_out.close()


