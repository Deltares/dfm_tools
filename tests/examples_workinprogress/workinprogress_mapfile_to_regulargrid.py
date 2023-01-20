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
import dfm_tools as dfmt

#TODO: also make it work for https://github.com/c-scale-community/use-case-hisea/blob/main/scripts/postprocessing/nc2regularGrid_listComprehension.py

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
    times_idx = None #one timestep in TTTZ model
    layers_idx = [45,46]
    varname_list = ['mesh2d_s1','mesh2d_sa1']#, 'mesh2d_ucx', 'mesh2d_ucy', 'mesh2d_tem1', 'mesh2d_sa1', 'mesh2d_water_quality_output_17', 'mesh2d_OXY']
dir_output = './data_output'
outname = os.path.basename(file_nc).replace('.nc','_regular.nc')

"""
file_nc : string, is the path to the 0'th partition of the DFlowFM map file output (use star for xugrid)
nx : integer, is the number of points in the x-direction (longitude) you wish to interpolate to. The points are evenly spaced.
ny : integer, is the number of points in the y-direction (latitude) you wish to interpolate to. The points are evenly spaced.
times_idx : numpy array or None or an array of times you want to do the interpolation for
layers_idx : numpy array, None or integer. The number of layers you wish to include. The script detect if there are layers or not. 
"""

if not os.path.exists(dir_output):
    os.makedirs(dir_output)
file_nc_out = os.path.join(dir_output, outname)

reg_x_vec = np.arange(lon_min,lon_max,lon_res)
reg_y_vec = np.arange(lat_min,lat_max,lat_res)

vars_pd = dfmt.get_ncvarproperties(file_nc=file_nc)

data_frommap_merged = dfmt.open_partitioned_dataset(file_nc.replace('_0000_','_0*_'))
if 'mesh2d_nLayers' in data_frommap_merged.dims: #TODO: make more generic
    dimn_layer = 'mesh2d_nLayers'
elif 'mesh2d_nLayers' in data_frommap_merged.dims:
    dimn_layer = 'mesh2d_nLayers'
else:
    dimn_layer = None

#select on time/layers
if times_idx is not None:
    data_frommap_merged = data_frommap_merged.isel(time=times_idx)
if layers_idx is not None and dimn_layer is not None:
    data_frommap_merged = data_frommap_merged.isel({dimn_layer:layers_idx}) #TODO: selecting 2/47 layers still gives 48 mesh2d_nInterfaces, is incorrect?
data_vars_list = list(data_frommap_merged.data_vars)

#setup Dataset with constant variables
data_xr_out = xr.Dataset()
data_xr_out['time'] = data_frommap_merged.time

lonvar_attrs = {'axis':'X',
                'reference':'geographical coordinates, WGS84 projection',
                'units':'degrees_east',
                '_CoordinateAxisType':'Lon',
                'long_name':'longitude',
                'valid_max':'180',
                'valid_min':'-180'}
lonvar = xr.DataArray(reg_x_vec, dims=('longitude'), attrs=lonvar_attrs)
data_xr_out['longitude'] = lonvar

latvar_attrs = {'axis':'Y',
                'reference':'geographical coordinates, WGS84 projection',
                'units':'degrees_north',
                '_CoordinateAxisType':'Lat',
                'long_name':'latitude',
                'valid_max':'90',
                'valid_min':'-90'}
latvar = xr.DataArray(reg_y_vec, dims=('latitude'), attrs=latvar_attrs)
data_xr_out['latitude'] = latvar

layervar_attrs = {'axis':'Z',
                  'long_name':'layer'}
layervar = xr.DataArray(layers_idx, dims=('layer'), attrs=layervar_attrs)
data_xr_out['layer'] = layervar

#data_frommap_z = data_frommap_merged['mesh2d_layer_z']
#if len(data_frommap_z.shape)>1: #TODO: re-enable this again
#    raise Exception('converting to regulargrid is currently only possible for sigmalayers')
depthvar = xr.DataArray(data_frommap_merged['mesh2d_layer_z'].to_numpy(), dims=('layer'), attrs=data_frommap_merged['mesh2d_layer_z'].attrs) #if not z layers, make exception
data_xr_out['depth'] = depthvar
data_xr_out = data_xr_out.set_coords('depth') #TODO: depth value klopt niet in top layers

#loop over other datavariables with time and faces dimensions
for varname in data_vars_list:
    if varname_list is not None: #optionally only use selected variables
        if varname not in varname_list:
            continue
    var_xr = data_frommap_merged[varname]
    if not 'time' in var_xr.dims: #process only variables with time varying values
        continue
    if not 'mesh2d_nFaces' in var_xr.dims: #process only variables with values on faces #TODO: make more generic
        continue
    print(varname)
    
    if dimn_layer in var_xr.dims:
        ntimes = var_xr.shape[var_xr.dims.index('time')]
        nlayers = var_xr.shape[var_xr.dims.index(dimn_layer)]
        field_array = np.empty((len(reg_y_vec),len(reg_x_vec),ntimes,nlayers))
        field_array[:] = np.nan
        for iT in range(ntimes):
            for iL in range(nlayers):
                var_xr_sel = var_xr.isel(time=iT).isel({dimn_layer:iL})
                field_array[:,:,iT,iL] = dfmt.scatter_to_regulargrid(xcoords=var_xr_sel.mesh2d_face_x, ycoords=var_xr_sel.mesh2d_face_y, reg_x_vec=reg_x_vec, reg_y_vec=reg_y_vec,
                                                                     values=var_xr_sel, method=method, maskland_dist=maskland_dist)[2]
        fieldvar_attrs = data_frommap_merged[varname].attrs
        fieldvar = xr.DataArray(field_array, dims=('latitude', 'longitude', 'time', 'layer'), attrs=fieldvar_attrs, name=varname)
        data_xr_out[varname] = fieldvar
    else:
        ntimes = var_xr.shape[var_xr.dims.index('time')]
        field_array = np.empty((len(reg_y_vec),len(reg_x_vec),ntimes))
        field_array[:] = np.nan
        for iT in range(ntimes):
            var_xr_sel = var_xr.isel(time=iT)
            field_array[:,:,iT] = dfmt.scatter_to_regulargrid(xcoords=var_xr_sel.mesh2d_face_x, ycoords=var_xr_sel.mesh2d_face_y, reg_x_vec=reg_x_vec, reg_y_vec=reg_y_vec,
                                                              values=var_xr_sel, method=method, maskland_dist=maskland_dist)[2]
        
        fieldvar_attrs = data_frommap_merged[varname].attrs
        fieldvar = xr.DataArray(field_array, dims=('latitude', 'longitude', 'time'), attrs=fieldvar_attrs, name=varname)
        data_xr_out[varname] = fieldvar

    print('done with variable %s' % varname)
    

time_elapsed = tm.time() - time_start
print('Duration: %f s' %time_elapsed)

for varname in varname_list:
    fig, ax = plt.subplots()#figsize=(12, 6))
    if 'layer' in data_xr_out[varname].dims:
        data_xr_out[varname].isel(time=0,layer=-1).plot(ax=ax,x='longitude',y='latitude',alpha=0.8,edgecolor='face')
    else:
        data_xr_out[varname].isel(time=0).plot(ax=ax,x='longitude',y='latitude',alpha=0.8,edgecolor='face')
    #source = ctx.providers.Esri.WorldImagery # ctx.providers.Stamen.Terrain (default), ctx.providers.CartoDB.Voyager, ctx.providers.NASAGIBS.ViirsEarthAtNight2012, ctx.providers.Stamen.Watercolor
    #ctx.add_basemap(ax, attribution=False, crs='epsg:4326', source=source)
    fig.savefig(os.path.join(dir_output, f'{model}_{varname}'))
data_xr_out.to_netcdf(file_nc_out)
#data_xr_out.close()


