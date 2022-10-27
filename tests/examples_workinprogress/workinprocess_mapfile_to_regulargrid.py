# -*- coding: utf-8 -*-
"""
Created on Wed Oct 26 13:59:19 2022

@author: veenstra
"""

#TODO: sturing resolutie en begin/eindcoordinaat
#TODO: sturing variabelenaam, tijdrange, layers

from dfm_tools.get_nc import get_ncmodeldata
from dfm_tools.get_nc_helpers import get_ncvarproperties #TODO BJORN: get_ncvarproperties instead of get_ncvardimlist
from dfm_tools.regulargrid import scatter_to_regulargrid

import os
import numpy as np
import time as tm
import xarray as xr
import matplotlib.pyplot as plt
plt.close('all')


time_start = tm.time()
"""change line below to the name of the 0th partition of your NetCDF map output"""

file_nc = r'p:\11200463-mixitin-mixotrophs\24.NS_D3DFM\Model_runs\NS_54\DFM_OUTPUT_DCSM-FM_4nm_waq\DCSM-FM_4nm_waq_0000_map.nc'
file_nc = r'p:\11202428-hisea\03-Model\Greece_model\waq_model\run2_DYNAMO_20200727\DFM_OUTPUT_tttz_waq\tttz_waq_0000_map.nc'
nx = 80
ny = 80

"""times and layers can also be specified by
"""
times_idx = 'all' #np.arange(40,44,2)
layers_idx = [45,46]
varname_list = ['mesh2d_s1', 'mesh2d_ucx', 'mesh2d_ucy', 'mesh2d_tem1', 'mesh2d_sa1', 'mesh2d_water_quality_output_17', 'mesh2d_OXY']

"""
file_nc : string, is the path to the 0'th partition of the DFlowFM map file output
xpts : integer, is the number of points in the x-direction (longitude) you wish to interpolate to. The points are evenly spaced.
ypts : integer, is the number of points in the y-direction (latitude) you wish to interpolate to. The points are evenly spaced.
tms : numpy array or 'all', an array of times you want to do the interpolation for
lrs : numpy array, 'all', or integer. The number of layers you wish to include. The script detect if there are layers or not. 
"""


dir_output = './data_output'
if not os.path.exists(dir_output):
    os.makedirs(dir_output)

vars_pd = get_ncvarproperties(file_nc=file_nc)

data_frommap_x = get_ncmodeldata(file_nc=file_nc, varname='mesh2d_face_x')
data_frommap_y = get_ncmodeldata(file_nc=file_nc, varname='mesh2d_face_y')
data_frommap_z = get_ncmodeldata(file_nc=file_nc, varname='mesh2d_layer_z', layer=layers_idx)
if len(data_frommap_z.shape)>1:
    raise Exception('converting to regulargrid is currently only possible for sigmalayers')
data_time = get_ncmodeldata(file_nc=file_nc, varname='time', timestep=times_idx)

outname = '%s_regular.nc' % os.path.split(file_nc)[1][0:-3]
file_nc_reg = os.path.join(dir_output, outname)

data_xr_source = xr.open_dataset(file_nc)

data_vars_list = list(data_xr_source.data_vars)

#setup Dataset with constant variables
data_xr_out = xr.Dataset()
x_grid, y_grid, dummy = scatter_to_regulargrid(xcoords=data_frommap_x, ycoords=data_frommap_y, ncellx=nx, ncelly=ny, values=data_frommap_x)
first_read = False
lon = x_grid[0, :]
lat = y_grid[:, 0]

lonvar_attrs = {'axis':'X', #TODO: get these from the original dataset after moving to xarray
                'reference':'geographical coordinates, WGS84 projection',
                'units':'degrees_east',
                '_CoordinateAxisType':'Lon',
                'long_name':'longitude',
                'valid_max':'180',
                'valid_min':'-180'}
lonvar = xr.DataArray(lon, dims=('lon'), attrs=lonvar_attrs)
data_xr_out['longitude'] = lonvar
data_xr_out = data_xr_out.set_coords('longitude')

latvar_attrs = {'axis':'Y', #TODO: get these from the original dataset after moving to xarray
                'reference':'geographical coordinates, WGS84 projection',
                'units':'degrees_east',
                '_CoordinateAxisType':'Lat',
                'long_name':'latitude',
                'valid_max':'90',
                'valid_min':'-90'}
latvar = xr.DataArray(lat, dims=('lat'), attrs=latvar_attrs)
data_xr_out['latitude'] = latvar
data_xr_out = data_xr_out.set_coords('latitude')

layervar_attrs = {'axis':'Z', #TODO: get these from the original dataset after moving to xarray
                  #'reference':'geographical coordinates, WGS84 projection',
                  'units':'m',
                  '_CoordinateZisPositive':'down',
                  '_CoordinateAxisType':'Height',
                  'long_name':'Depth'}
layervar = xr.DataArray(data_frommap_z, dims=('layer'), attrs=layervar_attrs)
data_xr_out['layer'] = layervar

timevar_attrs = {'units':'seconds since 2015-01-01 00:00:00',
                 'calendar':'standard',
                 'long_name':'time',
                 '_CoordinateAxisType':'Time'}
timevar = xr.DataArray(data_time, dims=('time'), attrs=timevar_attrs)
data_xr_out['time'] = timevar


skippedvars = []
for varname in data_vars_list:
    if varname_list is not None: #optionally only use selected variables
        if varname not in varname_list:
            continue
    var_xr = data_xr_source[varname]
    if not 'mesh2d_nFaces' in var_xr.dims: #skip variables on edges
        continue
    if 'mesh2d_nEdges' in var_xr.dims: #skip variables on edges
        continue
    
    print(varname)
    
    cell_size_x = (data_frommap_x.max() - data_frommap_x.min()) / nx
    cell_size_y = (data_frommap_x.max() - data_frommap_x.min()) / nx
    print('Cell size x: ' + str(cell_size_x))
    print('Cell size y: ' + str(cell_size_y))
    
    if ('time' in var_xr.dims) and ('mesh2d_nLayers' in var_xr.dims):
        data_frommap_var = get_ncmodeldata(file_nc=file_nc, varname=varname, timestep=times_idx, layer=layers_idx, return_xarray=True)
        ntimes = data_frommap_var.shape[var_xr.dims.index('time')]
        nlayers = data_frommap_var.shape[var_xr.dims.index('mesh2d_nLayers')]
        field_array = np.empty((nx,ny,ntimes,nlayers))
        for t in range(ntimes):
            for l in range(nlayers):
                field_array[:,:,t,l] = scatter_to_regulargrid(xcoords=data_frommap_x, ycoords=data_frommap_y, ncellx=nx, ncelly=ny,
                                                              values=data_frommap_var.isel(time=t,mesh2d_nLayers=l), method= 'nearest', maskland_dist=0.01)[2]
        
        fieldvar_attrs = data_xr_source[varname].attrs
        fieldvar = xr.DataArray(field_array, dims=('lon', 'lat', 'time', 'layer'), attrs=fieldvar_attrs) #, fill_value=-999. #TODO: why this fillvalue?
        data_xr_out[varname] = fieldvar
    elif 'time' in var_xr.dims:
        data_frommap_var = get_ncmodeldata(file_nc=file_nc, varname=varname, timestep=times_idx, return_xarray=True)
        ntimes = data_frommap_var.shape[var_xr.dims.index('time')]
        field_array = np.empty((nx,ny,ntimes))
        for t in range(ntimes):
            field_array[:,:,t] = scatter_to_regulargrid(xcoords=data_frommap_x, ycoords=data_frommap_y, ncellx=nx, ncelly=ny,
                                                        values=data_frommap_var.isel(time=t), method= 'nearest', maskland_dist=0.01)[2]
        
        fieldvar_attrs = data_xr_source[varname].attrs
        fieldvar = xr.DataArray(field_array, dims=('lon', 'lat', 'time'), attrs=fieldvar_attrs) #, fill_value=-999. #TODO: why this fillvalue?
        data_xr_out[varname] = fieldvar
    else:
        skippedvars.append(varname)
        continue
    
    print('done with variable %s' % varname)
    


data_xr_out.to_netcdf(file_nc_reg)
data_xr_out.close()

time_elapsed = tm.time() - time_start
print('Duration: %f s' %time_elapsed) #check how much time the script needs to run.



fig, ax = plt.subplots()#figsize=(12, 6))
data_xr = xr.open_dataset(r'c:\DATA\dfm_tools\tests\examples_workinprogress\data_output\tttz_waq_0000_map_regular.nc')
data_xr['mesh2d_sa1'].isel(time=0,layer=0).plot(ax=ax)
data_xr.close()



