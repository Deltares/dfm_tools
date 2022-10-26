# -*- coding: utf-8 -*-
"""
Created on Wed Oct 26 13:59:19 2022

@author: veenstra
"""

#sturing: resolutie en begin/eindcoordinaat
#sturing: variabelenaam, tijdrange, layers

from dfm_tools.get_nc import get_netdata, get_ncmodeldata
from dfm_tools.get_nc_helpers import get_ncvarproperties #TODO BJORN: get_ncvarproperties instead of get_ncvardimlist
from dfm_tools.regulargrid import scatter_to_regulargrid

import os
import numpy as np
from netCDF4 import Dataset
import time as tm
import xarray as xr


"""original model located at
p:\11202428-hisea\03-Model\Greece_model\waq_model\
"""

time_start = tm.time()
"""change line below to where your DFlowFM output directory is located"""
# dir_input = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'input', 'DFM_OUTPUT_tttz_waq'))
dir_input = '/data/input'
#mapfile='tttz_waq_0000_map.nc'
"""change line below to the name of the 0th partition of your NetCDF map output"""

fp_in = r'p:\11202428-hisea\03-Model\Greece_model\waq_model\run2_DYNAMO_20200727\DFM_OUTPUT_tttz_waq\tttz_waq_0000_map.nc'
nx = 80
ny = 80

"""times and layers can also be specified by
tms = np.arange(40,44,2)
lrs = [0,10,20,30,46]
"""
treg = 'all'
#lrs = 'all'
#tms = np.arange(40,44,2)
lreg = [46]
"""
fileNC : string, is the path to the 0'th partition of the DFlowFM map file output
xpts : integer, is the number of points in the x-direction (longitude) you wish to interpolate to. The points are evenly spaced.
ypts : integer, is the number of points in the y-direction (latitude) you wish to interpolate to. The points are evenly spaced.
tms : numpy array or 'all', an array of times you want to do the interpolation for
lrs : numpy array, 'all', or integer. The number of layers you wish to include. The script detect if there are layers or not. 
"""





#dir_output = os.path.abspath(os.path.join(os.path.dirname(__file__),'data','output'))
#dir_output = os.getcwd() + '\data\output'
dir_output = './data_output'
if not os.path.exists(dir_output):
    os.makedirs(dir_output)
print(dir_output)
file_nc = fp_in
input_nc = Dataset(file_nc, 'r', format='NetCDF4')
time_old = input_nc.variables['time'][:]
if treg != 'all':
    time_old = np.take(time_old, treg)

vars_pd = get_ncvarproperties(file_nc=file_nc) #TODO BJORN: get_ncvarproperties instead of get_ncvardimlist


df = vars_pd

key_values = ['mesh2d_tem1','time', 'mesh2d_s1', 'mesh2d_ucx', 'mesh2d_ucy', 'mesh2d_tem1', 'mesh2d_sa1', 'mesh2d_water_quality_output_17', 'mesh2d_OXY', 'mesh2d_face_x', 'mesh2d_face_y', 'mesh2d_layer_z']

df = df.loc[df.index.isin(key_values)] #TODO BJORN: df.index instead of df['nc_varkeys']

"""
####################################################################################################################
#   Regularise all files with 3 dimensions (time, nFaces, layers). 
#   This will be equal to four dimensions in the regular grid format since nFaces is the x- and y- dimension.
####################################################################################################################
"""
df2 = df.loc[df['shape'].apply(len) == 3] #TODO BJORN: df['shape'].apply(len) instead of df['ndims']

data_frommap_x = get_ncmodeldata(file_nc=file_nc, varname='mesh2d_face_x')
data_frommap_y = get_ncmodeldata(file_nc=file_nc, varname='mesh2d_face_y')
data_frommap_z = get_ncmodeldata(file_nc=file_nc, varname='mesh2d_layer_z', layer='all')
time = get_ncmodeldata(file_nc=file_nc, varname='time', timestep=treg)
outname = '%s_regular.nc' % os.path.split(fp_in)[1][0:-3]
file_nc_reg = os.path.join(dir_output, outname)
dataset_xr = xr.Dataset()
#root_grp = Dataset(file_nc_reg, 'w', format='NETCDF4')
#root_grp.description = 'Example simulation data'
first_read = True
i = 0
for index, row in df2.iterrows():
    print(row.name) #TODO BJORN: row.name instead of row['nc_varkeys']
    if row['dimensions'][1] == 'mesh2d_nEdges':
        continue
    data_frommap_var = get_ncmodeldata(file_nc=file_nc, varname=row.name, timestep=treg, layer=lreg) #TODO BJORN: row.name instead of row['nc_varkeys']
    data_frommap_var = data_frommap_var.filled(np.nan)
    field_array = np.empty((data_frommap_var.shape[0], ny, nx, data_frommap_var.shape[-1]))
    tms = data_frommap_var.shape[0]
    lrs = data_frommap_var.shape[-1]
    trange = range(0, tms)
    lrange = range(0, lrs)
    
    cell_size_x = (data_frommap_x.max() - data_frommap_x.min()) / nx
    cell_size_y = (data_frommap_x.max() - data_frommap_x.min()) / nx
    print('Cell size x: ' + str(cell_size_x))
    print('Cell size y: ' + str(cell_size_y))
    
    A = np.empty((tms, lrs))
    A = np.array([scatter_to_regulargrid(xcoords=data_frommap_x, ycoords=data_frommap_y, ncellx=nx, ncelly=ny,
                                          values=data_frommap_var[t, :, l].flatten(), method= 'nearest', maskland_dist=0.01) for t in
                  trange for l in lrange])
  
    x_grid = A[0][0]
    y_grid = A[0][1]
    A = A[:, 2, :, :]
    A = np.moveaxis(A, [0], [2])
    subs = np.split(A, tms, axis=2)

    field_array[:, :, :, 0:lrs] = [subs[tn] for tn in trange]
    field_array = np.ma.masked_invalid(field_array)
    #field_array = field_array.filled(-999.)
    
    print('done with variable %s' % row.name) #TODO BJORN: row.name instead of row['nc_varkeys']

    if first_read:
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
        dataset_xr['longitude'] = lonvar
        dataset_xr = dataset_xr.set_coords('longitude')
        
        latvar_attrs = {'axis':'Y', #TODO: get these from the original dataset after moving to xarray
                        'reference':'geographical coordinates, WGS84 projection',
                        'units':'degrees_east',
                        '_CoordinateAxisType':'Lat',
                        'long_name':'latitude',
                        'valid_max':'90',
                        'valid_min':'-90'}
        latvar = xr.DataArray(lat, dims=('lat'), attrs=latvar_attrs)
        dataset_xr['latitude'] = latvar
        dataset_xr = dataset_xr.set_coords('latitude')
        
        layervar_attrs = {'axis':'Z', #TODO: get these from the original dataset after moving to xarray
                          #'reference':'geographical coordinates, WGS84 projection',
                          'units':'m',
                          '_CoordinateZisPositive':'down',
                          '_CoordinateAxisType':'Height',
                          'long_name':'Depth'}
        layervar = xr.DataArray(data_frommap_z[lreg], dims=('layer'), attrs=layervar_attrs)
        dataset_xr['layer'] = layervar
        
        timevar_attrs = {'units':'seconds since 2015-01-01 00:00:00',
                         'calendar':'standard',
                         'long_name':'time',
                         '_CoordinateAxisType':'Time'}
        timevar = xr.DataArray(time_old, dims=('time'), attrs=timevar_attrs)
        dataset_xr['time'] = timevar
        

    fieldName = row.name #TODO BJORN: row.name instead of row['nc_varkeys']
    fieldvar_attrs = {}
    for ncattr in input_nc.variables[fieldName].ncattrs():
        if ncattr == "_FillValue":
            continue
        fieldvar_attrs[ncattr] = input_nc.variables[fieldName].getncattr(ncattr)
    fieldvar = xr.DataArray(field_array, dims=('time', 'lat', 'lon', 'layer'), attrs=fieldvar_attrs) #, fill_value=-999. #TODO: why this fillvalue?
    dataset_xr[fieldName] = fieldvar
    

    
    first_read = False
    i += 1


"""
####################################################################################################################
#   Regularise all files with 2 dimensions (time, nFaces, layers).
#   This will be equal to 3 dimensions in the regular grid format since nFaces is the x- and y- dimension.
####################################################################################################################
"""
print('STARTING 2D')
df2 = df.loc[df['shape'].apply(len) == 2] #TODO BJORN: df['shape'].apply(len) instead of df['ndims']

excludeList = ['edge', 'face', 'x', 'y']
for index, row in df2.iterrows():
    test = any(n in str(row.name) for n in excludeList)#TODO BJORN: row.name instead of row['nc_varkeys']
    if not test:
        if row['dimensions'][1] == 'mesh2d_nEdges':
            continue
        ntimes = row['shape'][0]
        data_frommap_var = get_ncmodeldata(file_nc=file_nc, varname=row.name, timestep=treg)#TODO BJORN: row.name instead of row['nc_varkeys']
        data_frommap_var = data_frommap_var.filled(np.nan)
        field_array = np.empty((data_frommap_var.shape[0], ny, nx))
        trange = range(0, data_frommap_var.shape[0])
        tms = data_frommap_var.shape[0]
        A = np.array([scatter_to_regulargrid(xcoords=data_frommap_x, ycoords=data_frommap_y, ncellx=nx, ncelly=ny,
                                             values=data_frommap_var[t, :].flatten(), method='nearest', maskland_dist=0.01) for t in
                      trange])

        A = A[:, 2, :, :]
        field_array[:, :, :] = A
        field_array = np.ma.masked_invalid(field_array)
        #field_array = field_array.filled(-999.) #TODO: why this fillvalue?
        """write data to new netcdf"""
        fieldName = row.name#TODO BJORN: row.name instead of row['nc_varkeys']
        fieldName = row.name #TODO BJORN: row.name instead of row['nc_varkeys']
        fieldvar_attrs = {}
        for ncattr in input_nc.variables[fieldName].ncattrs():
            if ncattr == "_FillValue":
                continue
            fieldvar_attrs[ncattr] = input_nc.variables[fieldName].getncattr(ncattr)
        fieldvar = xr.DataArray(field_array, dims=('time', 'lat', 'lon'), attrs=fieldvar_attrs) #, fill_value=-999. #TODO: why this fillvalue?
        dataset_xr[fieldName] = fieldvar
        
dataset_xr.to_netcdf(file_nc_reg)
dataset_xr.close()

time_elapsed = tm.time() - time_start
print('Duration: %f s' %time_elapsed) #check how much time the script needs to run.




import xarray as xr
data_xr = xr.open_dataset(r'c:\DATA\dfm_tools\tests\examples_workinprogress\data_output\tttz_waq_0000_map_regular.nc')
data_xr['mesh2d_sa1'].isel(time=0,layer=0).plot()
data_xr.close()



