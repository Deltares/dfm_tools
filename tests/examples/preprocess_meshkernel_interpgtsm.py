# -*- coding: utf-8 -*-
"""
Created on Wed Apr 26 12:25:27 2023

@author: veenstra
"""

# interpolation of bathymetry to gtsm grid.
# before may 2023, this was done with uds.ugrid.sel(xslice) which subsets grid on face_coordinates. However, "around the back" cells were selected in for instance slice(-60,-65). This caused the nodes in this case to range from -180 to 150, which resulted in active interpolation of almost the entire gebco grid/dataset.
# now we are selecting the node_x coordinates between slices, the resulting list of xy coordinates is used for the gebco interpolation, effectively requiring only a lon-slice of the gebco dataset at the time.

import xarray as xr
import xugrid as xu
import datetime as dt

ds_gebco = xr.open_dataset('p:\\metocean-data\\open\\GEBCO\\2022\\GEBCO_2022.nc')

lat_max = -174 #TODO: 180

#extend gebco to lon=-180 and lat=90 to make interpolation on all left nodes and the top node possible
ds_gebco = ds_gebco.reset_index(['lat','lon'])
ds_gebco['lon'].values[0] = -180 #replace -179.99791667 with -180
ds_gebco['lat'].values[-1] = 90 #replace 89.99791667 with 90
ds_gebco = ds_gebco.set_index({'lat':'lat','lon':'lon'})

file_net = r'p:\1230882-emodnet_hrsm\global_tide_surge_model\trunk\gtsm4.1\step11_global_1p25eu_net.nc'
uds = xu.open_dataset(file_net)
nnodes = uds.dims[uds.grid.node_dimension]

stepsize = 5 #20 seems optimal, 5 takes 900-1000 sec, 10 takes 500-600 sec, 20 takes 95 sec, 30 takes 377 sec
print(f'interpolating GEBCO to {nnodes} nodes in 360/{stepsize}={360/stepsize} steps:')
dtstart = dt.datetime.now()
for i in range(-180, lat_max, stepsize):
    xslice = slice(i,i+stepsize)
    
    bool_nodeinslice = (uds.grid.node_coordinates[:,0] >= xslice.start) & (uds.grid.node_coordinates[:,0] < xslice.stop)
    print(xslice, ':', bool_nodeinslice.sum(), 'nodes')
    
    x_sel, y_sel = uds.grid.node_coordinates[bool_nodeinslice].T
    x_sel_ds = xr.DataArray(x_sel,dims=(uds.grid.node_dimension))
    y_sel_ds = xr.DataArray(y_sel,dims=(uds.grid.node_dimension))
    z_sel = ds_gebco.interp(lon=x_sel_ds, lat=y_sel_ds)
    uds['NetNode_z'][bool_nodeinslice] = z_sel.elevation

print('plot data grid')
print(uds['NetNode_z'].isnull().sum()) #check if there are missing z-values
#uds.NetNode_z.ugrid.plot(center=False)

print(f'{(dt.datetime.now()-dtstart).total_seconds():.2f} sec')





