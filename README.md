dfm_tools
=========

dfm_tools are Python pre- and post-processing tools for D-FlowFM model files

Features
--------
- currently, we are updating dfm_tools to depend on [hydrolib-core](https://github.com/Deltares/HYDROLIB-core) for more robust and maintained core functionalities
- also, we are moving towards xarray for better handling of netCDF files
- post-processing and plotting:
	- support for net/map/fou/rst partitions (merge partitions and delete ghostcells) (and many other netCDF files)
	- support for flexible meshes (containing triangles, squares, pentagons, hexagons or any other shape)
	- select all data based on variable, timestep/datetime, layer, station (not yet on depth)
	- retrieve lists of variables, timesteps/datetimes, stations, cross sections, general structures
	- data selection/plotting by polyline/crossection (slicing the ugrid mapdata)
	- add satellite background to your plots
	- plotting z,t-plots
- pre-processing:
	- read and write almost all FM input data with [hydrolib-core](https://github.com/Deltares/HYDROLIB-core)
	- converting and plotting this data with helper functions in dfm_tools
	- e.g.: interpolating CMEMS data to model boundary and writing 2D/3D boundary condition files
- to get started:
	- html documentation: https://htmlpreview.github.io/?https://github.com/openearth/dfm_tools/blob/master/docs/dfm_tools/index.html
	- example scripts: https://github.com/openearth/dfm_tools/tree/master/tests/examples
	- examples of (mostly unformatted) figures created by this pytest testbank: n:\\Deltabox\\Bulletin\\veenstra\\info dfm_tools
	- want to get updates about dfm_tools? Send an email to jelmer.veenstra@deltares.nl
	
Installation
--------
- download and install Anaconda 64 bit (with Python 3.8 or later) from https://www.anaconda.com/distribution/#download-section
- install dfm_tools from github:
	- open anaconda prompt
	- ``conda create --name dfm_tools_env -c conda-forge python=3.8 git spyder -y`` (you can also install a newer python version, in python 3.7 the dependency hydrolib seemed not installable)
	- ``conda activate dfm_tools_env``
	- ``conda install -c conda-forge shapely cartopy pyepsg geopandas contextily xarray dask netcdf4 bottleneck cdsapi motuclient -y`` (installs conda-forge requirements)
	- ``python -m pip install git+https://github.com/openearth/dfm_tools`` (this command installs dfm_tools and all required non-conda packages)
	- ``python -m pip uninstall hydrolib-core -y`` (temporarily necessary to install the correct version of hydrolib-core)
	- ``python -m pip install git+https://github.com/deltares/hydrolib-core`` (temporarily necessary to install the correct version of hydrolib-core)
	- to remove venv when necessary: ``conda remove -n dfm_tools_env --all``
- what are all these packages for?:
	- shapely for slicing 2D/3D data (conda-forge channel is necessary since main channel version is 1.6.4, minimal requirement is 1.7.0)
	- cartopy for satellite imagery, coastlines etc on plots (conda-forge channel recommended by cartopy developers, and currently also necessary for correct shapely version)
	- pyepsg is necessary for cartopy and probably also for other libraries
	- geopandas for shapefile related operations
	- contextily for satellite imagery on plots, seems faster than cartopy
	- xarray developers advise to install dependecies dask/netCDF4/bottleneck with conda-forge also: https://docs.xarray.dev/en/v0.8.0/installing.html
	- cdsapi/motuclient: to download ERA5 and CMEMS data
- launch Spyder:
	- open 'Spyder(dfm_tools_env)' via your windows start menu (not 'Spyder' or 'Spyder(Anaconda3)', since dfm_tools was installed in the dfm_tools_env environment only)
	- copy the code from [Example usage](#example-usage) to your own scripts to get started
	- Qt error upon launching Spyder?: remove the system/user environment variable 'qt_plugin_path' set by an old Delft3D4 installation procedure.
	- netCDF4 DLL error upon import in Spyder?: remove Anaconda paths from the Path user environment variable (https://github.com/spyder-ide/spyder/issues/19220)
- to update dfm_tools:
	- open anaconda prompt
	- ``conda activate dfm_tools_env``
	- ``python -m pip install --upgrade git+https://github.com/openearth/dfm_tools.git``


Example usage
--------

```python
#data retrieval is easy, just use get_ncmodeldata() with file_nc argument
#then use the feedback in the error messages to set other arguments like varname, timestep, station and layer
from dfm_tools.get_nc import get_ncmodeldata
data_fromnc = get_ncmodeldata(file_nc='yourfile.nc')
```

```python
#this example includes plotting and using the metadata of the retrieved data
#import statements
import os
import xarray as xr
import matplotlib.pyplot as plt
plt.close('all')
from dfm_tools.get_nc import get_netdata, get_ncmodeldata, plot_netmapdata
from dfm_tools.get_nc_helpers import get_ncvardimlist, get_timesfromnc, get_hisstationlist
dir_testinput = r'c:\DATA\dfm_tools_testdata'

#uncomment the line below, copy data locally and change this path to increase performance
#dir_testinput = os.path.join(r'n:\Deltabox\Bulletin\veenstra\info dfm_tools\test_input')
file_nc_map = os.path.join(dir_testinput,'DFM_sigma_curved_bend','DFM_OUTPUT_cb_3d','cb_3d_map.nc')
file_nc_his = os.path.join(dir_testinput,'DFM_sigma_curved_bend','DFM_OUTPUT_cb_3d','cb_3d_his.nc')
data_xr_his = xr.open_dataset(file_nc_his)
stations_pd = data_xr_his.station_name.astype(str).to_pandas()

#get lists with vars/dims, times, station/crs/structures
vars_pd, dims_pd = get_ncvardimlist(file_nc=file_nc_map)
times_pd = get_timesfromnc(file_nc=file_nc_map)
statlist_pd = get_hisstationlist(file_nc=file_nc_his, varname='station_name')

#retrieve his data
#data_fromhis_wl = get_ncmodeldata(file_nc=file_nc_his, varname='waterlevel', station='all', timestep= 'all')
fig, ax = plt.subplots(1,1,figsize=(10,5))
for iS, station in enumerate(stations_pd):
    data_fromhis_wl = data_xr_his.waterlevel.isel(stations=iS)
    ax.plot(data_fromhis_wl.time,data_fromhis_wl,'-', label=station)
ax.legend()
ax.set_ylabel('%s (%s)'%(data_fromhis_wl.attrs['long_name'], data_fromhis_wl.attrs['units']))

#plot net/grid
ugrid_all = get_netdata(file_nc=file_nc_map)#,multipart=False)
fig, ax = plt.subplots()
pc = plot_netmapdata(ugrid_all.verts, values=None, ax=None, linewidth=0.5, color="crimson", facecolor="None")
ax.set_aspect('equal')

#plot water level on map
data_frommap_wl = get_ncmodeldata(file_nc=file_nc_map, varname='mesh2d_s1', timestep=3)#, multipart=False)
fig, ax = plt.subplots()
pc = plot_netmapdata(ugrid_all.verts, values=data_frommap_wl[0,:], ax=None, linewidth=0.5, cmap="jet")
pc.set_clim([-0.5,1])
fig.colorbar(pc, ax=ax)
ax.set_title('%s (%s)'%(data_frommap_wl.var_varname, data_frommap_wl.var_ncattrs['units']))
ax.set_aspect('equal')

#plot salinity on map
data_frommap_sal = get_ncmodeldata(file_nc=file_nc_map, varname='mesh2d_sa1', timestep=2, layer=5)#, multipart=False)
fig, ax = plt.subplots()
pc = plot_netmapdata(ugrid_all.verts, values=data_frommap_sal[0,:,0], ax=None, linewidth=0.5, cmap="jet")
fig.colorbar(pc, ax=ax)
ax.set_title('%s (%s)'%(data_frommap_sal.var_varname, data_frommap_sal.var_ncattrs['units']))
ax.set_aspect('equal')

#print contents of retrieved data withing data_frommap_sal variable
print_var = data_frommap_sal
print('++++++\nthe data in the variable %s is:\n%s\n'%(print_var.var_varname, print_var))
print('++++++\nthe time indices and times in the variable %s are:\n%s\n'%(print_var.var_varname, print_var.var_times))
#print('++++++\nthe station indices and station names in the variable %s are:\n%s\n'%(print_var.var_varname, print_var.var_stations))
print('++++++\nthe layer indices in the variable %s are:\n%s\n'%(print_var.var_varname, print_var.var_layers))
print('++++++\nthe shape of the variable %s is:\n%s\n'%(print_var.var_varname, print_var.shape))
print('++++++\nthe dimensions of the variable %s are (copied from netCDF variable):\n%s\n'%(print_var.var_varname, print_var.var_dimensions))
print('++++++\nthe netCDF variable where the data in variable %s comes from is:\n%s\n'%(print_var.var_varname, print_var.var_ncvarobject))
print('++++++\nsome example contents of this netCDF variable:')
print('\tthe dimension names of the netCDF variable %s are:\n\t\t%s'%(print_var.var_varname, print_var.var_dimensions))
print('\tthe shape of the netCDF variable %s is:\n\t\t%s'%(print_var.var_varname, print_var.var_shape))
print('\tthe units of the netCDF variable %s are:\n\t\t%s'%(print_var.var_varname, print_var.var_ncattrs['units']))
print('\tthe long_name of the netCDF variable %s is:\n\t\t%s'%(print_var.var_varname, print_var.var_ncattrs['long_name']))
print('\tthe standard_name of the netCDF variable %s is:\n\t\t%s'%(print_var.var_varname, print_var.var_ncattrs['standard_name']))
```
