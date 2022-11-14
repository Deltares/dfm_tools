dfm_tools
=========

dfm_tools are Python pre- and post-processing tools for D-FlowFM model files

Features
--------
- currently, we are updating dfm_tools to depend on [hydrolib-core](https://github.com/Deltares/HYDROLIB-core) for more robust and maintained core functionalities
- also, we are moving towards xarray for better handling of netCDF files, significant improvements in map file reading are therefore expected in the near future
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
	- [pdf with dfm_tools features and examples](https://nbviewer.org/github/openearth/dfm_tools/raw/master/docs/dfm_tools.pdf)
	- [html documentation](https://htmlpreview.github.io/?https://github.com/openearth/dfm_tools/blob/master/docs/dfm_tools/index.html)
	- [example scripts](https://github.com/openearth/dfm_tools/tree/master/tests/examples)
	- want to get updates about dfm_tools? Send an email to jelmer.veenstra@deltares.nl
	
Installation
--------
- download and install Anaconda 64 bit (with Python 3.8 or later) from https://www.anaconda.com/distribution/#download-section
- create Python environment and install dfm_tools:
	- open anaconda prompt
	- ``conda create --name dfm_tools_env -c conda-forge python=3.8 git spyder -y`` (you can also install a newer python version, in python 3.7 the dependency hydrolib seemed not installable)
	- ``conda activate dfm_tools_env``
	- ``conda install -c conda-forge shapely cartopy pyepsg geopandas contextily xarray dask netcdf4 bottleneck cdsapi pydap" -y`` (installs conda-forge requirements)
	- ``python -m pip install git+https://github.com/openearth/dfm_tools`` (this command installs dfm_tools and all required non-conda packages)
	- long paths error? Check last comment in https://github.com/Deltares/HYDROLIB-core/issues/327
	- to remove environment when necessary: ``conda remove -n dfm_tools_env --all``
- what are all these packages for?:
	- shapely for slicing 2D/3D data (conda-forge channel is necessary since main channel version is 1.6.4, minimal requirement is 1.7.0)
	- cartopy for satellite imagery, coastlines etc on plots (conda-forge channel recommended by cartopy developers, and currently also necessary for correct shapely version)
	- pyepsg is necessary for cartopy and probably also for other libraries
	- geopandas for shapefile related operations
	- contextily for satellite imagery on plots, seems faster than cartopy
	- xarray developers advise to install dependecies dask/netCDF4/bottleneck with conda-forge also: https://docs.xarray.dev/en/v0.8.0/installing.html
	- cdsapi/pydap: to download ERA5 and CMEMS data. Minimal pydap version is 3.3.0 (only available via conda-force on 10-11-2022)
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
More example scripts available at: https://github.com/openearth/dfm_tools/tree/master/tests/examples
```python
import os
import xarray as xr
import matplotlib.pyplot as plt
plt.close('all')
from dfm_tools.get_nc import get_netdata, get_ncmodeldata, plot_netmapdata
from dfm_tools.xarray_helpers import preprocess_hisnc

dir_testinput = os.path.join(r'n:\Deltabox\Bulletin\veenstra\info dfm_tools\test_input')
file_nc_map = os.path.join(dir_testinput,'DFM_sigma_curved_bend','DFM_OUTPUT_cb_3d','cb_3d_map.nc')
file_nc_his = os.path.join(dir_testinput,'DFM_sigma_curved_bend','DFM_OUTPUT_cb_3d','cb_3d_his.nc')

data_xr_his = xr.open_mfdataset(file_nc_his, preprocess=preprocess_hisnc)
stations_pd = data_xr_his['stations'].to_dataframe()

#retrieve his data and plot
fig, ax = plt.subplots(1,1,figsize=(10,5))
data_xr_his.waterlevel.plot.line(ax=ax, x='time')
ax.legend(data_xr_his.stations.to_series(),loc=1) #optional, to change legend location
fig.tight_layout()

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
ax.set_title(data_frommap_wl.var_varname)
ax.set_aspect('equal')

#plot salinity on map
data_frommap_sal = get_ncmodeldata(file_nc=file_nc_map, varname='mesh2d_sa1', timestep=2, layer=5)#, multipart=False)
fig, ax = plt.subplots()
pc = plot_netmapdata(ugrid_all.verts, values=data_frommap_sal[0,:,0], ax=None, linewidth=0.5, cmap="jet")
fig.colorbar(pc, ax=ax)
ax.set_title(data_frommap_sal.var_varname)
ax.set_aspect('equal')
```
