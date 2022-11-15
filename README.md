[![generate-documentation](https://github.com/openearth/dfm_tools/actions/workflows/generate-documentation.yml/badge.svg)](https://github.com/openearth/dfm_tools/actions/workflows/generate-documentation.yml)

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
	- cdsapi/pydap: to download ERA5 and CMEMS data. Minimal pydap version is 3.3.0 (only available via conda-forge on 10-11-2022)
- launch Spyder:
	- open 'Spyder(dfm_tools_env)' via your windows start menu (not 'Spyder' or 'Spyder(Anaconda3)', since dfm_tools was installed in the dfm_tools_env environment only)
	- copy the code from [Example usage](#example-usage) to your own scripts to get started
	- Qt error upon launching Spyder?: remove the system/user environment variable 'qt_plugin_path' set by an old Delft3D4 installation procedure.
	- netCDF4 DLL error upon import in Spyder?: remove Anaconda paths from the Path user environment variable (https://github.com/spyder-ide/spyder/issues/19220)
- to update dfm_tools:
	- open anaconda prompt
	- ``conda activate dfm_tools_env``
	- ``python -m pip install --upgrade git+https://github.com/openearth/dfm_tools.git``


Getting started
--------
- [pdf with dfm_tools features and examples](https://nbviewer.org/github/openearth/dfm_tools/raw/pptx/docs/dfm_tools.pdf?flush_cache=true)
- [html documentation](https://htmlpreview.github.io/?https://github.com/openearth/dfm_tools/blob/master/docs/dfm_tools/index.html)
- [jupyter notebook with example code](https://github.com/openearth/dfm_tools/blob/master/notebooks/postprocessing_readme_example.ipynb)
- [example scripts](https://github.com/openearth/dfm_tools/tree/master/tests/examples)
- want to get updates about dfm_tools? Send an email to jelmer.veenstra@deltares.nl
