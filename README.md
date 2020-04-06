dfm_tools
=========

dfm_tools are Python post-processing tools for Delft3D FM model outputfiles and other netCDF files

- Free software: GNU General Public License v3


Features
--------
- supported formats:
	- D-Flow FM output data (net, map, his, fou, rst files)
	- almost any other netCDF (ERA5, hirlam, SFINCS map, SFINCS his, Sobek observation)
	- Delft3D netCDF output files (you can get netcdf output with keywords in your mdf)
	- converted Delft3D and waqua data (converted to netCDF with getdata.pl) (Delft3D conversion with getdata.pl is not flawless, preferably rerun with netCDF as outputformat instead)
- data handling:
	- select all data based on variable, timestep/datetime, layer, station (not yet on depth)
	- get feedback about available variables, timesteps/datetimes, layers, stations when you retrieve the wrong ones
	- retrieve lists of variables, timesteps/datetimes, stations, cross sections, general structures
	- selection/plotting by polyline/crossection (slicing the ugrid data)
	- merge partitions and delete ghostcells automatically
	- take over masks in original data
- plotting:
	- plot flexible mesh net/map variables as polycollections/patches
	- plot regular grid variables with pcolor
	- plot cartopy features (land, sea, landboundary, country borders, satellite background)
	- plotting z,t-plots (see wishlist section for inaccuracies)
	- plot anything you like and how you like it
- other io functions:
	- tekal (.tek, .pli, .pliz, .pol, .ldb) data
	- read Delft3D files (.grd, .dep)
	- read/write mdu file




Example usage(#example-usage)
--------
```python
#import statements
import os
import matplotlib.pyplot as plt
plt.close('all')
from dfm_tools.get_nc import get_netdata, get_ncmodeldata, plot_netmapdata
from dfm_tools.get_nc_helpers import get_ncvardimlist, get_timesfromnc, get_hisstationlist

#define files, uncomment the line below, copy data locally and change this path to increase performance
#dir_testinput = os.path.join(r'n:\Deltabox\Bulletin\veenstra\info dfm_tools\test_input')
file_nc_map = os.path.join(dir_testinput,'DFM_sigma_curved_bend','DFM_OUTPUT_cb_3d','cb_3d_map.nc')
file_nc_his = os.path.join(dir_testinput,'DFM_sigma_curved_bend','DFM_OUTPUT_cb_3d','cb_3d_his.nc')

#get lists with vars/dims, times, station/crs/structures
vars_pd, dims_pd = get_ncvardimlist(file_nc=file_nc_map)
times_pd = get_timesfromnc(file_nc=file_nc_map)
statlist_pd = get_hisstationlist(file_nc=file_nc_his, varname='station_name')

#retrieve his data
data_fromhis_bl = get_ncmodeldata(file_nc=file_nc_his, varname='bedlevel', station='all')
data_fromhis_wl = get_ncmodeldata(file_nc=file_nc_his, varname='waterlevel', station='all', timestep= 'all')
fig, axs = plt.subplots(2,1,figsize=(10,8))
axs[0].plot(data_fromhis_bl.var_stations['station_name'],data_fromhis_bl,'-')
axs[1].plot(data_fromhis_wl.var_times,data_fromhis_wl,'-')

#retrieve net/map data, plot map data on grid
ugrid = get_netdata(file_nc=file_nc_map)#, multipart=False)
data_frommap_bl = get_ncmodeldata(file_nc=file_nc_map, varname='mesh2d_flowelem_bl')
data_frommap_sal = get_ncmodeldata(file_nc=file_nc_map, varname='mesh2d_sa1', timestep='all', layer='all')
fig, axs = plt.subplots(2,1,figsize=(6,8))
pc = plot_netmapdata(ugrid.verts, values=data_frommap_bl, ax=axs[0], linewidth=0.5, cmap='jet')
pc = plot_netmapdata(ugrid.verts, values=data_frommap_sal[0,:,-1], ax=axs[1], linewidth=0.5, cmap='jet')

#print contents of retrieved data withing data_frommap_sal variable
print_var = data_frommap_sal
print('++++++\nthe data in the variable %s is:\n%s\n'%(print_var.var_varname, print_var))
print('++++++\nthe time indices and times in the variable %s are:\n%s\n'%(print_var.var_varname, print_var.var_times))
print('++++++\nthe station indices and station names in the variable %s are:\n%s\n'%(print_var.var_varname, print_var.var_stations))
print('++++++\nthe layer indices in the variable %s are:\n%s\n'%(print_var.var_varname, print_var.var_layers))
print('++++++\nthe shape of the variable %s is:\n%s\n'%(print_var.var_varname, print_var.shape))
print('++++++\nthe netCDF variable where the data in variable %s comes from is:\n%s\n'%(print_var.var_varname, print_var.var_object))
print('++++++\nsome example contents of this netCDF variable:')
print('\tthe dimension names of the netCDF variable %s are:\n\t\t%s'%(print_var.var_varname, print_var.var_object.dimensions))
print('\tthe shape of the netCDF variable %s is:\n\t\t%s'%(print_var.var_varname, print_var.var_object.shape))
print('\tthe units of the netCDF variable %s are:\n\t\t%s'%(print_var.var_varname, print_var.var_object.units))
print('\tthe long_name of the netCDF variable %s is:\n\t\t%s'%(print_var.var_varname, print_var.var_object.long_name))
print('\tthe standard_name of the netCDF variable %s is:\n\t\t%s'%(print_var.var_varname, print_var.var_object.standard_name))
```
- for more examples, check https://github.com/openearth/dfm_tools/tree/master/tests (this is also the pytest testbank)
- examples of (mostly unformatted) figures created by this pytest testbank: n:\\Deltabox\\Bulletin\\veenstra\\info dfm_tools
- please check the TODO sections for known inaccuracies or features that are not yet available
- please report bugs and feature requests at the developers or at https://github.com/openearth/dfm_tools/issues (include OS, dfm_tools version, reproduction steps)
- want to get updates about dfm_tools? Send an email to jelmer.veenstra@deltares.nl


How to install dfm_tools
--------
- download Anaconda 64 bit Python 3.7 from https://www.anaconda.com/distribution/#download-section (miniconda is probably also sufficient, but this is not yet tested)
- install it with the recommended settings, but check 'add Anaconda3 to my PATH enviroment variable' if you want to use conda from the windows command prompt instead of anaconda prompt
- install dfm_tools from github
	- open command window (or anaconda prompt)
	- ``conda create --name dfm_tools_env python=3.7 git spyder -y`` (creating a venv is recommended, but at least do ``conda install git`` if you choose not to)
	- ``conda activate dfm_tools_env``
	- ``python -m pip install git+https://github.com/openearth/dfm_tools.git`` (this command installs all required packages and it also updates dfm_tools to the latest version if you already installed it before)
	- optional: ``conda install -c conda-forge "shapely>=1.7.0" -y`` (for slicing 2D/3D data) (conda-forge channel is necessary since main channel version is 1.6.4. The correct version is available via pip, but then geos dll is not properly linked, this will probably be solved in the future https://github.com/Toblerity/Shapely/issues/844)
	- optional: ``conda install -c conda-forge cartopy -y`` (for satellite imagery on plots) (conda-forge channel recommended by cartopy developers, and currently also necessary for correct shapely version)
- launch Spyder:
	- open 'Spyder(dfm_tools_env)' via your windows start menu (not 'Spyder' or 'Spyder(Anaconda3)', since dfm_tools was installed in dfm_tools_env)
	- test by printing dfm_tools version number: ``import dfm_tools; print(dfm_tools.__version__)`` (to double check if you are working in the venv where dfm_tools_env was installed)
	- to get figures in separate windows: go to Tools > Preferences > IPython console > Graphics > change graphics backend to 'Automatic' and restart Spyder (or the kernel).
	- check the section '[Example usage](#example-usage)' to get started


TODO wishlist
--------
- retrieve station/crs/gs list corresponding to a variable with get_hisstationlist(), now already used in stations/gs/crs check of get_nc.get_ncmodeldata()
- merge station/layer/times checks, these parts of get_nc.py have a lot of overlap
- add retrieval via depth instead of layer number (then dflowutil.mesh can be removed?):
	- refer depth w.r.t. reference level, water level or bed level
	- see test_workinprogress.py
- retrieve correct depths:
	- add depth array (interfaces/centers) to his and map variables (z/sigma layer calculation is already in get_modeldata_onintersection function)
	- depths can be retrieved from mesh2d_layer_z/mesh2d_layer_sigma, but has no time dimension so untrue for sigma and maybe for z? (wrong in dflowfm?)
	- layerzfrombedlevel keyword in mdu changes how zlayering is set up. Catch this exception with a keyword if necessary
- simplify input of modplot.velovect() for curved vectors
- improve z,t-plots from hisfile:
	- example in test_get_nc.test_gethismodeldata()
	- keep cen2cor(time_cen) definition?
	- WARNING: part of the z interfaces/center data in dflowfm hisfile is currently wrong, check your figures carefully
	- layer argument now has to be provided when retrieving zcoordinate_c (centers) from hisfile, but not when retrieving zcoordinate_w (interfaces), align this.
	- check center/corner correctness, pcolormesh does not completely correspond with contours
- improve cartopy satellite/basemap background:
	- add test if cartopy is installed before importing it, since these are optional modules (also cartopy import in user script, so does not work)
	- add more settings for linewidth/facecolor/alpha/linecolor
	- load geotiffs with satellite imagery (or png's where limits are provided by user) (files provided by user or automatically downloaded from predifined or provided source)
	- load World Imagery data from arcgis mapserver (e.g. https://www.arcgis.com/home/item.html?id=10df2279f9684e4a9f6a7f08febac2a9)
	- https://stackoverflow.com/questions/12116050/how-to-plot-remote-image-from-http-url
	- https://scitools.org.uk/cartopy/docs/v0.15/_modules/cartopy/mpl/geoaxes.html (stock_img() en background_img())
	- https://github.com/SciTools/cartopy/blob/master/lib/cartopy/data/raster/natural_earth/images.json
	- https://github.com/SciTools/cartopy/blob/master/lib/cartopy/data/raster/natural_earth/50-natural-earth-1-downsampled.png
	- http://earthpy.org/cartopy_backgroung.html
- add more io-functions:
	- convert data to kml (google earth) or shp?
	- add tekal write functions
- add tidal analysis:
	- https://github.com/sam-cox/pytides
	- https://pypi.org/project/pytides/
	- https://pypi.org/project/tidepy/
	- https://github.com/pwcazenave/tappy
	- https://pypi.org/project/UTide/
	- https://github.com/moflaher/ttide_py
- dimn_time is now actually variable name which does not work if time dimname is not the same as time varname > define whether to retrieve dim or var, or make separate definitions
- make merc keyword always optional by testing for minmax all vertsx between -181 and 361 and minmax all vertsy (lat) between -91 and 91 (+range for overlap for e.g. gtsm model)
- improve get_ncmodeldata second intersect function
	- optimize with distance from line: get maximum cell area (and infer width) from lineblockbbound selection, then decide on distance from line for selection of cells for crossect calculation
	- optimize by only retrieving necessary layerdepths/bed/waterlevel information for crossection
	- remove hardcoded (layer/)bed/waterlevel varnames
- add inpolygon/inboundbox selection of data:
	- optimize_dist keyword now draws inpolygon around line
	- to optimize intersect function when retrieving bed level and water level (do that with len(firstlinepart) optional keyword)
	- to retrieve other mapdata data faster
- add polygon ginput function (click in plot) (already partly exists in intersect/slice testscript)
- merge existing dfm model setup functions (and other useful stuff):
	- dflowutil: https://github.com/openearth/dfm_tools/tree/master/dflowutil (and test scripts, contains e.g. read/write functions for general datafromats (tim, bc))
	- MBay scripts
	- https://github.com/openearth/delft3dfmpy (arthur van dam)	
	- https://svn.oss.deltares.nl/repos/openearthtools/trunk/python/applications/delft3dfm (fiat, sobek etc)
	- https://svn.oss.deltares.nl/repos/openearthtools/trunk/python/applications/delft3dfm/dflowfmpyplot/pyd3dfm/streamline_ug.py (streamline plotting for structured grids, but many settings)
- make grid reading more flexible:
	- raise understandable error when no mesh2d_edge_x var in netcdf, instead of keyerror none (e.g. with get_netdata on hirlam files)
	- if no ugrid in netfile, try to read provided xy variables and make meshgrid or convert cen2cor or cor2cen if necessary (how to test this?)
	- improve plots for structured grid (CMEMS, ERA5, hirlam, grd etc)
	- https://github.com/NOAA-ORR-ERD/gridded (https://github.com/pyugrid/pyugrid is merged into gridded) (ghostcells en mapmergen worden afgehandeld? meer dan 4 nodes per cel? support for stations?)
	- tests.test_get_nc.test_getnetdata_plotnet_regular() is eerste opzet voor hirlam/ERA5 data, werkt heel anders dan D-flow FM
	- how to plot properties on edges/nodes (scatter is slow), maybe create dual mesh and plot like faces. most relevant variables are also available on faces, so is this necessary?
	- improve support for rstfiles (now only scatter, since only cell centers present?)
	- https://svn.oss.deltares.nl/repos/openearthtools/trunk/python/OpenEarthTools/openearthtools/io/dflowfm/patch2tri.py (equivalent van MIA)
	- https://svn.oss.deltares.nl/repos/openearthtools/trunk/python/OpenEarthTools/openearthtools/io/netcdf
	- plotting edges/nodes/faces and legend: https://github.com/pyugrid/pyugrid/blob/master/notebook_examples/connectivity_example.ipynb
	- add retrieve_netdata argument to get_ncmodeldata() which causes griddata to be retrieved as the second return variable (do this based on coodtinates / cf-conventions)?
- interactive data retrieval and plotting by calling get_ncmodeldata() without arguments


TODO non-content
--------
- make links of header-references in README.md
- request modplot.velovect() (curved vectors) to be added to matplotlib
- request shapely>=1.7.0 op main channel instead of only at conda-forge? cartopy also recommends conda-forge, so would not make a huge difference yet
- why does cartopy has to come from conda-forge?
- add variable units to plots in test bench
- readme korter maken (developer info naar aparte file?)
- update/delete cookiecutter text files (HISTORY is not up to date, remove including links in other files?)
- add documentation in comments of functions
- create overview of scripts and functions, including future location of missing features
- put testdata on deltares shared location?
- put testdata and testoutput on github and create jupyter notebook instead of pptx?
- arrange auto-testing online (jarvis?): https://docs.pytest.org/en/latest/getting-started.html
- register on PyPI, for easier install via pip (easier for regular users):
	- https://the-hitchhikers-guide-to-packaging.readthedocs.io/en/latest/quickstart.html#register-your-package-with-the-python-package-index-pypi
	- https://packaging.python.org/tutorials/packaging-projects/
	- how to automate this process? (buildserver including testing?)
	- also add changelog besides commit comments?
- update license with Deltares terms
- write documentation as comments and generate automatically?
- paths to project folders in test scripts are ok?
- style guide: https://www.python.org/dev/peps/pep-0008/
- contributing environment method: environment.yml or requirements_dev.txt?


Developer information: how to contribute to this git repository
--------
- First request github rights to contribute with the current developers
	- Jelmer Veenstra <jelmer.veenstra@deltares.nl>
	- Lora Buckman
	- Julien Groenenboom
- Get a local checkout of the github repository:
	- Download git from https://git-scm.com/download/win, install with default settings
	- open command window in a folder where you want to clone the dfm_tools github repo, e.g. C:\\DATA
	- ``git clone https://github.com/openearth/dfm_tools.git`` (repos gets cloned in C:\\DATA\\dfm_tools, this is a checkout of the master branch)
	- open git bash window in local dfm_tools folder (e.g. C:\\DATA\\dfm_tools)
	- ``git config --global user.email [emailaddress]``
	- ``git config --global user.name [username]``
	- create a branch called work_yourname on https://github.com/openearth/dfm_tools
	- open git bash window in local dfm_tools folder (e.g. C:\\DATA\\dfm_tools)
	- ``git remote update origin --prune`` (update local branch list)
	- ``git checkout work_yourname`` (checkout your branch, never do anything while the master is selected)
- Set up the dfm_tools developer python virtual environment (necessary for developing/testing):
	- open command window (or anaconda prompt) and navigate to dfm_tools folder, e.g. C:\\DATA\\dfm_tools
	- ``conda env create -f environment.yml`` (sometimes you need to press enter if it hangs extremely long)
	- to list venvs:``conda info --envs``
	- to remove venv when necessary: ``conda remove -n dfm_tools_env --all``
	- ``conda activate dfm_tools_env``
	- ``conda install spyder``
	- ``python -m pip install -e .`` (pip developer mode, any updates to the local folder by github (with ``git pull``) are immediately available in your python. It also installs all required packages)
	- ``conda install -c conda-forge "shapely>=1.7.0" cartopy``(conda-forge channel is necessary since main channel version is 1.6.4. The correct version is available via pip, but then geos dll is not properly linked, this will probably be solved in the future https://github.com/Toblerity/Shapely/issues/844. cartopy also recommends conda-forge channel)
	- test if dfm_tools is properly installed by printing the version number: ``python -c "import dfm_tools; print(dfm_tools.__version__)"``
	- open 'Spyder(dfm_tools_env)' via your windows start menu (not 'Spyder' or 'Spyder(Anaconda3)', since dfm_tools was installed in dfm_tools_env)
- Make your local changes to dfm_tools
- Work with your branch:
	- open git bash window in local dfm_tools folder (e.g. C:\\DATA\\dfm_tools)
	- ``git checkout work_yourname`` (checkout your branch, never do anything while the master is selected)
	- to update: ``git pull`` (combination of git fetch and git merge)
	- get clean checkout again (overwrite local changes):
		- ``git fetch --all`` (fetches changes)
		- ``git reset --hard`` (resets local checkout of repos branch to server version)
		- ``git pull`` (fetches and merges changes, local checkout of repos branch is now updated again)
	- ``git pull origin master`` (gets edits from master to current local branch, might induce conflicts. maybe better to just push to your branch and then handle pull request on github website)
- run test bank:
	- open command window (or anaconda prompt) in local dfm_tools folder (e.g. C:\\DATA\\dfm_tools)
	- ``conda activate dfm_tools_env``
	- ``pytest -v --tb=short`` (runs all tests)
	- ``pytest -v --tb=short -m unittest``
	- ``pytest -v --tb=short -m systemtest``
	- ``pytest -v --tb=short -m acceptance``
	- ``pytest -v --tb=short tests\test_get_nc.py::test_getplotmapWAQOS``
- Commit and push your changes to your branch:
	- open git bash window in local dfm_tools folder (e.g. C:\\DATA\\dfm_tools)
	- ``git checkout work_yourname`` (checkout your branch, never do anything while the master is selected)
	- ``git add .``
	- ``git commit -m "message to be included with your commit"``
	- ``git push`` (pushes changes to server, do not do this in while working in the master)
- increasing the version number after you committed all changes:
	- open cmd window in local dfm_tools folder (e.g. C:\\DATA\\dfm_tools)
	- optional: ``conda activate dfm_tools_env``
	- ``bumpversion major`` or ``bumpversion minor`` or ``bumpversion patch`` (changes version numbers in files and commits changes)
	- push this change in version number with ``git push`` (from git bash window or cmd also ok?)
- Request merging of your branch on https://github.com/openearth/dfm_tools/branches

