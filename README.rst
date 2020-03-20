=========
dfm_tools
=========


.. image:: https://img.shields.io/pypi/v/dfm_tools.svg
        :target: https://pypi.python.org/pypi/dfm_tools

.. image:: https://img.shields.io/travis/openearth/dfm_tools.svg
        :target: https://travis-ci.org/openearth/dfm_tools

.. image:: https://readthedocs.org/projects/dfm-tools/badge/?version=latest
        :target: https://dfm-tools.readthedocs.io/en/latest/?badge=latest
        :alt: Documentation Status

.. image:: https://pyup.io/repos/github/openearth/dfm_tools/shield.svg
        :target: https://pyup.io/repos/github/openearth/dfm_tools/
        :alt: Updates


dfm_tools are Python post-processing tools for Delft3D FM model outputfiles (netCDF) and more


* Free software: GNU General Public License v3
* Documentation: https://dfm-tools.readthedocs.io.



Features
--------
- read net data
- read map and his data
- plot net data with map data
- select all data based on variable, timestep/datetime, layer, station (not yet on depth)
- plotting zt-plots (see wishlist section for inaccuracies)
- merge partitions and delete ghostcells automatically
- take over masks in original data
- selection/plotting by polyline/crossection (slicing the ugrid data)
- pytest testbank
- examples of unformatted plots created by testbank scripts in tests folder: n:\\Deltabox\\Bulletin\\veenstra\\info dfm_tools
- please check the TODO sections for known inaccuracies or features that are not yet available
- tekal (tek, pli, pliz, pol, ldb) data
- read Delft3D files (grd, dep), although grid handling is still work in progress
- read almost any netcdf (ERA5, hirlam, SFINCS map, SFINCS his, Sobek observation), although grid handling and plotting is still work in progress

How to work with this git repository
--------
- Install Python:
	- Download the newest anaconda 64 bit: https://repo.anaconda.com/archive/Anaconda3-2019.10-Windows-x86_64.exe
	- install, including PATH checkbox

- Install the code from github via pip:
	- download and install git from https://git-scm.com/download/win
	- optional (but recommended): create and activate a separate Python virtual environment (see related information for a possible method)
	- open command window
	- ``python -m pip install git+https://github.com/openearth/dfm_tools.git`` (this also installs all required packages) (this also updates it to the latest version if you already installed it before)
	- quick test:
		- ``python -c "import dfm_tools; print(dfm_tools.__version__)"`` (print version number of the installed dfm_tools package)
		- test if you can import shapely.geometry: ``python -c "import shapely.geometry"`` (if not, look at the known bugs section in this readme. You will need this when slicing data)
	- Note: we will try to distribute dfm_tools via PyPI soon, so installation is slightly easier
	
- Use it in your scripts:
	- from dfm_tools.get_nc import get_netdata, get_ncmodeldata, plot_netmapdata
	- check scripts in tests folder on github for examples

Known bugs
--------
- the line ``import shapely.geometry`` does not work, while ``import shapely`` does (OSError: [WinError 126] The specified module could not be found), solution:
	- find geos.py in your environment (eg %userprofile%\\AppData\\Local\\Continuum\\anaconda3\\envs\\dfm_tools_env\\Lib\\site-packages\\shapely\\geos.py)
	- replace ``if os.getenv('CONDA_PREFIX', ''):`` with ``if 0:`` on line 143 (this disables this if statement and redirects to else)


TODO wishlist
--------
- select/check functions in dflowutil folder and merge with dfm_tools:
	- including dflowutil_examples/test_dflowutil.py and other test scripts
	- dflowutil contains eg readwrite functions for general datafromats (tim, bc)
- add retrieval via depth instead of layer number (then dflowutil.mesh can be removed?):
	- refer depth wrt reference level, water level or bed level, z variable is not correct in dfm-mapfile yet)
	- see test_workinprogress.py
- retrieve correct depths:
	- add depth array (interfaces/centers) to his and map variables (z/sigma layer calculation is already in get_modeldata_onintersection function)
	- depths can be retrieved from mesh2d_layer_z/mesh2d_layer_sigma, but has no time dimension so untrue for sigma and maybe for z? (wrong in dflowfm?)
	- layerzfrombedlevel keyword in mdu changes how zlayering is set up. Catch this exception with a keyword if necessary
- improve zt plots from hisfile:
	- example in test_get_nc.test_gethismodeldata()
	- WARNING: part of the z interfaces/center data in dflowfm hisfile is currently wrong, check your figures carefully
	- layer argument now has to be provided when retrieving zcoordinate_c (centers) from hisfile, but not when retrieving zcoordinate_w (interfaces), align this.
	- check center/corner correctness, pcolormesh does not completely correspond with contours
- add tekal write functions
- expand Delft3D read and plot options
- add sattelite basemap (cartopy/basemap), get latlon projection for axis
- expand general netcdf read and plot options (Sobek, ERA5, hirlam, SFINCS)
- remove hardcoded 'stations' dimension lookup
- raise understandable error when no mesh2d_edge_x var in netcdf, instead of keyerror none (eg with get_netdata on hirlam files)
- dimn_time is now actually variable name which does not work if time dimname is not the same as time varname
- make merc keyword always optional by testing for minmax all vertsx between -181 and 361 and minmax all vertsy (lat) between -91 and 91 (+range for overlap for eg gtsm model)
- optimize get_ncmodeldata for layerdepths/bedlevel/waterlevel (second intersect function), only retrieve necessary information for crossection
- add inpolygon/inboundbox selection of data:
	- optimize_dist keyword now draws inpolygon around line
	- to optimize intersect function when retrieving bed level and water level (do that with len(firstlinepart) optional keyword)
	- to retrieve other mapdata data faster
- add polygon ginput function (click in plot) (already partly exists in intersect/slice testscript)
- pyugrid (ghostcells en mapmergen worden afgehandeld?), voorbeelden in ieder geval als inspiratie voor plotopties):
	- https://github.com/pyugrid/pyugrid/blob/master/notebook_examples/COMT_example.ipynb
	- https://github.com/pyugrid/pyugrid/blob/master/notebook_examples/Delft3D%20examples.ipynb
	- https://github.com/pyugrid/pyugrid/blob/master/notebook_examples/connectivity_example.ipynb
	- https://github.com/pyugrid/pyugrid/blob/master/notebook_examples/plotting_example.ipynb
	- https://github.com/pyugrid/pyugrid/blob/master/notebook_examples/vector_plotting_example.ipynb
- make grid reading more flexible:
	- improve plots for structured grid (CMEMS, ERA5, hirlam, grd etc)
	- https://github.com/NOAA-ORR-ERD/gridded
	- tests.test_get_nc.test_gethirlam() is eerste opzet voor hirlam/ERA5 data, werkt heel anders dan D-flow FM
	- how to plot properties on edges (scatter is slow), maybe create dual mesh and plot like faces. most relevant variables are also available on faces, so is this necessary?
	- add support for rstfiles (different way of storing grid data, only face nodes present?)
	- https://svn.oss.deltares.nl/repos/openearthtools/trunk/python/OpenEarthTools/openearthtools/io/dflowfm/patch2tri.py
	- https://svn.oss.deltares.nl/repos/openearthtools/trunk/python/OpenEarthTools/openearthtools/io/netcdf
	- see test_workinprogress.py


TODO non-content
--------
- update install/venv manuals, venv manual is now not in line with user-install
- register on PyPI, for easier install via pip (easier for regular users):
	- https://the-hitchhikers-guide-to-packaging.readthedocs.io/en/latest/quickstart.html#register-your-package-with-the-python-package-index-pypi
	- https://packaging.python.org/tutorials/packaging-projects/
	- how to automate this process?
	- also add changelog besides commit comments?
- publish some example figures online, maybe py-notebook and example data?
- arrange auto-testing online (jarvis?): https://docs.pytest.org/en/latest/getting-started.html
- put testdata on deltares shared location?
- put testdata and testoutput on github and create jupyter notebook instead of pptx?
- update license with Deltares terms
- update all text files and documentations
- write documentation as comments and generate automatically
- create overview tree of all functions, also add missing functions here
- paths to project folders in test scripts are ok?
- add minimal version numbers to requirements.txt (maybe also to environment.yml)
- style guide: https://www.python.org/dev/peps/pep-0008/
- contributing method: environment.yml (README.rst) or requirements_dev.txt (CONTRIBUTING.rst)?


Related information
--------
- Create a separate python environment and link from Spyder:
	- open command line and navigate to dfm_tools github folder, eg C:\\DATA\\dfm_tools
	- ``conda env create -f environment.yml`` (sometimes you need to press enter if it hangs extremely long)
	- ``conda info --envs`` (shows dfm_tools_env virtual environment)
	- ``conda activate dfm_tools_env``
	- ``python -c "import sys; print(sys.executable)"`` (the resulting path you need some steps later, eg C:\\Users\\[user]\\AppData\\Local\\Continuum\\anaconda3\\envs\\dfm_tools_env\\python.exe)
	- ``conda deactivate``
	- open spyder from start menu or anaconda or anything
	- Go to Tools >> Preferences >> Python interpreter >> point to dfm_tools_env python.exe (print of sys.executable)
	- restart IPython console
	- optional: ``conda remove -n dfm_tools_env --all`` (to remove it again when necessary)
- how to contribute to this git repository
	- First request rights to contribute with the current developers
	- Get a local checkout of the github repository:
		- Download git from https://git-scm.com/download/win, install with default settings
		- open command line in a folder where you want to clone the dfm_tools github repo, eg C:\\DATA
		- ``git clone https://github.com/openearth/dfm_tools.git`` (repos gets cloned to local drive, checkout of master branch)
		- to update: navigate to dfm_tools folder in git bash window and ``git pull`` (combination of git fetch and git merge)
	- Install your local github clone via pip (developer mode):
		- open command window, navigate to dfm_tools folder, eg C:\\DATA\\dfm_tools
		- optional: create and activate a separate Python virtual environment (see related information for a possible method)
		- ``python -m pip install -e .`` (pip developer mode, any updates to the local folder by github (with ``git pull``) are immediately available in your python. It also installs all required packages)
		- ``python -c "import dfm_tools; print(dfm_tools.__version)"`` (print version number of the installed dfm_tools package)
	- Branching:
		- open git bash window in local dfm_tools folder (eg C:\\DATA\\dfm_tools)
		- ``git config --global user.email [emailaddress]``
		- ``git config --global user.name [username]``
		- Create your own branch option 1:
			- manually create a branch on https://github.com/openearth/dfm_tools
			- open git bash window in local dfm_tools folder (eg C:\\DATA\\dfm_tools)
			- ``git remote update origin --prune`` (update local branch list)
			- ``git checkout branchname`` (checkout branch)
		- Create your own branch option 2:
			- open git bash window in local dfm_tools folder (eg C:\\DATA\\dfm_tools)
			- ``git checkout --branch branchname`` (create new branch and checkout, combination of git branch and git checkout commands)
		- get clean checkout again (overwrite local changes):
			- ``git fetch --all`` (fetches changes)
			- ``git reset --hard`` (resets local checkout of repos branch to server version)
			- ``git pull`` (fetches and merges changes, local checkout of repos branch is now updated again)

	- Commit and push your changes to your online branch:
		- open git bash window in local dfm_tools folder (eg C:\\DATA\\dfm_tools)
		- optional: ``git pull origin master`` (gets edits from master to current local branch, might induce conflicts. maybe better to just push to your branch and then handle pull request on github website)
		- ``git add .``
		- ``git commit -m "message to be included with your commit"``
		- ``git push`` (pushes changes to server, do not do this in while working in the master)
	- increasing the version number (with bumpversion):
		- open cmd window in local dfm_tools folder (eg C:\\DATA\\dfm_tools)
		- optional: ``conda activate dfm_tools_env``
		- ``bumpversion major`` or ``bumpversion minor`` or ``bumpversion patch`` (changes version numbers in files and commits changes)
		- push your changes with ``git push`` (from git bash window or cmd also ok?)
	- Request merging of your branch on https://github.com/openearth/dfm_tools/branches
- run test bank:
	- create python virtual environment with environment.yml (developer/test dependencies are there)
	- fix the bug related to geos.py (section 'known bugs')
	- open command line in local dfm_tools folder
	- ``pytest -v --tb=short`` (runs all tests)
	- ``pytest -v --tb=short -m unittest``
	- ``pytest -v --tb=short -m systemtest``
	- ``pytest -v --tb=short -m acceptance``
	- ``pytest -v --tb=short tests\test_get_nc.py::test_getplotmapWAQOS``

Credits
-------

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage


