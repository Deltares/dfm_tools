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


dfm_tools are post-processing tools for Delft3D FM model outputfiles and more


* Free software: GNU General Public License v3
* Documentation: https://dfm-tools.readthedocs.io.



Features
--------
- read net data
- read map and his data
- plot net data with map data
- select all data based on variable, timestep/datetime, layer, station (not yet on depth)
- merge partitions and delete ghostcells automatically
- take over masks in original data
- selection/plotting by polyline/crossection (slicing the ugrid data)
- pytest testbank

Terms of use
--------
- this toolbox is now only available via github, it will soon be registered on pypi so pip installing and updating is possible without github checkout
- it might be that function names and arguments names will slightly change in the upcoming weeks. This would mean the script you create now might need some edits to work with the toolbox in the future. this is mainly the case for the get_modeldata_onintersection() function
- please check the TODO sections for known inaccuracies
- please do not use the dflowutil/dflowutil_examples scripts if you did not before, this will be phased out eventually. All you probably need is dfm_tools

Known bugs
--------
- the line ``import shapely.geometry`` does not work, while ``import shapely`` does, solution:
	- find geos.py in your environment (eg %userprofile%\\AppData\\Local\\Continuum\\anaconda3\\envs\\github_env\\Lib\\site-packages\\shapely\\geos.py)
	- replace ``if os.getenv('CONDA_PREFIX', ''):`` with ``if 0:`` on line 143
	
How to work with this git repository
--------
- Install Github:
	- Download git from https://git-scm.com/download/win, install with default settings
	- open command line in a folder where you want to clone the dfm_tools github repo, eg C:\\DATA\\GitHub
	- ``git clone https://github.com/openearth/dfm_tools.git`` (repos gets cloned to local drive, checkout of master branch)
	- to update: navigate to dfm_tools folder in git bash window and ``git pull`` (combination of git fetch and git merge)
	- NOTE: in the near future (hopefully within a week), this package should be installable via pip, after registering on PyPI. then users do not need github anymore, only developers do

- Install Python:
	- Download the newest anaconda 64 bit
	- install, including PATH checkbox

- Install your local github clone via pip (developer mode):
	- open command window, navigate to dfm_tools folder, eg C:\\DATA\\GitHub\\dfm_tools
	- optional: create and activate a separate Python virtual environment (see related information for a possible method)
	- ``python -m pip install -e .`` (pip developer mode, any updates to the local folder by github (with ``git pull``) are immediately available in your python. It also installs all required packages)
	- ``python -c "import dfm_tools; print(dfm_tools.__version)"`` (print version number of the installed dfm_tools package)
	
- Use it in your scripts:
	- from dfm_tools.get_nc import get_netdata, get_ncmodeldata, plot_netmapdata
	- check scripts in tests folder on github for examples

TODO high priority (before launch)
--------
- register on PyPI, for easier install via pip (for regular users, not developers):
	- https://the-hitchhikers-guide-to-packaging.readthedocs.io/en/latest/quickstart.html#register-your-package-with-the-python-package-index-pypi 
	- also add version numbers (only master branch), git commit automatic minor numbers? (bumpversion comes with cookiecutter?)
	- also add changelog besides commit comments?
- exclude dflowutil from pypi package?

TODO
--------
- update license with Deltares terms
- paths to project folders in test scripts are ok?
- add retrieval via depth instead of layer number (then dflowutil.mesh can be removed?) (refer depth wrt reference level, water level or bed level, z variable is not correct in dfm-mapfile yet)
- retrieve correct depths:
	- add depth array (interfaces/centers) to his and map variables (z/sigma layer calculation is already in get_modeldata_onintersection function)
	- depths can be retrieved from mesh2d_layer_z/mesh2d_layer_sigma, but has no time dimension so untrue for sigma and maybe for z? (wrong in dflowfm?)
	- layerzfrombedlevel keyword in mdu changes how zlayering is set up. Catch this exception with a keyword if necessary
- remove hardcoded 'stations' dimension lookup
- dimn_time is now actually variable name which does not work if time dimname is not the same as time varname
- contributing method: environment.yml (README.rst) or requirements_dev.txt (CONTRIBUTING.rst)?
- perform actions by dimension names instead of ndims (eg station_name variable has two dimensions but no time)
- make merc keyword always optional by testing for minmax all vertsx between -181 and 361 and minmax all vertsy (lat) between -91 and 91 (+range for overlap for eg gtsm model)
- optimize get_ncmodeldata for layerdepths/bedlevel/waterlevel (second intersect function), only retrieve necessary information for crossection
- add inpolygon/inboundbox selection of data:
	- to optimize intersect function when retrieving bed level and water level (do that with len(firstlinepart) optional keyword)
	- to retrieve other mapdata data faster
	- https://stackoverflow.com/questions/31542843/inpolygon-for-python-examples-of-matplotlib-path-path-contains-points-method
- make patched zt plots from hisfile (careful, z interfaces data in hisfile is wrong)
- as user: get stationlist, dimensionlist, variablelist, more? (partly internally available)
- add polygon read/write function (also ldb files)
- add polygon ginput function (click in plot) (already partly exists in intersect/slice testscript)
- style guide: https://www.python.org/dev/peps/pep-0008/
- pyugrid (ghostcells en mapmergen worden afgehandeld?), voorbeelden in ieder geval als inspiratie voor plotopties):
	- https://github.com/pyugrid/pyugrid/blob/master/notebook_examples/Delft3D%20examples.ipynb
	- https://github.com/pyugrid/pyugrid/blob/master/notebook_examples/connectivity_example.ipynb
	- https://github.com/pyugrid/pyugrid/blob/master/notebook_examples/plotting_example.ipynb
	- https://github.com/pyugrid/pyugrid/blob/master/notebook_examples/vector_plotting_example.ipynb
- any grid: https://github.com/NOAA-ORR-ERD/gridded
- how to plot properties on edges (scatter is slow), maybe create dual mesh and plot like faces. most relevant variables are also available on faces, so is this necessary?
- add (look for) readwrite functions for general datafromats (tim, tekal etc)
- add plot of structured grid (CMEMS etc)
- add foufiles, rstfiles? (partitioned but with different dimensions, should already partially work)
- add minimal version numbers to requirements.txt (maybe also to environment.yml)
- create overview tree of all functions, also add missing functions here
- write documentation as comments and generate automatically
- improve testbank:
	- parametrize test_grid_gethismodeldata
	- arrange auto-testing online (jarvis?): https://docs.pytest.org/en/latest/getting-started.html
- add comparable functions for sobek and Delft3D

Related information
--------
- Create a separate python environment and link from Spyder:
	- open command line and navigate to dfm_tools folder, eg C:\\DATA\\GitHub\\dfm_tools
	- ``conda env create -f environment.yml`` (sometimes you need to press enter if it hangs extremely long)
	- ``conda info --envs`` (shows github_env virtual environment)
	- ``conda activate github_env``
	- ``python -c "import sys; print(sys.executable)"`` (the resulting path you need some steps later, eg C:\\Users\\[user]\\AppData\\Local\\Continuum\\anaconda3\\envs\\github_env\\python.exe)
	- ``conda deactivate``
	- open spyder from start menu or anaconda or anything
	- Go to Tools >> Preferences >> Python interpreter >> point to github_env python.exe (print of sys.executable)
	- restart IPython console
	- optional: ``conda remove -n github_env --all`` (to remove it again when necessary)
- how to contribute to this git repository
	- First request rights to contribute
	- Branching:
		- open git bash window in local dfm_tools folder (eg C:\\DATA\\GitHub\\dfm_tools)
		- ``git config --global user.email [emailaddress]``
		- ``git config --global user.name [username]``
		- Create your own branch option 1:
			- manually create a branch on https://github.com/openearth/dfm_tools
			- open git bash window in local dfm_tools folder (eg C:\\DATA\\GitHub\\dfm_tools)
			- ``git remote update origin --prune`` (update local branch list)
			- ``git checkout branchname`` (checkout branch)
		- Create your own branch option 2:
			- open git bash window in local dfm_tools folder (eg C:\\DATA\\GitHub\\dfm_tools)
			- ``git checkout --branch branchname`` (create new branch and checkout, combination of git branch and git checkout commands)
		- get clean checkout again (overwrite local changes):
			- ``git fetch --all`` (fetches changes)
			- ``git reset --hard`` (resets local checkout of repos branch to server version)
			- ``git pull`` (fetches and merges changes, local checkout of repos branch is now updated again)

	- Commit and push your changes to your online branch:
		- optional: ``git pull origin master`` (gets edits from master to current local branch, might induce conflicts. maybe better to just push to your branch and then handle pull request on github website)
		- ``git add .``
		- ``git commit -m "message to be included with your commit"``
		- ``git push`` (pushes changes to server, do not do this in while working in the master)
	- Request merging of your branch on https://github.com/openearth/dfm_tools/branches
- run test bank:
	- create python virtual environment with environment.yml (developer/test dependencies are there)
	- fix the bug related to geos.py (section 'known bugs')
	- open command line in local dfm_tools folder
	- ``pytest -v --tb=short`` (runs all tests)
	- ``pytest -v --tb=short -m unittest``
	- ``pytest -v --tb=short -m systemtest``
	- ``pytest -v --tb=short -m acceptance``
	- ``pytest -v --tb=short tests\test_grid.py::test_mapOS``

Credits
-------

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage


