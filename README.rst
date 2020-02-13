=========
dfm_tools
=========

dfm_tools are post-processing tools for Delft3D FM


Features
--------
- read net data
- read map and his data
- plot net data with map data
- select all data based on variable, timestep/datetime, layer, station (not yet on depth)
- merge partitions and delete ghostcells automatically
- take over masks in original data
- selection/plotting by polygon/crossection, so slicing the ugrid data (distance calculation based on coordinates, latlon will be off but plan is to implement correct conversion, testcase needed)

License
--------
- this toolbox is free to use, but completely at your own risk (a proper licence file will be added soon)
- this toolbox is now only available via github, once the structure is decided upon, it will be registered on pypi so pip installing is possible without checkout
- there is no warranty whatsoever, users are responsible to check the output themselves
- it might be that script names, function names, argument names and more will slightly change in the upcoming weeks. This would mean the script you create now might need some edits to work with the toolbox in the future.

TODO high priority (before launch)
--------
- discuss dfm_tools structure:
	- ugrid class naar ugrid functie en rest van functies in get_dfm script? Intersect naar los script of bij get_dfm?
	- which functions in which class/script
	- function names (eg rename get_hismapmodeldata to get_netcdfdata)
	- argument names (eg rename lay?)
- add ownrisk-license
- register on PyPI, for easier install via pip (for regular users, not developers):
	- https://the-hitchhikers-guide-to-packaging.readthedocs.io/en/latest/quickstart.html#register-your-package-with-the-python-package-index-pypi 
	- also add version numbers (only master branch), git commit automatic numbers?
	- also add changelog besides commit comments?

TODO high priority
--------
- convert latlon to merc in intersect function, testcase needed
- perform actions by dimension names instead of ndims (eg station_name variable has two dimensions but no time)
- add retrieval via depth instead of layer number (then dflowutil.mesh can be removed?) (refer depth wrt reference level, water level or bed level, z variable is not correct in dfm-mapfile yet)
- retrieve correct depths:
	- add depth array (interfaces/centers) to his and map variables (z/sigma layer calculation is already in get_modeldata_onintersection function)
	- depths can be retrieved from mesh2d_layer_z/mesh2d_layer_sigma, but has no time dimension so untrue for sigma and maybe for z? (wrong in dflowfm?)
	- layerzfrombedlevel keyword in mdu changes how zlayering is set up. Catch this exception with a keyword if necessary

TODO longterm
--------
- make patched zt plots from hisfile (careful, z interfaces data in hisfile is wrong)
- as user: get stationlist, dimensionlist, variablelist, more? (partly internally available)
- add minimal version numbers to requirements.txt
- add polygon read/write function, add ginput polygon function (click in plot) (already partly exists in intersect/slice testscript)
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
- add foufiles, rstfiles? (partitioned but with different dimensions)
- add inpolygon selection of data
- create overview tree of all functions, also add missing functions here
- write documentation as comments and generate automatically
- create testbank (keep example codelines) and setup auto-testing online (jarvis?): https://docs.pytest.org/en/latest/getting-started.html
- add comparable functions for sobek and Delft3D

How to work with this git repository
--------
- Install Github:
	- Download git from https://git-scm.com/download/win, install with default settings
	- open command line in a folder where you want to clone the dfm_tools github repo, eg C:\\DATA\\GitHub
	- ``git clone https://github.com/openearth/dfm_tools.git`` (repos gets cloned to local drive, checkout of master branch)
	- to update: navigate to dfm_tools folder and ``git pull``
	- NOTE: it is also possible to download the zip from https://github.com/openearth/dfm_tools, but this is not recommended since getting the updates is easier this way
	- NOTE: in the near future (hopefully within a week), this package should be installable via pip, after registering on PyPI. then users do not need github anymore, only developers do

- Install Python:
	- Download the newest anaconda 64 bit
	- install, including PATH checkbox

- Install your local github clone via pip (developer mode):
	- open command window, navigate to dfm_tools folder, eg C:\\DATA\\GitHub\\dfm_tools
	- optional: create and activate a separate Python virtual environment
	- ``python -m pip install -e .``
	- (pip developer mode, any updates to folder by github will be available)
	- (also install all packages in requirements.txt)

- Use it in your scripts:
	- from dfm_tools.grid import get_netdata, get_hismapmodeldata, plot_netmapdata
	- check tests folder for examples


How to contribute to this git repository
--------
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

- Commit and push your changes to your online branch:
	- optional: ``git pull origin master`` (gets edits from master to current local branch, might induce conflicts. maybe better to just push to your branch and then handle pull request on github website)
	- ``git add .``
	- ``git commit -m "message to be included with your commit"``
	- ``git push`` (pushes changes to server, do not do this in while working in the master)
- Request merging of your branch on https://github.com/openearth/dfm_tools/branches

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
