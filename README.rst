=========
dfm_tools
=========

dfm_tools are post-processing tools for Delft3D FM


Features
--------
- read net data
- read map and his data
- plot net data with map data
- select data based on variable, timestep, layer

TODO
--------
- move round datetimes to get_timesfromnc (should be pandas Series first)
- discuss the (future) structure: which functions in which class/script, function names, argument names?
- allowed to be freely available? add ownrisk-license
- register on PyPI, for easier install via pip (for regular users, not developers):
	- https://the-hitchhikers-guide-to-packaging.readthedocs.io/en/latest/quickstart.html#register-your-package-with-the-python-package-index-pypi 
	- also add version numbers (only master branch?)
	- also add changelog besides commit comments?
- how to plot properties on edges (scatter is slow)
- add retrieval via depth instead of layer number (then dflowutil.mesh can be removed?)
- hisinfo per station opvragen (now only all), also add requested stations to output as values_all.stations
- perform actions by dimension names instead of ndims (station_name variable has two dimensions but no time)
- add requested variablename(?), times and layers to outputdata (necessary for plotting his and more), 
- add polygon/crossection
- add inpolygon selection of data
- Dfm_tools slicing ugrid:
	- https://stackoverflow.com/questions/47868134/how-to-slice-2d-grid-from-3d-irrigular-data
	- https://github.com/pyugrid/pyugrid/tree/master/notebook_examples
	- https://stackoverflow.com/questions/15748767/interpolation-subsampling-of-3d-data-in-python-without-vtk/15753011#15753011
- write documentation as comments and generate automatically
- create testbank (keep example codelines) and setup auto-testing online (jarvis?): https://docs.pytest.org/en/latest/getting-started.html
- add comparable functions for sobek and Delft3D
- collect more functions/scripts from other users and bundle/couple


How to work with this git repository
--------
- Install Github:
	- Download git from https://git-scm.com/download/win, install with default settings
	- open command line in a folder where you want to clone the dfm_tools github repo, eg C:\DATA\GitHub
	- ``git clone https://github.com/openearth/dfm_tools.git`` (repos gets cloned to local drive, checkout of master branch)
	- to update: ``git pull`` (?)
	- NOTE: it is also possible to download the zip from https://github.com/openearth/dfm_tools, but this is not recommended since getting the updates is easier this way

- Install Python:
	- Download the newest anaconda 64 bit
	- install, including PATH checkbox

- Optional: create separate python environment and link from Spyder:
	- open command line and navigate to dfm_tools folder, eg C:\DATA\GitHub\dfm_tools
	- ``conda env create -f environment.yml`` (sometimes you need to press enter if it hangs extremely long)
	- ``conda info --envs`` (shows github_env virtual environment)
	- ``conda activate github_env``
	- ``python -c "import sys; print(sys.executable)"`` (the resulting path you need some steps later, eg C:\Users\[user]\AppData\Local\Continuum\anaconda3\envs\github_env\python.exe)
	- ``conda deactivate``
	- open spyder from start menu or anaconda or anything
	- Go to Tools >> Preferences >> Python interpreter >> point to github_env python.exe (print of sys.executable)
	- restart IPython console
	- optional: ``conda remove -n github_env --all`` (to remove it again when necessary)

- Install your local github clone via pip:
	- open command window, navigate to dfm_tools folder, eg C:\DATA\GitHub\dfm_tools
	- optional: >> activate github_env
	- >> python -m pip install -e .
	- (pip developer mode, any updates to folder by github will be available)
	- (also install all packages in requirements.txt)

- Use it in your scripts:
	- from dfm_tools.grid import get_netdata, get_hismapmodeldata, plot_netmapdata
	- check tests folder for examples


How to contribute to this git repository
--------
- First request rights to contribute
- Branching:
	- open git bash window in local dfm_tools folder (eg C:\DATA\GitHub\dfm_tools)
	- ``git config --global user.email [emailaddress]``
	- ``git config --global user.name [username]``
	- Create your own branch option 1:
		- manually create a branch on https://github.com/openearth/dfm_tools
		- open git bash window in local dfm_tools folder (eg C:\DATA\GitHub\dfm_tools)
		- ``git remote update origin --prune`` (update local branch list)
		- ``git checkout branchname`` (checkout branch)
	- Create your own branch option 2:
		- open git bash window in local dfm_tools folder (eg C:\DATA\GitHub\dfm_tools)
		- ``git checkout --branch branchname`` (create new branch and checkout, combination of git branch and git checkout commands)

- Commit and push your changes to your online branch:
	- optional: ``git pull origin master`` (gets edits from master to current local branch, might induce conflicts. maybe better to just push to your branch and then handle pull request on github website)
	- ``git add .``
	- ``git commit -m "message to be included with your commit"``
	- ``git push`` (pushes changes to server, do not do this in while working in the master)
- Request merging of your branch on https://github.com/openearth/dfm_tools/branches
