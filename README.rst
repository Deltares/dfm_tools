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
- select data based on datetime (now only timestep/index): get_timeid_fromdatetime
- add retrieval via depth instead of layer number (then dflowutil.mesh can be removed?)
- hisinfo per station opvragen (now only all)
- add requested variablename(?), times and layers to outputdata (necessary for plotting his and more), also multiple variables with different dimensions? (class?)     
- add polygon/crossection
- add inpolygon selection of data
- Dfm_tools slicing ugrid:
	https://stackoverflow.com/questions/47868134/how-to-slice-2d-grid-from-3d-irrigular-data
	https://github.com/pyugrid/pyugrid/tree/master/notebook_examples
	https://stackoverflow.com/questions/15748767/interpolation-subsampling-of-3d-data-in-python-without-vtk/15753011#15753011
- documentatie schrijven als comments en automatisch laten genereren
- testbank maken en online zetten

How to work with this git repository
--------
INSTALL GITHUB:
- Download git from https://git-scm.com/download/win, install with default settings
- open command line in a folder where you want to clone the dfm_tools github repo, eg C:\DATA\GitHub
- git clone https://github.com/openearth/dfm_tools.git (repos gets cloned to local drive, checkout of master branch)
- to update: ??

INSTALL PYTHON:
- Download the newest anaconda 64 bit
- install, including PATH checkbox

OPTIONAL: CREATE SEPARATE PYTHON ENVIRONMENT AND LINK FROM SPYDER:
- open command line and navigate to dfm_tools folder, eg C:\DATA\GitHub\dfm_tools
- ``conda env create -f environment.yml``
(sometimes you need to press enter if it hangs extremely long)
- >> conda info --envs (shows github_env virtual environment)
>> conda activate github_env
	python -c "import sys; print(sys.executable)"
	(the resulting path you need some steps later)
- >> conda deactivate
- open spyder from start menu or anaconda or anything
- Go to Tools >> Preferences >> Python interpreter >> point to github_env python.exe (print of sys.executable)
- restart IPython console
- >> conda remove -n github_env --all
- (to remove it again when necessary)

LINK GITHUB CLONE VIA PIP:
- open command window, navigate to dfm_tools folder, eg C:\DATA\GitHub\dfm_tools
- optional: >> activate github_env
- >> python -m pip install -e .
- (pip developer mode, any updates to folder by github will be available)
- (also install all packages in requirements.txt)

USE IT IN YOUR SCRIPTS:
- from dfm_tools.grid import get_netdata, get_hismapmodeldata, plot_netmapdata
- check tests folder for examples

How to contribute with this git repository (brancing, committing and pushing edits)
--------
BRANCHING:
- First request rights to contribute
- open GIT bash window
- navigate to local dfm_tools folder (not to folder above, then git does not work)
- git config --global user.email [emailaddress]
- git config --global user.name [username]
- Create your own branch option 1:
	manually create a branch on website
	git remote update origin --prune
	(update local branch list)
	git checkout branchname
	(checkout branch)
- Create your own branch 2:
	git checkout --branch branchname
	(create new branch and checkout, combination of git branch and git checkout commands)

REALLY COMMIT AND PUSH:
- Optional: git pull origin master
- (gets edits from master to current local branch, might induce conflicts. maybe better to just push to your branch and then handle pull request on github website)
- git add .
- git commit -m "message to be included with your commit"
- git push
- (pushes changes to server, do not do this in while working in the master)
- (local changes are now also visible under branchname on github, there you can request merging with master)

