Information for developers
--------

Create python environment dfm_tools_devenv and install dfm_tools in developer mode:

- download and install Anaconda 64 bit Python 3.7 (or higher) from https://www.anaconda.com/distribution/#download-section (miniconda should also be sufficient, but this is not yet tested). Install it with the recommended settings, but check 'add Anaconda3 to my PATH environment variable' if you want to use conda from the windows command prompt instead of anaconda prompt
- download git from https://git-scm.com/download/win, install with default settings
- create a branch called work_yourname on https://github.com/openearth/dfm_tools
- open git bash window where you want to clone the dfm_tools github repository (e.g. C:\\DATA\\)
- optional: ``git config --global user.email [emailaddress]``
- optional: ``git config --global user.name [username]``
- optional: ``git remote update origin --prune`` (update local branch list)
- ``git clone -b work_yourname https://github.com/openearth/dfm_tools dfm_tools`` (repo gets cloned in C:\\DATA\\dfm_tools, this is a checkout of the work_yourname branch)
- update your branch if main has been updated: add+commit+push everything in branch first, ``git checkout main``, ``git pull``, ``git checkout development``, ``git merge main -m ''``, ``git push``
- open anaconda prompt and navigate to dfm_tools local folder, e.g. ``C:\\DATA\\dfm_tools``
- ``conda env create -f environment.yml`` (creates an environment called dfm_tools_devenv) TODO: yml now contains almost the same as requirements.txt, with additionally pdoc3/pytest/bump2version. Update this manual according to this
- ``conda info --envs`` (should show dfm_tools_devenv virtual environment in the list)
- ``conda activate dfm_tools_devenv``
- ``conda install -c conda-forge spyder shapely cartopy pyepsg geopandas contextily xarray dask netcdf4 bottleneck -y``
- ``python -m pip install -e .`` (pip developer mode, any updates to the local folder are immediately available in your python. It also installs all required non-conda packages)
- test if dfm_tools is properly installed by printing the version number: ``python -c "import dfm_tools; print(dfm_tools.__version__); import netCDF4"``
- ``conda deactivate``
- to remove dfm_tools_devenv when necessary: ``conda remove -n dfm_tools_devenv --all``
- open 'Spyder(dfm_tools_devenv)' via your windows start menu (not 'Spyder' or 'Spyder(Anaconda3)', since dfm_tools was installed in dfm_tools_devenv)
- Make your local changes to dfm_tools scripts

Work with your branch:

- open git bash window in local dfm_tools folder (e.g. C:\\DATA\\dfm_tools)
- ``git checkout work_yourname`` (checkout your branch, never do anything while the master is selected)
- to update: ``git pull`` (combination of git fetch and git merge)
- get clean checkout again (overwrite local changes):
	- ``git fetch --all`` (fetches changes)
	- ``git reset --hard`` (resets local checkout of repos branch to server version)
	- ``git pull`` (fetches and merges changes, local checkout of repos branch is now updated again)
- ``git pull origin master`` (gets edits from master to current local branch, might induce conflicts. maybe better to just push to your branch and then handle pull request on github website)

Running the testbank:

- open anaconda prompt in local dfm_tools folder (e.g. C:\\DATA\\dfm_tools)
- ``conda activate dfm_tools_devenv``
- ``pytest`` (runs all tests)
- ``pytest -m unittest``
- ``pytest -m systemtest``
- ``pytest -m "not acceptance"``
- ``pytest -m acceptance``(runs the acceptance tests, which are the scripts in [the examples folder](https://github.com/openearth/dfm_tools/tree/master/tests/examples))
- ``pytest -m "not slow"``
- ``pytest tests\test_get_nc.py::test_getplotmapWAQOS``
- the following arguments are automatically provided via pytest.ini: ``-v --tb=short``, add ``--cov=dfm_tools`` for a coverage summary

Regenerate html documentation:

- open anaconda prompt in local dfm_tools folder (e.g. C:\\DATA\\dfm_tools)
- ``conda activate dfm_tools_devenv``
- ``pdoc --html dfm_tools -o docs --force``

Commit and push your changes to your branch:

- open git bash window in local dfm_tools folder (e.g. C:\\DATA\\dfm_tools)
- ``git checkout work_yourname`` (checkout your branch, never do anything while the master is selected)
- ``git add .``
- ``git commit -m "message to be included with your commit"``
- ``git push`` (pushes changes to server, do not do this in while working in the master)

Increasing the dfm_tools version number:

- commit all changes via git
- open anaconda prompt in local dfm_tools folder (e.g. C:\\DATA\\dfm_tools)
- optional?: ``conda activate dfm_tools_devenv``
- ``bumpversion major`` or ``bumpversion minor`` or ``bumpversion patch`` (changes version numbers in files and commits changes)
- push this change in version number with ``git push`` (from git bash window or cmd also ok?)
- request merging of your branch on https://github.com/openearth/dfm_tools/branches
