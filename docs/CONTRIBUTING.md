# Information for developers

## Setup local developer environment

- download and install Anaconda 64 bit Python 3.8 (or higher) from [anaconda.com](https://www.anaconda.com/distribution/#download-section) (miniconda should also be sufficient, but this is not yet tested). Install it with the recommended settings.
- download git from [git-scm.com](https://git-scm.com/download/win), install with default settings
- open git bash window where you want to clone the dfm_tools github repository (e.g. ``C:\DATA\``)
- git clone https://github.com/deltares/dfm_tools (creates a folder dfm_tools with the checked out repository)
- ``cd dfm_tools``
- optional: ``git config --global user.email [emailaddress]``
- optional: ``git config --global user.name [username]``
- open anaconda prompt and navigate to dfm_tools local folder, e.g. ``C:\DATA\dfm_tools``
- ``conda env create -f environment.yml`` (creates an environment called dfm_tools_env) TODO: yml now contains almost the same as requirements.txt. Update this manual according to this
- ``conda activate dfm_tools_env``
- optional: ``conda install spyder -c conda-forge`` (installs spyder, using conda-forge channel since it was also used in environment.yml)
- ``python -m pip install -r requirements_dev.txt`` (installs pytest and other developer dependencies)
- ``python -m pip install -e .`` (pip developer mode, any updates to the local folder are immediately available in your python. It also installs all required non-conda packages) >> maybe add ``test`` to install also test requirements [like this](https://stackoverflow.com/questions/15422527/best-practices-how-do-you-list-required-dependencies-in-your-setup-py)
- ``conda deactivate``
- to remove dfm_tools_env when necessary: ``conda remove -n dfm_tools_env --all``

## Contributing

- open an existing issue or create a new issue at https://github.com/Deltares/dfm_tools/issues
- create a branch via ``Development`` on the right. This branch is now linked to the issue and the issue will be closed once the branch is merged with main again
- open git bash window in local dfm_tools folder (e.g. ``C:\DATA\dfm_tools``)
- ``git fetch origin`` followed by ``git checkout [branchname]``
- make your local changes to the dfm_tools code (also create unittests), after each subtask do ``git commit -am 'description of what you did'`` (``-am`` adds all changed files to the commit)
- check if all edits were committed with ``git status``, if there are new files created also do ``git add [path-to-file]`` and commit again
- ``git push`` to push your committed changes your branch on github
- open a pull request at the branch on github, there you can see what you just pushed and the automated checks will show up (testbank and code quality analysis).
- optionally make additional local changes (+commit+push) untill you are done with the issue and the automated checks have passed
- optionally increase the dfm_tools version with: ``bumpversion patch``
- request a review on the pull request

## Running the testbank (also partly runs on github automatically)

- open anaconda prompt in local dfm_tools folder (e.g. ``C:\DATA\dfm_tools``)
- ``conda activate dfm_tools_env``
- ``pytest`` (runs all tests)
- ``pytest -m "not acceptance"``
- ``pytest -m acceptance`` (runs the acceptance tests, which are the scripts in [the examples folder](https://github.com/Deltares/dfm_tools/tree/main/tests/examples)) and [the examples_workinprogress folder](https://github.com/Deltares/dfm_tools/tree/main/tests/examples_workinprogress))
- ``pytest -m "not requireslocaldata"`` (this is what runs on github)

## Generate documentation with mkdocs (automatically runs via Github Actions upon push to main)

- open anaconda prompt in local dfm_tools folder (e.g. ``C:\DATA\dfm_tools``)
- ``conda activate dfm_tools_env``
```
cp README.md docs
mkdocs build
```

## Increase the dfm_tools version number

- commit all changes via git
- open anaconda prompt in local dfm_tools folder (e.g. ``C:\DATA\dfm_tools``)
- ``conda activate dfm_tools_env``
- ``bumpversion major`` or ``bumpversion minor`` or ``bumpversion patch`` (changes version numbers in files and commits changes)
- push changes with ``git push`` (from git bash window)

## Create release

- make sure the ``main`` branch is up to date (important issues solved, all pullrequests closed, the versionnumber is correct)
- copy the dfm_tools version from https://github.com/Deltares/dfm_tools/blob/main/setup.cfg (e.g. ``0.11.0``)
- go to https://github.com/Deltares/dfm_tools/releases/new
- click ``choose a tag`` and type v+versionnumber (e.g. ``v0.11.0``), click ``create new tag: v0.11.0 on publish``
- set the release title to the tagname (e.g. ``v0.11.0``)
- click ``Generate release notes`` and check if there is anything that should be removed
- if all is set, click ``Publish release``
- a release is created and the github action publishes it on PyPI (https://pypi.org/project/dfm-tools/)

## What are all these packages for?

- shapely for slicing 2D/3D data
- cartopy for satellite imagery, coastlines etc on plots (can only be installed via conda)
- pyepsg is necessary for cartopy and probably also for other libraries
- geopandas for shapefile related operations
- contextily for satellite imagery on plots, seems faster than cartopy
- xarray developers advise to install dependecies dask/netCDF4/bottleneck with conda-forge also: https://docs.xarray.dev/en/v0.8.0/installing.html
- xugrid: wrapper around xarray by Huite Bootsma, for ugrid support
- cdsapi/pydap: to download ERA5 and CMEMS data. Minimal pydap version is 3.3.0?

## Potential spyder issues

- Qt error upon launching Spyder?: remove the system/user environment variable 'qt_plugin_path' set by an old Delft3D4 installation procedure.
- netCDF4 DLL error upon import in Spyder?: remove Anaconda paths from the Path user environment variable (https://github.com/spyder-ide/spyder/issues/19220)
