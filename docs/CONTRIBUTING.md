# Contributing

## Checkout git repository

- this is just a suggestion, feel free to work with VScode or any other git-compatible workflow
- download git from [git-scm.com](https://git-scm.com/download/win), install with default settings
- open git bash window where you want to clone the dfm_tools github repository (e.g. `C:\DATA\`)
- `git clone https://github.com/deltares/dfm_tools` (creates a folder dfm_tools with the checked out repository)
- `cd dfm_tools`
- optional: `git config --global user.email [emailaddress]`
- optional: `git config --global user.name [username]`

## Setup local developer environment

- python 3.10 to 3.13 are supported
- download Miniforge3 from [conda-forge.org](https://conda-forge.org/miniforge) and install it with the recommended settings.
- open Miniforge Prompt and navigate to the local checkout folder of the repository
- `conda create --name dfm_tools_env python=3.12 git spyder -y` (`git` and `spyder` are optional)
- `conda activate dfm_tools_env`
- `python -m pip install -e .[dev,docs,examples]` (pip developer mode, any updates to the local folder are immediately available in your python. It also installs all requirements via pip, square brackets are to install optional dependency groups)
- `conda deactivate`
- to remove dfm_tools_env when necessary: `conda remove -n dfm_tools_env --all`

## Contributing

- open an existing issue or create a new issue at [the issues page](https://github.com/Deltares/dfm_tools/issues)
- create a branch via `Development` on the right. This branch is now linked to the issue and the issue will be closed once the branch is merged with main again
- alternatively fork the repository and do your edits there
- open git bash window in local dfm_tools folder (e.g. `C:\DATA\dfm_tools`)
- `git fetch origin` followed by `git checkout [branchname]`
- make your local changes to the dfm_tools code (including docstrings and unittests for functions), after each subtask do `git commit -am 'description of what you did'` (`-am` adds all changed files to the commit)
- check if all edits were committed with `git status`, if there are new files created also do `git add [path-to-file]` and commit again
- `git push` to push your committed changes your branch on github
- open a pull request at the branch on github, there you can see what you just pushed and the automated checks will show up (testbank and code quality analysis).
- optionally make additional local changes (+commit+push) untill you are done with the issue and the automated checks have passed
- optionally increase the dfm_tools version with: `bumpversion patch`
- request a review on the pull request
- after review, squash+merge the branch into main

## Running the testbank

- open Miniforge Prompt and navigate to the local checkout folder of the repository
- `conda activate dfm_tools_env`
- `pytest` (runs all tests)
- `pytest -m "not acceptance"`
- `pytest -m acceptance` (runs the acceptance tests, which are the scripts in [the examples folder](https://github.com/Deltares/dfm_tools/tree/main/tests/examples))
- `pytest -m "not requireslocaldata"` (this is what runs on github)
- this workflow automatically runs via Github Actions upon push and pullrequest to main

## Generate documentation with sphinx

- open Miniforge Prompt and navigate to the local checkout folder of the repository
- `conda activate dfm_tools_env`
- `sphinx-build docs docs/_build`
- this workflow automatically runs via Github Actions upon push to main

## Increase the version number

- commit all changes via git
- open Miniforge Prompt and navigate to the local checkout folder of the repository
- `conda activate dfm_tools_env`
- `bumpversion major` or `bumpversion minor` or `bumpversion patch` (changes version numbers in files and commits changes)
- push changes with `git push` (from git bash window)

## Create release

- make sure the `main` branch is up to date (check pytest warnings, important issues solved, all pullrequests and branches closed)
- create and checkout branch for release
- bump the versionnumber with `bumpversion minor`
- update `docs/whats-new.md` and add a date to the current release heading
- run local testbank with `pytest -m "not requireslocaldata"`
- local check with: `python -m build` and `twine check dist/*` ([does not work on WCF](https://github.com/pypa/setuptools/issues/4133))
- commit+push to branch and merge PR
- copy the dfm_tools version from [pyproject.toml](https://github.com/Deltares/dfm_tools/blob/main/pyproject.toml) (e.g. `0.11.0`)
- create a [new release](https://github.com/Deltares/dfm_tools/releases/new)
- click `choose a tag` and type v+versionnumber (e.g. `v0.11.0`), click `create new tag: v0.11.0 on publish`
- set the release title to the tagname (e.g. `v0.11.0`)
- click `Generate release notes` and replace the `What's Changed` info by a tagged link to `docs/whats-new.md`
- if all is set, click `Publish release`
- a release is created and the github action publishes it [on PyPI](https://pypi.org/project/dfm-tools)
- post-release: commit+push `bumpversion patch` and `UNRELEASED` header in `docs/whats-new.md` to distinguish between release and dev version

## What are all these packages for?

- shapely for slicing 2D/3D data
- geopandas for shapefile related operations
- contextily for satellite imagery on plots, faster than cartopy
- xarray developers advise to install dependecies dask/netCDF4/bottleneck
- xugrid: wrapper around xarray by Huite Bootsma, for ugrid support
- cdsapi/pydap: to download ERA5 and CMEMS data

## Potential spyder issues

- Qt error upon launching Spyder?: remove the system/user environment variable 'qt_plugin_path' set by an old Delft3D4 installation procedure.
- netCDF4 DLL error upon import in Spyder?: [remove Anaconda paths from the Path user environment variable](https://github.com/spyder-ide/spyder/issues/19220)
