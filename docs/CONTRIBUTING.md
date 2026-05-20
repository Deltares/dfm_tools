# Developer's guide

With the following steps you can install dfm_tools as a developer so you have access to latest version on github and can make changes to the code.

## Install

First, install Pixi: https://pixi.prefix.dev/latest/installation

Then, clone the dfm_tools ``git`` repo from [github](https://github.com/Deltares/dfm_tools), then navigate into the code folder (where the pyproject.toml is located):
```
git clone https://github.com/Deltares/dfm_tools.git
cd dfm_tools
```

If git is not recognized as a command, first install Git from https://git-scm.com/install. VSCode and PyCharm should be bundled with a git extension so a manual installation is not always necessary.

Then, create and activate a new pixi environment. This includes a developers installation of dfm_tools and the default environment as specified in pyproject.toml:

```
pixi install
```

Any other environments will be installed/updated when you use them for instance when running the tests.

## Running the tests

Running the tests can be done within any of the environment specified in the pyproject.toml (defaults to `default`):

```
pixi run -e default test
pixi run -e default pytest
```

Most of the time it is smart to skip the acceptance tests (the scripts in [the examples folder](https://github.com/Deltares/dfm_tools/tree/main/tests/examples)), since these are quite slow:
```
pixi run test -m "not acceptance"
```

Or even better, to exclude all tests that require local (P-drive) data (the examples/acceptancetests and a bit more tests), which might also be slow. This is also the selection of tests that runs on github:
```
pixi run test -m "not requireslocaldata"
```

Additionally, it might be smart to also skip the slow ERA5 tests since they will crash the entire test workflow on windows when taking too long:
```
pixi run test -m "not requireslocaldata and not era5slow"
```

## Updating the lockfile

If you add any dependencies or change anything in the package configuration, you have to update the lockfile. This is also done (if needed) when running other pixi commands:

```
pixi lock
```

## Generating the docs

Generating the docs is added as a pixi task, which executes sphinx-build:

```
pixi run docs-build
```

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

## Increase the version number

- commit all changes via git
- `pixi run bumpversion major` (or `minor` or `patch`) (changes version numbers in files and commits changes)
- push changes with `git push` (from git bash window/terminal)

## Create release

- make sure the `main` branch is up to date (check pytest warnings, important issues solved, all pullrequests and branches closed)
- create and checkout branch for release
- bump the versionnumber with `pixi run bumpversion minor`
- update the lockfile with `pixi lock`
- update `docs/whats-new.md` and add a date to the current release heading
- run local testbank with `pixi run pytest -m "not requireslocaldata"`
- local check with: `pixi run python -m build` and `pixi run twine check dist/*` ([does not work on WCF](https://github.com/pypa/setuptools/issues/4133))
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
- xarray developers advise to install dependecies dask/netCDF4/bottleneck for io and performance
- xugrid: wrapper around xarray by Huite Bootsma, for ugrid support
- cdsapi/pydap: to download ERA5 and CMEMS data
