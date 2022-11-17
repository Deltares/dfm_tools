[![generate-documentation](https://github.com/openearth/dfm_tools/actions/workflows/generate-documentation.yml/badge.svg)](https://github.com/openearth/dfm_tools/actions/workflows/generate-documentation.yml)

dfm_tools
=========

A Python package for pre- and postprocessing D-FlowFM model input and output files. Contains convenience functions built on top of other packages like [xarray](https://github.com/pydata/xarray), [hydrolib-core](https://github.com/Deltares/HYDROLIB-core) and many more.

Information and examples
--------
- [pdf](https://nbviewer.org/github/openearth/dfm_tools/raw/pptx/docs/dfm_tools.pdf?flush_cache=true) with dfm_tools information, features and examples
- [online documentation](https://htmlpreview.github.io/?https://github.com/openearth/dfm_tools/blob/master/docs/dfm_tools/index.html) generated from docstrings
- [jupyter notebooks](https://github.com/openearth/dfm_tools/blob/master/notebooks) with example code
- [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/openearth/dfm_tools/HEAD) to run these notebooks interactively (loading takes a while)
- [github folder](https://github.com/openearth/dfm_tools/tree/master/tests/examples) with more example scripts


Installation
--------
- download and install Anaconda 64 bit (with Python 3.8 or later) from https://www.anaconda.com/distribution/#download-section
- open Anaconda prompt
- ``conda create --name dfm_tools_env -c conda-forge python=3.8 spyder -y`` (you can also install a newer python version)
- ``conda activate dfm_tools_env``
- ``conda install -c conda-forge git shapely cartopy pyepsg geopandas contextily xarray dask netcdf4 bottleneck cdsapi pydap -y`` (installs conda-forge requirements)
- ``python -m pip install git+https://github.com/openearth/dfm_tools`` (this command installs dfm_tools and all required non-conda packages, also use to update)
- long paths error? Check last comment in https://github.com/Deltares/HYDROLIB-core/issues/327
- to remove environment when necessary: ``conda remove -n dfm_tools_env --all``
