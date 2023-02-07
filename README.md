[![generate-documentation](https://github.com/Deltares/dfm_tools/actions/workflows/generate-documentation.yml/badge.svg)](https://github.com/Deltares/dfm_tools/actions/workflows/generate-documentation.yml)
[![Quality Gate Status](https://sonarcloud.io/api/project_badges/measure?project=Deltares_dfm_tools&metric=alert_status)](https://sonarcloud.io/summary/new_code?id=Deltares_dfm_tools)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/Deltares/dfm_tools/HEAD)

dfm_tools
=========

A Python package for pre- and postprocessing D-FlowFM model input and output files. Contains convenience functions built on top of other packages like [xarray](https://github.com/pydata/xarray), [xugrid](https://github.com/Deltares/xugrid), [hydrolib-core](https://github.com/Deltares/HYDROLIB-core) and many more.

Information and examples
--------
- [pdf](https://nbviewer.org/github/Deltares/dfm_tools/raw/pptx/docs/dfm_tools.pdf?flush_cache=true) with dfm_tools information, features and examples
- [online documentation](https://htmlpreview.github.io/?https://github.com/Deltares/dfm_tools/blob/main/docs/dfm_tools/index.html) generated from docstrings
- [jupyter notebooks](https://github.com/Deltares/dfm_tools/blob/main/notebooks) with example code
- [use binder](https://mybinder.org/v2/gh/Deltares/dfm_tools/HEAD) to run these notebooks interactively (loading takes a while)
- [github folder](https://github.com/Deltares/dfm_tools/tree/main/tests/examples) with more example scripts


Installation
--------
- download and install Anaconda 64 bit (with Python 3.8 or later) from https://www.anaconda.com/distribution/#download-section
- open Anaconda prompt
- ``conda create --name dfm_tools_env -c conda-forge python=3.8 spyder -y`` (you can also install a newer python version)
- ``conda activate dfm_tools_env``
- ``conda install -c conda-forge git shapely cartopy pyepsg geopandas contextily xarray dask netcdf4 bottleneck xugrid cdsapi pydap -y`` (installs conda-forge requirements)
- ``python -m pip install git+https://github.com/Deltares/dfm_tools`` (this command installs dfm_tools and all required non-conda packages, also use to update)
- long paths error? Check [this Github issue](https://github.com/Deltares/HYDROLIB-core/issues/327#issuecomment-1266534032)
- OpenSSL error? Fix your conda base env by doing [this](https://github.com/conda/conda/issues/11795#issuecomment-1335666474) or maybe [this](https://github.com/conda/conda/issues/11982#issuecomment-1285538983). Let us know if you encounter this issue.
- to remove environment when necessary: ``conda remove -n dfm_tools_env --all``
