#Introduction

There are a few dfm_tools tutorials available as [Jupyter Notebooks](https://github.com/Deltares/dfm_tools/blob/main/docs/notebooks) which are embedded on the following pages. There are also [example scripts](https://github.com/Deltares/dfm_tools/tree/main/tests/examples) available with more specific applications. The [pdf slides](https://nbviewer.org/github/Deltares/dfm_tools/raw/main/docs/dfm_tools.pdf) contain an overview of dfm_tools information, features and examples.

dfm_tools is built on top of other well documentated packages like [xarray](https://docs.xarray.dev/en/stable/getting-started-guide/quick-overview.html), [xugrid](https://deltares.github.io/xugrid/user_guide.html), [pandas](https://pandas.pydata.org/docs/getting_started/index.html), [geopandas](https://geopandas.org/en/stable/getting_started/introduction.html) and [HYDROLIB-core](https://deltares.github.io/HYDROLIB-core). Checking the `type()` of a variable will help to find the relevant documentation.

#Pre-processing
For pre-processing (e.g. the modelbuilder), the focus is often on data conversion. Raw data is being read with pandas, geopdandas or xarray, then processed (e.g. interpolated) and then written to model input files with xarray or HYDROLIB-core. For mesh generation, the [MeshKernelPy](https://deltares.github.io/MeshKernelPy/examples/index.html) package is used.

#Post-processing
For post-processing, the main sources of information are xarray and xugrid. For instance, the function `dfmt.open_partitioned_dataset()` returns a variable of type xugrid.UgridDataset. The [xugrid user guide](https://deltares.github.io/xugrid/user_guide.html) contains many useful examples. The xugrid package in its turn wraps the underlying data as xarray datasets. This is a powerful package for lazy loading netcdf data (among others) and performing delayed operations on them. If you are unfamiliar with it, please read the [xarray in 45 minutes tutorial](https://tutorial.xarray.dev/overview/xarray-in-45-min.html).
