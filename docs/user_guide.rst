User Guide
==========

This user guide is example based, and should give you an idea how to accomplish common tasks.

Information on specific methods and classes can be found in the API Reference.

There are a few tutorials available as `Jupyter Notebooks <https://github.com/Deltares/dfm_tools/blob/main/docs/notebooks>`_, which are embedded here. You can use `Binder <https://mybinder.org/v2/gh/Deltares/dfm_tools/HEAD?urlpath=/tree/docs/notebooks>`_ to run these notebooks interactively (loading can take a while). There are also `example scripts <https://github.com/Deltares/dfm_tools/tree/main/tests/examples>`_ available with more specific applications. The `pdf slides <https://nbviewer.org/github/Deltares/dfm_tools/raw/main/docs/dfm_tools.pdf>`_ contain an overview of dfm_tools information, features and examples.


dfm_tools is built on top of other well documentated packages like `xarray <https://docs.xarray.dev/en/stable/getting-started-guide/quick-overview.html>`_, `xugrid <https://deltares.github.io/xugrid/user_guide.html>`_, `pandas <https://pandas.pydata.org/docs/getting_started/index.html>`_, `geopandas <https://geopandas.org/en/stable/getting_started/introduction.html>`_ and `HYDROLIB-core <https://deltares.github.io/HYDROLIB-core>`_. Checking the `type()` of a variable will help to find the relevant documentation.

For pre-processing (e.g. the modelbuilder), the focus is often on data conversion. Raw data is being read with pandas, geopdandas or xarray, then processed (e.g. interpolated) and then written to model input files with xarray or HYDROLIB-core. For mesh generation, the `MeshKernelPy <https://deltares.github.io/MeshKernelPy/examples/index.html>`_ package is used.

For post-processing, the main sources of information are xarray and xugrid. For instance, the function `dfmt.open_partitioned_dataset()` returns a variable of type xugrid.UgridDataset. The `xugrid user guide <https://deltares.github.io/xugrid/user_guide.html>`_ contains many useful examples. The xugrid package in its turn wraps the underlying data as xarray datasets. This is a powerful package for lazy loading netcdf data (among others) and performing delayed operations on them. If you are unfamiliar with it, please read the `xarray in 45 minutes tutorial <https://tutorial.xarray.dev/overview/xarray-in-45-min.html>`_.



.. toctree::
   :titlesonly:
   :hidden:
   :maxdepth: 2

   notebooks/postprocessing_example.ipynb
   notebooks/postprocessing_example_delft3d4.ipynb
   notebooks/preprocessing_example_hydrolib.ipynb
   notebooks/modelbuilder_example.ipynb
   notebooks/subset_retrieve_sealevel_observations.ipynb