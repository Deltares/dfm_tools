=======
History
=======
0.7.33 (2021-03-19)
------------------
* made ghostid removal more efficient, partitioned map retrieval is now significantly faster
* including repair in mask of data, this is now repaired and a unittest is added to catch this in the future

0.7.31 (2021-01-13)
------------------
* added contextily basemap installation procedure and to testbank, more efficient than cartopy basemap

0.7.30 (2021-01-12)
------------------
* added fix to avoid crash with cartopy bug

0.7.29 (2021-01-12)
------------------
* fixed pandas empty columns issue in get_ncvardimlist occurring since 1.2.0 (Dec 26)
* several cleanup actions in tekal, waqua netcdf, github files
* updated html documentation
* added some getdata waqua commands to comments.

0.7.26 (2020-08-20)
------------------
* implemented first version of a zt-plot function for data from D-Flow FM hisfile

0.7.24 (2020-07-13)
------------------
* bugfix in var_times property of retrieved data, when retrieving as daterange
* bugfix in scatter_to_regulargrid(), masked values were not replaced by nans
* removed some non-unicode characters

0.7.23 (2020-06-12)
------------------
* fixed bug with retrieving non-partitioned variables from all partitions
* added ncdump function
* fixed bug with map merging file search
* added html documentation from docstrings (not all filled in yet)
* added read/write functions for bc-files

0.7.20 (2020-05-06)
------------------
* added option to retrieve data from top/bottom layers of z-layer D-Flow FM model
* fixed bug with empty string as varname
* added merge netcdf time function
* fixed dependencies (now all conda packages come from conda-forge channel)

0.7.19 (2020-04-28)
------------------
* improved time variable reading (more efficient when retrieving only a time-subset of a variable from a netCDF with long time dimension)
* improved time variable/dimensions reading (arbitrary time variable/dimension names are suported, as well as multiple time variables/dimensions)
* times were previously recalculated to UTC/GMT, this is now fixed
* conversion of negative indices to positive, sorting them and make unique
* read/write noos (matroos) data
* write bc file for D-Flow FM
* added example code to export D-Flow FM results to shapefile
* retrieving varname was possible from variable keys, now also possible from long_name or standard_name

0.7.6 (2020-04-06)
------------------
* Simplified installation method (check readme on github, link below)
* Improved retrieval on index (eg first and last timestep with [0,-1])
* Improved insights in variable contents/dimensions/shapes, to make it easier to know how to plot what with what
* Added regular grid features (reading eg wave grids and meteo data grids, meshgrid from xy vectors, corner2center, center2corner, corner2bounds, some plotting)
* Convert regular grid data to polycollection (same as ugrid.verts), so slicing (side view through 3D data) of regular grid is almost possible (this is still under construction)
* Read SFINCS map and his files
* Read virtually any NetCDF (ERA5, hirlam or other meteo files)
* Read Delft3D output (if this is in NetCDF output, you can get this by adding two keywords the .mdf)
* Read converted WAQUA/TRIWAQ output (converted to NetCDF with getdata.pl on h6, which works really well, let me know if you need help with this)
* Testbank now contains some new plot features like quivers, curved quivers and streamlines
* Plotting basemaps with cartopy land/ocean/landboundary/countryborders and a basic backgroundmap (proper satellite images still to be added)

0.6.4 (2020-03-19)
------------------
* Slightly different syntax which is better understandable (updated a while ago, so you probably will not notice, but it might be that you have to update your script)
* A first version of zt-plots (for instance the development over time of salinity of a station over the entire waterdepth)
* Matching function for WAQ statistics variables
* More flexible dimension reading (so more variables can be read)
* More robust
* Added sobek observation reading functionality (also netcdf)
* Added Delft3D grid and dep reading functionality (copied from OET)
* Added tekal reading functionality (for tek, pli, pliz, pol and ldb files)

0.2.0 (2020-02-14)
------------------
* restructured scripts and functions
* added safeguard for shapely import bug

0.1.16 (2020-02-14)
------------------
* correction for test case

0.1.15 (2020-02-14)
------------------
* test bank now properly coupled
* found solution for shapely bug (fix is in readme)

0.1.14 (2020-02-13)
------------------
* made intersect function more robust with exception cases
* added possibility to make cross section of 2D variable (was only available for 3D)

0.1.13 (2020-02-13)
------------------
* increased performance of grid/line intersection function (only check for intersections within lineboundbox)
* optimized intersect performance, added mercator if latlon

0.1.11 (2020-02-12)
------------------
* improved distance calculation in get_modeldata_onintersection function (second intersection function)

0.1.10 (2020-02-12)
------------------
* final hisfile-station fixes and updated tests script

0.1.9 (2020-02-12)
------------------
* added station selection for hisfiles, including updated testcases
* improved stability of layer retrieval

0.1.8 (2020-02-12)
------------------
* implemented first version of grid/line intersection function
* improved hisfile reading and made netfile reading more robust

0.1.7 (2020-02-11)
------------------
* added retrieval by datetime
* worked on his support
* made domain check more robust

0.1.5 (2020-02-10)
------------------
* improved his reading and dimension handling, updated testscript
* added checks for time/layer selection, made more robust
* added checks for timesteps and layers, also all times are possible
* added his and all times functionality
* fixed some bugs and made code neater and more efficient
* fixed indexing bug that surfaced with RMM model data
* plotting grids and mapvalues is now possible on certain depths and certain times, still very buggy and a lot left to do
* added plotting options for grids, including values as colors, but no multidomain yet
* added some tests, fixed grid.py to work with older variable names by adding translate function
* add tests and OET useful files

0.1.3 (2020-02-04)
------------------
* updated requirements.txt, less elaborate

0.1.2 (2020-02-04)
------------------
* transfered dflowutil to dfm_tools. write to any destination in utils.py
* dflowutils: allows writing to any destination, not just p drive

0.1.0 (2020-01-29)
------------------
* creation of the repository dfm_tools
