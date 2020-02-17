=======
History
=======

0.4.1 (2020-02-17)
------------------

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
