## 0.12.0 (unreleased)

### Feat

- open_dataset_curvilinear() added support for curvilinear datasets like CMCC and WAQUA by @veenstrajelmer in https://github.com/Deltares/dfm_tools/pull/397
- open_dataset_delft3d4() added support for curvilinear Delft3D4 datasets by @veenstrajelmer in https://github.com/Deltares/dfm_tools/pull/398

### Fix

- proper selection for either CMEMS reanalysis/forecast product in download_CMEMS() by @JulienGroenenboom in https://github.com/Deltares/dfm_tools/pull/388

## 0.11.0 (2023-04-24)

### Feat

- dfmt.open_partitioned dataset() replaces get_ncmodeldata and related functions, now reading of ugrid (D-FlowFM, D-HYDRO) files is done with [xugrid](https://github.com/Deltares/xugrid) and [xarray](https://github.com/pydata/xarray).

- write to any destination in utils.py by @LBuckma in https://github.com/Deltares/dfm_tools/pull/2
- support for any number of headerlines by @LBuckma in https://github.com/Deltares/dfm_tools/pull/106
- prevented necessity for bedlevel in get_Dataset_atdepths(), if not explicitly needing it by @veenstrajelmer in https://github.com/Deltares/dfm_tools/pull/241
- removed additional grid plot method, is now in xugrid by @veenstrajelmer in https://github.com/Deltares/dfm_tools/pull/243
- 210 deprecate dfmtget ncmodeldata and related functions by @veenstrajelmer in https://github.com/Deltares/dfm_tools/pull/246
- 248 add support for mapformat=1 by @veenstrajelmer in https://github.com/Deltares/dfm_tools/pull/249
- 250 deprecate dfmtget ugrid verts by @veenstrajelmer in https://github.com/Deltares/dfm_tools/pull/251
- temporarily dropping edge/node/interface dims by @veenstrajelmer in https://github.com/Deltares/dfm_tools/pull/255
- 254 temporarily drop edgenodeinterface dims to avoid extra dimensions in returned dataset by @veenstrajelmer in https://github.com/Deltares/dfm_tools/pull/257
- 258 add ds rename funtion for waq variables by @veenstrajelmer in https://github.com/Deltares/dfm_tools/pull/259
- 208 improve dfmtget dataset atdepths by @veenstrajelmer in https://github.com/Deltares/dfm_tools/pull/261
- 262 only return depth sliced variables in get dataset atdepths by @veenstrajelmer in https://github.com/Deltares/dfm_tools/pull/263
- 252 convert map cross section polygon slice to xugrid alternative by @veenstrajelmer in https://github.com/Deltares/dfm_tools/pull/266
- 239 xarray writing mfdataset results in incorrect data when not using manual encoding by @veenstrajelmer in https://github.com/Deltares/dfm_tools/pull/268
- 269 recompute scalingoffset for dtypeint mfdataset ds by @veenstrajelmer in https://github.com/Deltares/dfm_tools/pull/270
- 271 cleanup dfmtregulargridrasterize ugrid input arguments by @veenstrajelmer in https://github.com/Deltares/dfm_tools/pull/272
- 252 convert map cross section polygon slice to xugrid alternative 2 by @veenstrajelmer in https://github.com/Deltares/dfm_tools/pull/273
- added test_intersect_edges() and test_intersect_edges_withsort() by @veenstrajelmer in https://github.com/Deltares/dfm_tools/pull/276
- cleaned up deprecated code including testcase and added new unittests by @veenstrajelmer in https://github.com/Deltares/dfm_tools/pull/278
- 279 remove deprecated testcases by @veenstrajelmer in https://github.com/Deltares/dfm_tools/pull/280
- 281 remove test import libraries by @veenstrajelmer in https://github.com/Deltares/dfm_tools/pull/282
- added rename_fouvars() function including testcase, removed deprecate… by @veenstrajelmer in https://github.com/Deltares/dfm_tools/pull/283
- added remove_periodic_cells() function for GTSM "around the back" cells by @veenstrajelmer in https://github.com/Deltares/dfm_tools/pull/286
- 285 add gtsm to interpolate tide to bc by @veenstrajelmer in https://github.com/Deltares/dfm_tools/pull/288
- 289 implement reconstruction of fullgrid depths with new depth variables by @veenstrajelmer in https://github.com/Deltares/dfm_tools/pull/291
- 292 neater reconstruction of zsigma and sigma via formula terms by @veenstrajelmer in https://github.com/Deltares/dfm_tools/pull/293
- cleaned old polygon from code smells by @veenstrajelmer in https://github.com/Deltares/dfm_tools/pull/298
- 297 solve sonarcloud code smells by @veenstrajelmer in https://github.com/Deltares/dfm_tools/pull/299
- 297 solve sonarcloud code smells by @veenstrajelmer in https://github.com/Deltares/dfm_tools/pull/300
- cleaned up get_xzcoords_onintersection and facevar selection instead of edge/node dim dropping by @veenstrajelmer in https://github.com/Deltares/dfm_tools/pull/303
- 304 add modelbuilder to dfm tools by @veenstrajelmer in https://github.com/Deltares/dfm_tools/pull/305
- 301 convert timmodel to pandasdataframe by @veenstrajelmer in https://github.com/Deltares/dfm_tools/pull/306
- 296 also compute zcc coordinates in reconstruct zw zcc by @veenstrajelmer in https://github.com/Deltares/dfm_tools/pull/308
- 296 also compute zcc coordinates in reconstruct zw zcc by @veenstrajelmer in https://github.com/Deltares/dfm_tools/pull/309
- moved rasterize_ugrid location by @veenstrajelmer in https://github.com/Deltares/dfm_tools/pull/314
- 310 deprecate get varnamefromattrs including unittest by @veenstrajelmer in https://github.com/Deltares/dfm_tools/pull/315
- updated minimal versions for meshkernel and hydrolib-core by @veenstrajelmer in https://github.com/Deltares/dfm_tools/pull/316
- 321 update cmems download urls and make more user friendly by @veenstrajelmer in https://github.com/Deltares/dfm_tools/pull/322
- 324 test pip install by @veenstrajelmer in https://github.com/Deltares/dfm_tools/pull/325
- 327 remove plot background cartopy function by @veenstrajelmer in https://github.com/Deltares/dfm_tools/pull/328
- 323 move meshkernel functions from modelbuilder to meshkernel helpers by @veenstrajelmer in https://github.com/Deltares/dfm_tools/pull/329
- fixed xarray version to ensure successful plotting of chunked ugrid datasets by @veenstrajelmer in https://github.com/Deltares/dfm_tools/pull/331
- 332 loosen fixed dependency versions after issues are solved by @veenstrajelmer in https://github.com/Deltares/dfm_tools/pull/333
- 334 move deprecated functions to deprecatedpy by @veenstrajelmer in https://github.com/Deltares/dfm_tools/pull/335
- 336 add pytest testbank for py38 py310 py311 by @veenstrajelmer in https://github.com/Deltares/dfm_tools/pull/337
- 326 create release and pypipip installable by @veenstrajelmer in https://github.com/Deltares/dfm_tools/pull/338
- updated CONTRIBUTING.md by @veenstrajelmer in https://github.com/Deltares/dfm_tools/pull/340
