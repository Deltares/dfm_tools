# What's new

## 0.44.0 (2025-12-22)

# Feat
- update to new RWS Waterwebservices (and ddlpy) in [#1301](https://github.com/Deltares/dfm_tools/pull/1301)


## 0.43.0 (2025-11-27)

### Fix
- robustly retrieving bbox in `dfmt.meshkernel_get_bbox()` in [#1292](https://github.com/Deltares/dfm_tools/pull/1292)
- simpler parsing of uhslc geojson after pyogrio update in [#1294](https://github.com/Deltares/dfm_tools/pull/1294)
- update for `dfmt.meshkernel_get_illegalcells()` to align with new meshkernel version in [#909](https://github.com/Deltares/dfm_tools/issues/909)
- only include illegalcells up to 6 edges in `dfmt.meshkernel_get_illegalcells()` in [#1296](https://github.com/Deltares/dfm_tools/pull/1296)


## 0.42.0 (2025-11-03)

### Fix
- removed references to soon to be discontinued CMEMS multiyear-interim datasets in [#1279](https://github.com/Deltares/dfm_tools/pull/1279)
- allow for polyline names starting with numeric characters in [#1281](https://github.com/Deltares/dfm_tools/pull/1281)
- avoid writing of mdu keywords that are obsolete in Delft3D FM 2026.01 in [#1283](https://github.com/Deltares/dfm_tools/issues/1283)

### Deprecated
- removed poorly written and hycom-specific `dfmt.download_OPeNDAP()` in [#1289](https://github.com/Deltares/dfm_tools/pull/1289)


## 0.40.0 (2025-10-08)

### Feat
- support for list of files in `dfmt.cmems_nc_to_ini()` for better performance in [#1258](https://github.com/Deltares/dfm_tools/pull/1258)
- added polygon argument to `dfmt.refine_basegrid()` in [#1262](https://github.com/Deltares/dfm_tools/pull/1262)
- performance improvement for `dfmt.get_Dataset_atdepths()` in [#1264](https://github.com/Deltares/dfm_tools/pull/1264)
- support for multiple 2D quantities in `dfmt.Dataset_to_TimeSeries()` and therefore in `dfmt.plipointsDataset_to_ForcingModel()` in [#1268](https://github.com/Deltares/dfm_tools/pull/1268)


## 0.39.0 (2025-08-21)

### Feat
- extrapolate depth first in `dfmt.cmems_nc_to_ini()` in [#1231](https://github.com/Deltares/dfm_tools/pull/1231)
- avoid longpaths error in VSCode by removing direct fiona import in [#1239](https://github.com/Deltares/dfm_tools/pull/1239)
- make all paths in dfm_tools settable in [#1235](https://github.com/Deltares/dfm_tools/pull/1235)
- reduce filesize of bc files by rounding data to six decimal places in [#1240](https://github.com/Deltares/dfm_tools/pull/1240)


## 0.38.0 (2025-06-26)

### Fix
- made SSC linked stations in `ssc_add_linked_stations()` case-insensitive in [#6b3949f](https://github.com/Deltares/dfm_tools/commit/6b3949fdac152016b4c5f043f514ffd72c9f5e6e)
- avoid usage of moved private hydrolib-core function in [#1210](https://github.com/Deltares/dfm_tools/pull/1210)
- fixed z-sigma layer reconstruction for waterlevels below the z-sigma-interface in [#1219](https://github.com/Deltares/dfm_tools/pull/1219)
- improved aspect ratio in `dfmt.generate_basegrid()` by updating meshkernel dependency in [#1223](https://github.com/Deltares/dfm_tools/pull/1223)

### Feat
- enable separate retrieval of UHSLC rqds/fast observationdata again in `dfmt.ssh_retrieve_data()` (via `uhslc_ssh_retrieve_data()`) in [#1205](https://github.com/Deltares/dfm_tools/pull/1205)
- separate ERA5 quantities (supported from Delft3D-FM 2024.02) in [#1211](https://github.com/Deltares/dfm_tools/pull/1211)
- updated docker run script to efficiently work with Delft3D-FM 2025.02 (and later) docker containers in [#1216](https://github.com/Deltares/dfm_tools/pull/1216)


## 0.37.0 (2025-05-14)

This release drops support for Python 3.9.

### Fix
- retain encoding and long_name attribute also for variables converted with `convert_meteo_units()` in `dfmt.merge_meteofiles()` in [#1164](https://github.com/Deltares/dfm_tools/pull/1164)
- updated GESLA3 datasource in [#1173](https://github.com/Deltares/dfm_tools/pull/1173)
- stricter merging of datasets in `dfmt.cmems_nc_to_ini()` in [#1174](https://github.com/Deltares/dfm_tools/pull/1174)
- support for copernicusmarine 2.1.0 in [#1180](https://github.com/Deltares/dfm_tools/pull/1180)
- also construct sigmalayers when sigma-variables are coordinates [#1183](https://github.com/Deltares/dfm_tools/pull/1183)
- improved performance of `dfmt.uda_to_faces()` by using xugrid alternative in [#1177](https://github.com/Deltares/dfm_tools/pull/1177)
- check if all requested variables are present in merged dataset in `dfmt.preprocess_merge_meteofiles_era5()` in [#1184](https://github.com/Deltares/dfm_tools/pull/1184)
- remove accents from UHSLC-fast insitu catalog so `S64` datatype can be applied in `_make_hydrotools_consistent()` in [#1186](https://github.com/Deltares/dfm_tools/pull/1186)

### Deprecated
- removed `dfmt.preprocess_woa` since WOA merging fails and is not used in [#1161](https://github.com/Deltares/dfm_tools/pull/1161)
- support for python 3.9 is dropped in [#1177](https://github.com/Deltares/dfm_tools/pull/1177)
- deprecated sources "uhslc-fast" and "uhslc-rqds" in `dfmt.ssh_catalog_subset()` and `dfmt.ssh_retrieve_data()` in favor of "uhslc" in [#1193](https://github.com/Deltares/dfm_tools/pull/1193)

### Feat
- added cross-referencing SSC/IOC/UHSLC stations and distances in `ssc_ssh_read_catalog()` in [#1191](https://github.com/Deltares/dfm_tools/pull/1191)
- set coords and reduce filesize of retrieved observations in `dfmt.ssh_retrieve_data()` in [#1189](https://github.com/Deltares/dfm_tools/pull/1189)
- auto-merging of uhslc fast and rqds datasets in `dfmt.ssh_retrieve_data()` in [#1193](https://github.com/Deltares/dfm_tools/pull/1193)


## 0.36.0 (2025-03-18)

### Feat
- improved interpolation and extrapolation in `dfmt.cmems_nc_to_ini()` in [#1152](https://github.com/Deltares/dfm_tools/pull/1152)
- added `gtsm3-era5-cds` data as a source of observations in `dfmt.ssh_catalog_subset()` and `dfmt.ssh_retrieve_data()` in [#1153](https://github.com/Deltares/dfm_tools/pull/1153).

### Fix
- improved performance of `dfmt.merge_meteofiles()` by adding xarray arguments and ensuring alignment in [#1148](https://github.com/Deltares/dfm_tools/pull/1148)
- allow for number coordinate variable in `dfmt.merge_meteofiles()` in [#1157](https://github.com/Deltares/dfm_tools/pull/1157)
- proper time slicing and out of bounds error of `dfmt.merge_meteofiles()` in [#1149](https://github.com/Deltares/dfm_tools/pull/1149)


## 0.35.0 (2025-02-20)

### Feat
- add constant waterlevel offset with `dfmt.constant_to_bc()` in [#1130](https://github.com/Deltares/dfm_tools/pull/1130)
- download GSHHS data from github instead in [#1132](https://github.com/Deltares/dfm_tools/pull/1132)
- support multiple grid refinements via `kwargs` in`dfmt.refine_basegrid()` in [#1136](https://github.com/Deltares/dfm_tools/pull/1136)

### Fix
- correct ssr conversion factor in `convert_meteo_units()` in [#1134](https://github.com/Deltares/dfm_tools/pull/1134)
- support for `pathlib.Path` in `file_to_list()` in [#1139](https://github.com/Deltares/dfm_tools/pull/1139)
- rename new (feb 2025) CDS varnames `avg_ie` and `avg_tprate` back to `mer` and `mtpr` in `dfmt.preprocess_ERA5()` in [#1141](https://github.com/Deltares/dfm_tools/pull/1141)


## 0.34.0 (2025-02-05)

### Feat
- usage of outside time buffer in `dfmt.cmems_nc_to_ini()` so noon-centered or monthly timestamps are also supported in [#1087](https://github.com/Deltares/dfm_tools/pull/1087)
- correct CMEMS daily mean data ("P1D-m") from midnight to noon by adding a 12-hour offset in `dfmt.download_CMEMS()` in [#1088](https://github.com/Deltares/dfm_tools/pull/1088)
- updated cdsapi request to new format in [#1103](https://github.com/Deltares/dfm_tools/pull/1103)
- added UHSLC backup for GSHHS data in [#1112](https://github.com/Deltares/dfm_tools/pull/1112)
- auto download CMEMS phyc reanalysis from monthly mean dataset and auto convert to freq=M in [#622](https://github.com/Deltares/dfm_tools/issues/622)

### Fix
- made p-drive paths for tide models and gesla3 work on linux also in [#1083](https://github.com/Deltares/dfm_tools/pull/1083) and [#1085](https://github.com/Deltares/dfm_tools/pull/1085)
- support for different face dimension names in `dfmt.enrich_rst_with_map` in [#1114](https://github.com/Deltares/dfm_tools/pull/1114)


## 0.33.0 (2025-01-20)

### Feat
- optimized performance for getting CMEMS time extents and spatial buffer in [#1059](https://github.com/Deltares/dfm_tools/pull/1059)
- replaced buffer and floor/ceil with copernicusmarine `coordinates_selection_method`, this deprecated the `buffer` argument for `dfmt.download_CMEMS()` [#1061](https://github.com/Deltares/dfm_tools/pull/1061)
- inclusive selection of outside timesteps in `open_prepare_dataset()` and thus in `cmems_nc_to_bc()` in [#1062](https://github.com/Deltares/dfm_tools/pull/1062)

### Fix
- fixed inexact latlon bbox in modelbuilder with `dfmt.meshkernel_get_bbox()` in [#1067](https://github.com/Deltares/dfm_tools/pull/1067)
- included polyfile basename in bc-files from tide and cmems to prevent overwrite with multiple polyfiles in [#1071](https://github.com/Deltares/dfm_tools/pull/1071)
- renamed `ext_bnd` argument to `ext_new` in `dfmt.cmems_nc_to_bc()` for consistency in [#1071](https://github.com/Deltares/dfm_tools/pull/1071)


## 0.32.0 (2025-01-14)

### Feat
- updated to copernicusmarine v2 in [#1046](https://github.com/Deltares/dfm_tools/pull/1046)
- optimized CMEMS download performance in [#1049](https://github.com/Deltares/dfm_tools/pull/1049)


## 0.31.0 (2024-10-28)

### Fix
- fixed CDS login via prompt for new users in `cds_credentials()` in [#1035](https://github.com/Deltares/dfm_tools/pull/1035)


## 0.30.0 (2024-10-20)

### Fix
- fixed cmems-nrt insitu again by dropping station with varying coordinates again in [#1023](https://github.com/Deltares/dfm_tools/pull/1023)
- update executable paths in batfile created by `dfmt.create_model_exec_files()` in [#1030](https://github.com/Deltares/dfm_tools/pull/1030)


## 0.29.0 (2024-09-27)

### Feat
- update from CDS-Beta to new CDS in [#1013](https://github.com/Deltares/dfm_tools/pull/1013)


## 0.28.0 (2024-09-26)

### Fix
- more robust CDS/ECMWF authentication check and support for new cdsapi versions in [#1004](https://github.com/Deltares/dfm_tools/pull/1004)
- set default copernicusmarine buffer to 0.5 instead of 0 in [#1009](https://github.com/Deltares/dfm_tools/pull/1009)


## 0.27.0 (2024-09-09)

### Fix
- simplified `dfmt.meshkernel_to_UgridDataset()` by using new xugrid version in [#991](https://github.com/Deltares/dfm_tools/pull/991)
- updated `dfmt.uda_to_faces()` to accomodate for new fill_value of xugrid (this version includes fix for contour/contourf) [#992](https://github.com/Deltares/dfm_tools/pull/992)


## 0.26.0 (2024-09-03)

### Fix
- properly assigning units attribute in `ds_apply_conversion_dict()` (it did not always stick) in [#965](https://github.com/Deltares/dfm_tools/pull/965)
- skipping initialwaterlevel in `dfmt.cmems_nc_to_ini()` [#970](https://github.com/Deltares/dfm_tools/pull/970)
- update to cdsapi 0.7.2 and properly catching error for dummy dataset in [#972](https://github.com/Deltares/dfm_tools/pull/972)
- deprecated `dfmt.open_dataset_extra()` (partly replaced by `dfmt.open_prepare_dataset()`) in [#974](https://github.com/Deltares/dfm_tools/pull/974)
- improved nan-conversion in `dfmt.forcinglike_to_Dataset()` in [#982](https://github.com/Deltares/dfm_tools/pull/982)
- improved performance of `dfmt.open_partitioned_dataset()` for datasets with many variables in [#984](https://github.com/Deltares/dfm_tools/pull/984)


## 0.25.0 (2024-08-16)

### Feat
- making `dfmt.open_dataset_extra()` more modular by partly moving it to separate private functions in [#913](https://github.com/Deltares/dfm_tools/pull/913)
- added station_id variable to dataset returned by `dfmt.interp_uds_to_plipoints()` in [#914](https://github.com/Deltares/dfm_tools/pull/914)
- update private functions under `dfmt.download_ERA5()` to CDS-Beta (requires ECMWF apikey instead) in [#925](https://github.com/Deltares/dfm_tools/pull/925)
- simplified prevention of int dtypes in `dfmt.preprocess_ERA5()` in [#943](https://github.com/Deltares/dfm_tools/pull/943)
- simplified `dfmt.open_dataset_extra()` by dropping multi-quantity support and converted to private function in [#946](https://github.com/Deltares/dfm_tools/pull/946)
- improved `dfmt.interp_uds_to_plipoints()` by supporting outofbound points and new xugrid version in [#948](https://github.com/Deltares/dfm_tools/pull/948)
- neater chunking warning filter in `dfmt.open_partitioned_dataset()` in [#952](https://github.com/Deltares/dfm_tools/pull/952)
- performance improvement of `dfmt.open_partitioned_dataset()` by setting `remove_edges=False` in [#960](https://github.com/Deltares/dfm_tools/pull/960)

### Fix
- also apply convert_360to180 to longitude variable in `dfmt.open_dataset_curvilinear()` in [#913](https://github.com/Deltares/dfm_tools/pull/913)
- prevented adding of time dimension and dropped 0-sized cells in `dfmt.open_dataset_curvilinear()` in [#929](https://github.com/Deltares/dfm_tools/pull/929)
- dropping lat/lon/verts variables to maintain single source of truth in `dfmt.open_dataset_curvilinear()` in [#931](https://github.com/Deltares/dfm_tools/pull/931)
- fix for 3D initial fields in `dfmt.cmems_nc_to_ini()` to avoid top layer values over the entire depth in [#955](https://github.com/Deltares/dfm_tools/pull/955)


## 0.24.0 (2024-07-12)

### Feat
- improved usability of `dfmt.LineBuilder()` in [#854](https://github.com/Deltares/dfm_tools/pull/854)
- added workaround for grids that are not orthogonal after cutting the land with `dfmt.meshkernel_get_illegalcells()` in [#866](https://github.com/Deltares/dfm_tools/pull/866)
- updated CMEMS bcg multiyear dataset name in [#880](https://github.com/Deltares/dfm_tools/pull/880)
- added CMEMS reananalysis-interim (myint) datasets to `dfmt.download_CMEMS()` in [#883](https://github.com/Deltares/dfm_tools/pull/883) and [#903](https://github.com/Deltares/dfm_tools/pull/903)
- avoid duplicate and empty polyline names in `dfmt.geodataframe_to_PolyFile()` in [#896](https://github.com/Deltares/dfm_tools/pull/896)
- support for multiple polylines per polyfile in `interpolate_tide_to_bc()` and `cmems_nc_to_bc()` in [#906](https://github.com/Deltares/dfm_tools/pull/906)
- add to ext as part of `interpolate_tide_to_bc()`, it now requires `ext_new` and has no return value anymore in [#906](https://github.com/Deltares/dfm_tools/pull/906)

### Fix
- cleanups for datasets retrieved with `dfmt.ssh_retrieve_data()` in [#867](https://github.com/Deltares/dfm_tools/pull/867)


## 0.23.0 (2024-05-15)

### Feat
- added `conversion_dict` parameter to `dfmt.cmems_nc_to_bc()` in [#841](https://github.com/Deltares/dfm_tools/pull/841)

### Fix
- improved `dfmt.download.copernicusmarine_reset()` by changing order of tasks in [#833](https://github.com/Deltares/dfm_tools/pull/833)
- improved performance of `dfmt.merge_meteofiles()` by using more generic chunking method in [#840](https://github.com/Deltares/dfm_tools/pull/840)
- improved `dfmt.open_dataset_delft3d4()` by dropping original `grid` variable and attrs in [#844](https://github.com/Deltares/dfm_tools/pull/844)


## 0.22.0 (2024-04-09)

### Feat
- generate kml file from ssh_catalog (`geopandas.GeoDataFrame`) with `dfmt.ssh_catalog_tokmlfile()` in [#809](https://github.com/Deltares/dfm_tools/pull/809)
- generate time availability overview and statistics from insitu netcdf files with `dfmt.ssh_netcdf_overview()` in [#807](https://github.com/Deltares/dfm_tools/pull/807)
- simplified `rwsddl_ssh_retrieve_data()` by using `ddlpy.dataframe_to_xarray()` in [#805](https://github.com/Deltares/dfm_tools/pull/815)
- updated cmems analysisforecast bio dataset ids in [#816](https://github.com/Deltares/dfm_tools/pull/816)
- from ftp to copernicusmarine files service for insitu data in [#818](https://github.com/Deltares/dfm_tools/pull/818)
- simplified `dfmt.download.copernicusmarine_credentials()` in [#819](https://github.com/Deltares/dfm_tools/pull/819)

### Fix
- fixed xy coordinates for cmems insitu netcdf files in [#807](https://github.com/Deltares/dfm_tools/pull/807)


## 0.21.0 (2024-03-14)

### Feat
- added ddl/rws insitu data to `dfmt.ssh_catalog_subset()` and `dfmt.ssh_retrieve_data()`  in [#791](https://github.com/Deltares/dfm_tools/pull/791) and [#796](https://github.com/Deltares/dfm_tools/pull/796)
- added cmems nrt insitu data to `dfmt.ssh_catalog_subset()` and `dfmt.ssh_retrieve_data()` in [#791](https://github.com/Deltares/dfm_tools/pull/791) and [#793](https://github.com/Deltares/dfm_tools/pull/793)
- support for meshkernel>=4.1.0 in [#801](https://github.com/Deltares/dfm_tools/pull/801)


## 0.20.0 (2024-02-22)

### Fix
- solve xarray FutureWarning about `ds.dims[dim]` by [@veenstrajelmer](https://github.com/veenstrajelmer) in [#779](https://github.com/Deltares/dfm_tools/pull/779)
- support copernicusmarine>=1.0.2 in `copernicusmarine_credentials()` and made `copernicusmarine_reset()` modular by [@veenstrajelmer](https://github.com/veenstrajelmer) in [#787](https://github.com/Deltares/dfm_tools/pull/787)


## 0.19.0 (2024-02-08)

### Feat
- automatically parse epsg from FM mapfile to crs in `dfmt.open_partitioned_dataset()` by [@veenstrajelmer](https://github.com/veenstrajelmer) in [#685](https://github.com/Deltares/dfm_tools/pull/685)
- added spatial/temporal subsetting and retrieving of insitu observation data with `dfmt.ssh_catalog_subset()` and `dfmt.ssh_retrieve_data()` by [@veenstrajelmer](https://github.com/veenstrajelmer) in [#711](https://github.com/Deltares/dfm_tools/pull/711)
- added conversion of absolute to relative paths for modelbuilder with `dfmt.make_paths_relative()` by [@veenstrajelmer](https://github.com/veenstrajelmer) in [#723](https://github.com/Deltares/dfm_tools/pull/723)
- add retrieval and plotting of borders with `dfmt.get_borders_gdb()` and `dfmt.plot_borders()` by [@veenstrajelmer](https://github.com/veenstrajelmer) in [#750](https://github.com/Deltares/dfm_tools/pull/750)
- improve z and sigmaz layer reconstruction (nans below bed and above waterlevel and centers as mean of interfaces) in `reconstruct_zw_zcc()` by [@veenstrajelmer](https://github.com/veenstrajelmer) in [#755](https://github.com/Deltares/dfm_tools/pull/755) and [#759](https://github.com/Deltares/dfm_tools/pull/759)
- simplified correctly merging of ERA5 int16 meteofiles in `merge_meteofiles()` by [@veenstrajelmer](https://github.com/veenstrajelmer) in [#738](https://github.com/Deltares/dfm_tools/pull/738)
- expecting da instead of ds in `refine_basegrid()` by [@veenstrajelmer](https://github.com/veenstrajelmer) in [#762](https://github.com/Deltares/dfm_tools/pull/762)
- added support for Python 3.12 by [@veenstrajelmer](https://github.com/veenstrajelmer) in [#763](https://github.com/Deltares/dfm_tools/pull/763)
- moved from `copernicus-marine-client` to `copernicusmarine` in `dfmt.download_CMEMS()` by [@veenstrajelmer](https://github.com/veenstrajelmer) in [#721](https://github.com/Deltares/dfm_tools/pull/721) and [#772](https://github.com/Deltares/dfm_tools/pull/772)
- added reset of Copernicus Marine Toolbox with `dfmt.download.copernicusmarine_reset()` by [@veenstrajelmer](https://github.com/veenstrajelmer) in [#714](https://github.com/Deltares/dfm_tools/pull/714)

### Fix
- prevent ValueError upon concatenation of only emtpy coastlines geodataframes by [@veenstrajelmer](https://github.com/veenstrajelmer) in [#717](https://github.com/Deltares/dfm_tools/pull/717)
- added dropped crs variable in network file by [@veenstrajelmer](https://github.com/veenstrajelmer) in [#719](https://github.com/Deltares/dfm_tools/pull/719)
- moved gtsmv4.1 tide interpolation from extrapolated to rasterized dataset by [@veenstrajelmer](https://github.com/veenstrajelmer) in [#731](https://github.com/Deltares/dfm_tools/pull/731)
- limited `dfmt.download_CMEMS()` freq argument to M/D only by [@veenstrajelmer](https://github.com/veenstrajelmer) in [#733](https://github.com/Deltares/dfm_tools/pull/733)
- suppres xarray chunking warning by [@veenstrajelmer](https://github.com/veenstrajelmer) in [#742](https://github.com/Deltares/dfm_tools/pull/742)
- remove mesh2d hardcoding from layer reconstruction by [@veenstrajelmer](https://github.com/veenstrajelmer) in [#747](https://github.com/Deltares/dfm_tools/pull/747)
- generalized finding of sigma layer/interface variables via attrs in `reconstruct_zw_zcc()` by [@veenstrajelmer](https://github.com/veenstrajelmer) in [#765](https://github.com/Deltares/dfm_tools/pull/765)

### Deprecated
- deprecated `dfmt.preprocess_hirlam()` since xarray>=2023.9.0 supports multidimensional coordinates by [@veenstrajelmer](https://github.com/veenstrajelmer) in [#770](https://github.com/Deltares/dfm_tools/pull/770)


## 0.18.0 (2023-12-08)

### Feat
- replaced `is_geographic` with `crs` argument in `dfmt.make_basegrid()` by [@veenstrajelmer](https://github.com/veenstrajelmer) in [#696](https://github.com/Deltares/dfm_tools/pull/696)
- moved from `OPeNDAP` to `copernicus-marine-client` in `dfmt.download_CMEMS()` including safer authentication by [@veenstrajelmer](https://github.com/veenstrajelmer) in [#691](https://github.com/Deltares/dfm_tools/pull/691)
- added timerange check of resulting dataset of `dfmt.open_dataset_extra()` by [@veenstrajelmer](https://github.com/veenstrajelmer) in [#708](https://github.com/Deltares/dfm_tools/pull/708)


## 0.17.0 (2023-11-17)
Support for Python 3.8 was dropped ([more info](https://github.com/Deltares/dfm_tools/issues/267))

### Feat
- pop np.nan `_FillValue` attrs in `dfmt.open_partitioned_dataset()` by [@veenstrajelmer](https://github.com/veenstrajelmer) in [#662](https://github.com/Deltares/dfm_tools/pull/662)
- interpolation of edge/node variables to faces with `dfmt.uda_to_faces()` (deprecates `dfmt.uda_edges_to_faces()`) by [@veenstrajelmer](https://github.com/veenstrajelmer) in [#651](https://github.com/Deltares/dfm_tools/pull/651) and [#644](https://github.com/Deltares/dfm_tools/pull/644)
- meshkernel 3.0.0 support by [@veenstrajelmer](https://github.com/veenstrajelmer) in [#621](https://github.com/Deltares/dfm_tools/pull/621) and [#665](https://github.com/Deltares/dfm_tools/pull/665)
- removed support for Python 3.8 by [@veenstrajelmer](https://github.com/veenstrajelmer) in [#666](https://github.com/Deltares/dfm_tools/pull/666)

### Fix
- get `is_geographic` from crs instead of hardcoded in `add_crs_to_dataset` by [@veenstrajelmer](https://github.com/veenstrajelmer) in [#657](https://github.com/Deltares/dfm_tools/pull/657)
- performance improvements of `dfmt.uda_to_faces()` by [@veenstrajelmer](https://github.com/veenstrajelmer) in [#652](https://github.com/Deltares/dfm_tools/pull/652)


## 0.16.0 (2023-11-03)

### Feat
- more robust support for CMEMS/CDS credentials including environment variables by [@veenstrajelmer](https://github.com/veenstrajelmer) in [#633](https://github.com/Deltares/dfm_tools/issues/633)
- enrich rst file with topology from corresponding mapfile with `dfmt.enrich_rst_with_map()` by [@veenstrajelmer](https://github.com/veenstrajelmer) in [#626](https://github.com/Deltares/dfm_tools/issues/626)
- add cellinfo to minimal 2D networks (with 1D topology) with `dfmt.add_network_cellinfo()` by [@veenstrajelmer](https://github.com/veenstrajelmer) in [#624](https://github.com/Deltares/dfm_tools/issues/624)
- xugrid feature `uds.ugrid.to_nonperiodic()` deprecates `dfmt.remove_periodic_cells()` by [@veenstrajelmer](https://github.com/veenstrajelmer) in [#503](https://github.com/Deltares/dfm_tools/issues/503)
- support for initial fields for variables other than salinity/temperature with `dfmt.cmems_nc_to_ini()` by [@veenstrajelmer](https://github.com/veenstrajelmer) in [#619](https://github.com/Deltares/dfm_tools/pull/619)

### Fix
- fix initial fields for salinity/temperature with `dfmt.cmems_nc_to_ini()` by [@veenstrajelmer](https://github.com/veenstrajelmer) in [#619](https://github.com/Deltares/dfm_tools/pull/619)
- increased buffer in `dfmt.download_ERA5()` by [@veenstrajelmer](https://github.com/veenstrajelmer) in [#612](https://github.com/Deltares/dfm_tools/pull/612)
- support for Polygon geometries in `dfmt.geodataframe_to_PolyFile()` by [@veenstrajelmer](https://github.com/veenstrajelmer) in [#610](https://github.com/Deltares/dfm_tools/pull/610)
- fill nan-values in initial salinity/temperature netcdf dataset in `dfmt.preprocess_ini_cmems_to_nc()` by [@veenstrajelmer](https://github.com/veenstrajelmer) in [#617](https://github.com/Deltares/dfm_tools/pull/617)
- skip all-nan boundary support points instead of converting to zeros in `dfmt.plipointsDataset_to_ForcingModel()` by [@veenstrajelmer](https://github.com/veenstrajelmer) in [#637](https://github.com/Deltares/dfm_tools/pull/637)


## 0.15.0 (2023-10-19)

### Feat
- create `dimr_config.xml` and `run_parallel.bat` with `dfmt.create_model_exec_files()` by [@veenstrajelmer](https://github.com/veenstrajelmer) in [#599](https://github.com/Deltares/dfm_tools/pull/599)
- conversion of `hcdfm.ForcingModel` to `xr.Dataset` with `dfmt.ForcingModel_to_plipointsDataset()` by [@veenstrajelmer](https://github.com/veenstrajelmer) in [#567](https://github.com/Deltares/dfm_tools/pull/567)
- support for reading asc files with `dfmt.read_asc()` by [@veenstrajelmer](https://github.com/veenstrajelmer) in [#571](https://github.com/Deltares/dfm_tools/pull/571)
- added `dfmt.meshkernel_delete_withshp()` to delete parts of a mesh with a shapefile, while only reading the shapefile within the bounding box of the mesh by [@rqwang](https://github.com/rqwang) in [#548](https://github.com/Deltares/dfm_tools/pull/548) and [#566](https://github.com/Deltares/dfm_tools/pull/566)
- improved spatial interpolation in `dfmt.interp_regularnc_to_plipoints()` by combining linear with nearest interpolation by [@veenstrajelmer](https://github.com/veenstrajelmer) in [#561](https://github.com/Deltares/dfm_tools/pull/561)
- added `GTSMv4.1` and `GTSMv4.1_opendap` datasets for tide interpolation with `dfmt.interpolate_tide_to_bc()` by [@veenstrajelmer](https://github.com/veenstrajelmer) in [#544](https://github.com/Deltares/dfm_tools/pull/544)
- support for `preprocess` argument in `dfmt.open_partitioned_dataset()` by [@veenstrajelmer](https://github.com/veenstrajelmer) in [#530](https://github.com/Deltares/dfm_tools/pull/530)

### Fix
- increased buffer for WAQ cmems variables by [@veenstrajelmer](https://github.com/veenstrajelmer) in [#601](https://github.com/Deltares/dfm_tools/pull/601)
- prevent concatenation of datasets with slightly different coordinates by [@veenstrajelmer](https://github.com/veenstrajelmer) in [#595](https://github.com/Deltares/dfm_tools/pull/595)
- aligned ncbnd dataset with FM/FEWS conventions by [@veenstrajelmer](https://github.com/veenstrajelmer) in [#567](https://github.com/Deltares/dfm_tools/pull/567)
- documented minimal dependency versions to avoid issues with installation in existing environments containing old dependency versions by [@veenstrajelmer](https://github.com/veenstrajelmer) in [#581](https://github.com/Deltares/dfm_tools/pull/581)


## 0.14.0 (2023-09-15)

### Feat
- decoding of default fillvalues with `dfmt.decode_default_fillvals()` by [@veenstrajelmer](https://github.com/veenstrajelmer) in [#496](https://github.com/Deltares/dfm_tools/pull/496)
- interpolation of edge variables to faces with `dfmt.uda_edges_to_faces()` by [@veenstrajelmer](https://github.com/veenstrajelmer) in [#482](https://github.com/Deltares/dfm_tools/pull/482)
- interpolation of variables on layer interfaces to layer centers with `dfmt.uda_interfaces_to_centers()` by [@veenstrajelmer](https://github.com/veenstrajelmer) in [#482](https://github.com/Deltares/dfm_tools/pull/482)
- generate polyline only for open boundary with `dfmt.generate_bndpli_cutland()` (deprecates `dfmt.generate_bndpli()`) by [@veenstrajelmer](https://github.com/veenstrajelmer) in [#466](https://github.com/Deltares/dfm_tools/pull/466)
- automatically remove unassociated edges with `dfmt.remove_unassociated_edges()` by [@veenstrajelmer](https://github.com/veenstrajelmer) in [#470](https://github.com/Deltares/dfm_tools/pull/470) and [#494](https://github.com/Deltares/dfm_tools/pull/494)
- clarifying error message when providing incorrectly formatted cds api key by [@JulienGroenenboom](https://github.com/JulienGroenenboom) in [#514](https://github.com/Deltares/dfm_tools/pull/514)

### Fix
- removed `matplotlib._api` import to also support `matplotlib<3.4.0` by [@veenstrajelmer](https://github.com/veenstrajelmer) in [#465](https://github.com/Deltares/dfm_tools/pull/465)


## 0.13.1 (2023-07-12)

### Fix
- proper fix for FM-compatible network file (1-based face_node_connectivity) by [@veenstrajelmer](https://github.com/veenstrajelmer) in [#461](https://github.com/Deltares/dfm_tools/pull/461)


## 0.13.0 (2023-07-12)

### Feat
- exposed opendap testdata with `dfmt.data.fm_grevelingen_map()` (and others) by [@veenstrajelmer](https://github.com/veenstrajelmer) in [#439](https://github.com/Deltares/dfm_tools/pull/439)
- made coastlines portable with ` dfmt.data.gshhs_coastlines_shp()` by [@veenstrajelmer](https://github.com/veenstrajelmer) in [#443](https://github.com/Deltares/dfm_tools/pull/443)
- improved CMEMS and ERA5 authentication via [getpass](https://docs.python.org/3/library/getpass.html) by [@veenstrajelmer](https://github.com/veenstrajelmer) in [#449](https://github.com/Deltares/dfm_tools/pull/449)
- aligned mesh deletion with `dfmt.meshkernel_delete_withgdf()` by [@veenstrajelmer](https://github.com/veenstrajelmer) in [#455](https://github.com/Deltares/dfm_tools/pull/455)

### Fix
- added `dfmt.uds_to_1based_ds()` to ensure FM-compatible network file (1-based face_node_connectivity) by [@veenstrajelmer](https://github.com/veenstrajelmer) in [#456](https://github.com/Deltares/dfm_tools/pull/456)


## 0.12.0 (2023-07-07)

### Feat
- added support for curvilinear datasets like CMCC and WAQUA with `open_dataset_curvilinear()` by [@veenstrajelmer](https://github.com/veenstrajelmer) in [#397](https://github.com/Deltares/dfm_tools/pull/397)
- added support for curvilinear Delft3D4 datasets with `open_dataset_delft3d4()` by [@veenstrajelmer](https://github.com/veenstrajelmer) in [#398](https://github.com/Deltares/dfm_tools/pull/398)
- retrieving [GSHHS landboundary](https://www.ngdc.noaa.gov/mgg/shorelines/) with `get_coastlines_gdb()` and cartopy-alternative `plot_coastlines()` with this dataset by [@veenstrajelmer](https://github.com/veenstrajelmer) in [#342](https://github.com/Deltares/dfm_tools/pull/342)
- added interpolation from UgridDataset (D-FlowFM mapfile) to plipoints with `interp_uds_to_plipoints()` by [@veenstrajelmer](https://github.com/veenstrajelmer) [#387](https://github.com/Deltares/dfm_tools/pull/387)
- added geodataframe support from polyfile with `PolyFile_to_geodataframe_points()` by [@veenstrajelmer](https://github.com/veenstrajelmer) in [#387](https://github.com/Deltares/dfm_tools/pull/387)

### Fix
- proper selection for either CMEMS reanalysis/forecast product in `download_CMEMS()` by [@JulienGroenenboom](https://github.com/JulienGroenenboom) in [#308](https://github.com/Deltares/dfm_tools/pull/388)
- slightly extend spatial and time domain when downloading ERA5 and CMEMS data by [@JulienGroenenboom](https://github.com/JulienGroenenboom) in [#390](https://github.com/Deltares/dfm_tools/pull/390)
- renamed uxuy to uxuyadvectionvelocitybnd in `open_dataset_extra()` and related example script by [@veenstrajelmer](https://github.com/veenstrajelmer) in [#366](https://github.com/Deltares/dfm_tools/pull/366)
- add warning when interpolating to a boundary that contains multiple polylines by [@Tammo-Zijlker-Deltares](https://github.com/Tammo-Zijlker-Deltares) in [#372](https://github.com/Deltares/dfm_tools/pull/372)
- improved netcdf writing from meshkernel+xugrid generated network with `add_crs_to_dataset()` by [@veenstrajelmer](https://github.com/veenstrajelmer) in [#422](https://github.com/Deltares/dfm_tools/pull/422)


## 0.11.0 (2023-04-24)

Replaced large parts of the code with [HYDROLIB-core](https://github.com/Deltares/hydrolib-core), [xugrid](https://github.com/Deltares/xugrid) and [xarray](https://github.com/pydata/xarray). This means a non-backwards compatible change in the API but improved experience. xugrid gives loads of new features like coordinate conversion, shapefile export, cropping a unstructured dataset, regridding, rasterizing.

### Feat
- improved support for ugrid datasets (D-FlowFM, D-HYDRO) with `open_partitioned_dataset()` by [@veenstrajelmer](https://github.com/veenstrajelmer) in several PRs, this deprecates `get_ncmodeldata()`, `get_netdata()`, `plot_netmapdata()` and `get_ugrid_verts()`
- improved rasterization of ugrid datasets with `rasterize_ugrid()` by [@veenstrajelmer](https://github.com/veenstrajelmer) in [#314](https://github.com/Deltares/dfm_tools/pull/314), this deprecates `scatter_to_regulargrid()`
- depth slicing of ugrid datasets with `get_Dataset_atdepths()` by [@veenstrajelmer](https://github.com/veenstrajelmer) in [#261](https://github.com/Deltares/dfm_tools/pull/261)
- polygon slicing of ugrid datasets with `polyline_mapslice()` by [@veenstrajelmer](https://github.com/veenstrajelmer) in [#266](https://github.com/Deltares/dfm_tools/pull/266) and [#273](https://github.com/Deltares/dfm_tools/pull/273)
- reconstruction of fullgrid ugrid variable from sigma, zsigma and z-layers with `reconstruct_zw_zcc()` by [@veenstrajelmer](https://github.com/veenstrajelmer) in [#291](https://github.com/Deltares/dfm_tools/pull/291) and [#293](https://github.com/Deltares/dfm_tools/pull/293)
- added support for D-FlowFM mapformat 1 by [@veenstrajelmer](https://github.com/veenstrajelmer) in [#249](https://github.com/Deltares/dfm_tools/pull/249)
- support for plotting of periodic grid with `remove_periodic_cells()` by [@veenstrajelmer](https://github.com/veenstrajelmer) in [#286](https://github.com/Deltares/dfm_tools/pull/286)
- renaming water quality variables in dataset with `rename_waqvars()` by [@veenstrajelmer](https://github.com/veenstrajelmer) in [#259](https://github.com/Deltares/dfm_tools/pull/259), this deprecates `get_varnamefromattrs()`
- renaming variables in fourier dataset with `rename_fouvars()` by [@veenstrajelmer](https://github.com/veenstrajelmer) in [#283](https://github.com/Deltares/dfm_tools/pull/283)
- simplified background plotting with contextily and cartopy by [@veenstrajelmer](https://github.com/veenstrajelmer) in [#328](https://github.com/Deltares/dfm_tools/pull/328), this deprecates `plot_background()`
- create release and pypi/pip installable by [@veenstrajelmer](https://github.com/veenstrajelmer) in [#338](https://github.com/Deltares/dfm_tools/pull/338)
- using [HYDROLIB-core](https://github.com/Deltares/hydrolib-core) for reading/writing of D-FlowFM files, this deprecates `write_bcfile()`, `read_bcfile()`, `Polygon()`, `write_timfile()` and `read_timfile()`
- interpolation of regularly gridded dataset (CMEMS, HYCOM etc) to plipoints with `interp_regularnc_to_plipoints()` by [@veenstrajelmer](https://github.com/veenstrajelmer)
- interpolation of regularly gridded tide datasets to plipoints with `interpolate_tide_to_bc()` by [@veenstrajelmer](https://github.com/veenstrajelmer) in [#288](https://github.com/Deltares/dfm_tools/pull/288)
- added modelbuilder by [@veenstrajelmer](https://github.com/veenstrajelmer) in [#305](https://github.com/Deltares/dfm_tools/pull/305)
- added meshkernelpy helper functions like `make_basegrid()`, `refine_basegrid()` and `meshkernel_to_UgridDataset()` by [@veenstrajelmer](https://github.com/veenstrajelmer) in [#329](https://github.com/Deltares/dfm_tools/pull/329)

### Fix
- made bedlevel variable optional in `get_Dataset_atdepths()`, if not explicitly needing it by [@veenstrajelmer](https://github.com/veenstrajelmer) in [#241](https://github.com/Deltares/dfm_tools/pull/241)
