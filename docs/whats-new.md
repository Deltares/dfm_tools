## 0.12.0 (unreleased)

### Feat

- added support for curvilinear datasets like CMCC and WAQUA with `open_dataset_curvilinear()` by [@veenstrajelmer](https://github.com/veenstrajelmer) in [#397](https://github.com/Deltares/dfm_tools/pull/397)
- added support for curvilinear Delft3D4 datasets with `open_dataset_delft3d4()` by [@veenstrajelmer](https://github.com/veenstrajelmer) in [#398](https://github.com/Deltares/dfm_tools/pull/398)
- retrieving [GSHHS landboundary](https://www.ngdc.noaa.gov/mgg/shorelines/) with `get_coastlines_gdb()` and cartopy-alternative `plot_coastlines()` with these dataset by [@veenstrajelmer](https://github.com/veenstrajelmer) in [#342](https://github.com/Deltares/dfm_tools/pull/342)
- added interpolation from UgridDataset (D-FlowFM mapfile) to plipoints with `interp_uds_to_plipoints()` by [@veenstrajelmer](https://github.com/veenstrajelmer) [#387](https://github.com/Deltares/dfm_tools/pull/387)
- added geodataframe support from polyfile with `PolyFile_to_geodataframe_points()` by [@veenstrajelmer](https://github.com/veenstrajelmer) in [#387](https://github.com/Deltares/dfm_tools/pull/387)

### Fix

- proper selection for either CMEMS reanalysis/forecast product in `download_CMEMS()` by [@JulienGroenenboom](https://github.com/JulienGroenenboom) in [#308](https://github.com/Deltares/dfm_tools/pull/388)
- renamed uxuy to uxuyadvectionvelocitybnd in `open_dataset_extra()` and related example script by [@veenstrajelmer](https://github.com/veenstrajelmer) in [#366](https://github.com/Deltares/dfm_tools/pull/366)
- add warning when interpolating to a boundary that contains multiple polylines by [@Tammo-Zijlker-Deltares](https://github.com/Tammo-Zijlker-Deltares)[#372](https://github.com/Deltares/dfm_tools/pull/372)

## 0.11.0 (2023-04-24)

Replaced large parts of the code with [HYDROLIB-core](https://github.com/Deltares/hydrolib-core), [xugrid](https://github.com/Deltares/xugrid) and [xarray](https://github.com/pydata/xarray). This means a non-backwards compatible change in the API but improved experience.

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
- create release and pypipip installable by [@veenstrajelmer](https://github.com/veenstrajelmer) in [#338](https://github.com/Deltares/dfm_tools/pull/338)
- using [HYDROLIB-core](https://github.com/Deltares/hydrolib-core) for reading/writing of D-FlowFM files, this deprecates `write_bcfile()`, `read_bcfile()`, `Polygon()`, `write_timfile()` and `read_timfile()`
- interpolation of regularly gridded dataset (CMEMS, HYCOM etc) to plipoints with `interp_regularnc_to_plipoints()` by [@veenstrajelmer](https://github.com/veenstrajelmer)
- interpolation of regularly gridded tide datasets to plipoints with `interpolate_tide_to_bc()` by [@veenstrajelmer](https://github.com/veenstrajelmer) in [#288](https://github.com/Deltares/dfm_tools/pull/288)
- added modelbuilder by [@veenstrajelmer](https://github.com/veenstrajelmer) in [#305](https://github.com/Deltares/dfm_tools/pull/305)
- added meshkernelpy helper functions like `make_basegrid()`, `refine_basegrid()` and `meshkernel_to_UgridDataset()` by [@veenstrajelmer](https://github.com/veenstrajelmer) in [#329](https://github.com/Deltares/dfm_tools/pull/329)

### Fix

- made bedlevel variable optional in `get_Dataset_atdepths()`, if not explicitly needing it by [@veenstrajelmer](https://github.com/veenstrajelmer) in [#241](https://github.com/Deltares/dfm_tools/pull/241)
