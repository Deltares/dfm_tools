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

### Feat

- `open_partitioned_dataset()` replaces get_ncmodeldata and related functions, now reading of ugrid (D-FlowFM, D-HYDRO) files is done with [xugrid](https://github.com/Deltares/xugrid) and [xarray](https://github.com/pydata/xarray).
- create release and pypipip installable by [@veenstrajelmer](https://github.com/veenstrajelmer) in [#338](https://github.com/Deltares/dfm_tools/pull/338)
- [unfiltered release notes](https://github.com/Deltares/dfm_tools/releases/tag/v0.11.0)
