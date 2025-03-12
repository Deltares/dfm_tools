# -*- coding: utf-8 -*-
"""
Created on Tue Jan 23 17:12:48 2024

@author: veenstra
"""

import os
import pytest
import dfm_tools as dfmt
from dfm_tools.errors import OutOfRangeError
from dfm_tools.xarray_helpers import file_to_list
import pandas as pd
import xarray as xr
import numpy as np
from pathlib import Path
from dfm_tools.xarray_helpers import interpolate_na_multidim


@pytest.mark.unittest
def test_preprocess_era5_valid_time(ds_era5_empty):
    ds_era5_empty = ds_era5_empty.rename(time="valid_time")

    ds = dfmt.preprocess_ERA5(ds_era5_empty)
    assert "valid_time" in ds_era5_empty.dims
    assert "time" not in ds_era5_empty.dims
    assert "valid_time" not in ds.dims
    assert "time" in ds.dims


@pytest.mark.unittest
def test_preprocess_era5_expver_coord(ds_era5_empty):
    ds = dfmt.preprocess_ERA5(ds_era5_empty)
    assert "expver" not in ds_era5_empty.coords
    assert "expver" in ds.coords


@pytest.mark.unittest
def test_preprocess_era5_expver_dim(ds_era5_empty):
    ntimes = len(ds_era5_empty.time)
    data_dummy = np.ones(shape=(ntimes,2))
    ds_era5_empty['dummy'] = xr.DataArray(data_dummy, dims=('time','expver'))

    ds = dfmt.preprocess_ERA5(ds_era5_empty)
    assert "expver" in ds_era5_empty.dims
    assert "expver" not in ds.dims


@pytest.mark.unittest
def test_preprocess_era5_mer_mtpr(ds_era5_empty):
    ds_era5_empty['avg_tprate'] = xr.DataArray()
    ds_era5_empty['avg_ie'] = xr.DataArray()
    ds = dfmt.preprocess_ERA5(ds_era5_empty)
    assert "avg_tprate" in ds_era5_empty.data_vars
    assert "avg_ie" in ds_era5_empty.data_vars
    assert "mtpr" in ds.data_vars
    assert "mer" in ds.data_vars


@pytest.mark.unittest
def test_preprocess_era5_int32(ds_era5_empty):
    ds_era5_empty['dummy_int'] = xr.DataArray()
    ds_era5_empty['dummy_int'].encoding['dtype'] = 'int32'
    ds_era5_empty['dummy_int'].encoding['scale_factor'] = 1
    ds_era5_empty['dummy_int'].encoding['add_offset'] = 1
    
    # assertions fail after preprocess_ERA5 (attrs are popped from input ds)
    assert "dtype" in ds_era5_empty['dummy_int'].encoding.keys()
    assert "scale_factor" in ds_era5_empty['dummy_int'].encoding.keys()
    assert "add_offset" in ds_era5_empty['dummy_int'].encoding.keys()
    
    ds = dfmt.preprocess_ERA5(ds_era5_empty)
    
    assert ds['dummy_int'].encoding['dtype'] == 'float32'
    assert "scale_factor" not in ds['dummy_int'].encoding.keys()
    assert "add_offset" not in ds['dummy_int'].encoding.keys()


@pytest.mark.unittest
@pytest.mark.requiressecrets
@pytest.mark.era5slow # temporarily skip these on github
@pytest.mark.timeout(60) # useful since CDS downloads are terribly slow sometimes, so skip in that case
def test_merge_meteofiles(file_nc_era5_pattern):
    # file_nc_era5_pattern comes from file_nc_era5_pattern() in conftest.py
    # deliberately take time_slice.stop as a non-existing timestep to check
    # outside bounds
    ds = dfmt.merge_meteofiles(
        file_nc=file_nc_era5_pattern, 
        preprocess=dfmt.preprocess_ERA5, 
        time_slice=slice("2010-01-30","2010-02-01 22:30")
        )
    assert ds.sizes["time"] == 72
    assert ds.time.to_pandas().iloc[0] == pd.Timestamp('2010-01-30')
    assert ds.time.to_pandas().iloc[-1] == pd.Timestamp('2010-02-01 23:00')
    assert "msl" in ds.data_vars


@pytest.mark.unittest
def test_merge_meteofiles_outofrange_times(ds_era5_empty, tmp_path):
    file_nc = os.path.join(tmp_path, "era5_msl_empty.nc")
    ds_era5_empty.to_netcdf(file_nc)
    file_nc_pat = os.path.join(tmp_path, "*.nc")

    date_min = "2030-01-01"
    date_max = "2030-02-01"

    # merge meteo
    with pytest.raises(OutOfRangeError) as e:
        _ = dfmt.merge_meteofiles(
            file_nc=file_nc_pat,
            preprocess=dfmt.preprocess_ERA5, 
            time_slice=slice(date_min, date_max),
            )
    assert "requested tstop 2030-02-01 00:00:00 outside" in str(e.value)


@pytest.mark.unittest
def test_merge_meteofiles_duplicated_times(ds_era5_empty, tmp_path):
    file_nc = os.path.join(tmp_path, "era5_msl_empty.nc")
    ds = xr.concat([ds_era5_empty,ds_era5_empty], dim='time')
    ds.to_netcdf(file_nc)
    file_nc_pat = os.path.join(tmp_path, "*.nc")

    date_min = "2010-01-31"
    date_max = "2010-02-01"

    # merge meteo
    ds_merged = dfmt.merge_meteofiles(
        file_nc=file_nc_pat,
        preprocess=dfmt.preprocess_ERA5, 
        time_slice=slice(date_min, date_max),
        )
    
    assert len(ds.time) == 18
    assert len(ds_merged.time) == 9


@pytest.mark.unittest
def test_merge_meteofiles_times_gap(ds_era5_empty, tmp_path):
    file_nc = os.path.join(tmp_path, "era5_msl_empty.nc")
    ds_era5_empty = ds_era5_empty.isel(time=[0,1,2,3,6,7,8])
    ds_era5_empty.to_netcdf(file_nc)
    file_nc_pat = os.path.join(tmp_path, "*.nc")

    date_min = "2010-01-31"
    date_max = "2010-02-01"

    # merge meteo
    with pytest.raises(ValueError) as e:
        _ = dfmt.merge_meteofiles(
            file_nc=file_nc_pat,
            preprocess=dfmt.preprocess_ERA5, 
            time_slice=slice(date_min, date_max),
            )
    assert "time gaps found in selected dataset" in str(e.value)


@pytest.mark.unittest
def test_merge_meteofiles_rename_latlon(ds_era5_empty, tmp_path):
    date_min = "2010-01-31"
    date_max = "2010-02-01"
    
    # lat/lon latitude/longitude vars
    ds = ds_era5_empty.copy()
    ds = ds.rename({'longitude':'lon', 'latitude':'lat'})
    file_nc = os.path.join(tmp_path, "era5_lonlat_empty.nc")
    ds.to_netcdf(file_nc)
    file_nc_pat = file_nc.replace(".nc", "*.nc")
    ds_merged = dfmt.merge_meteofiles(
        file_nc=file_nc_pat,
        preprocess=dfmt.preprocess_ERA5, 
        time_slice=slice(date_min, date_max),
        )
    assert "longitude" in ds_merged.data_vars
    assert "latitude" in ds_merged.data_vars
    
    # x/y latitude/longitude vars
    ds = ds_era5_empty.copy()
    ds = ds.rename({'longitude':'x', 'latitude':'y'})
    file_nc = os.path.join(tmp_path, "era5_xy_empty.nc")
    ds.to_netcdf(file_nc)
    file_nc_pat = file_nc.replace(".nc", "*.nc")
    ds_merged = dfmt.merge_meteofiles(
        file_nc=file_nc_pat,
        preprocess=dfmt.preprocess_ERA5, 
        time_slice=slice(date_min, date_max),
        )
    assert "longitude" in ds_merged.data_vars
    assert "latitude" in ds_merged.data_vars

    # no latitude/longitude vars
    ds = ds_era5_empty.copy()
    ds = ds.drop_vars(['longitude', 'latitude'])
    file_nc = os.path.join(tmp_path, "era5_none_empty.nc")
    ds.to_netcdf(file_nc)
    file_nc_pat = file_nc.replace(".nc", "*.nc")
    with pytest.raises(KeyError) as e:
        _ = dfmt.merge_meteofiles(
            file_nc=file_nc_pat,
            preprocess=dfmt.preprocess_ERA5, 
            time_slice=slice(date_min, date_max),
            )
    assert "no longitude/latitude, lon/lat or x/y variables" in str(e.value)


@pytest.mark.unittest
def test_merge_meteofiles_convert360to180(ds_era5_empty, tmp_path):
    file_nc = os.path.join(tmp_path, "era5_msl_empty.nc")
    lon_vals = np.arange(0, 360, 0.5) # from -180 to 179.5
    ds_era5_empty = ds_era5_empty.drop_vars(['longitude'])
    ds_era5_empty['longitude'] = xr.DataArray(lon_vals, dims='longitude')
    ds_era5_empty.to_netcdf(file_nc)
    file_nc_pat = os.path.join(tmp_path, "*.nc")
    
    date_min = "2010-01-31"
    date_max = "2010-02-01"

    # merge meteo
    ds_merged = dfmt.merge_meteofiles(
        file_nc=file_nc_pat,
        preprocess=dfmt.preprocess_ERA5, 
        time_slice=slice(date_min, date_max),
        )
    assert np.isclose(ds_era5_empty["longitude"][0], 0)
    assert np.isclose(ds_era5_empty["longitude"][-1], 359.5)
    assert np.isclose(ds_merged["longitude"][0], -180.0)
    assert np.isclose(ds_merged["longitude"][-1], 179.5)
 
    
@pytest.mark.unittest
def test_merge_meteofiles_number_coordinate(ds_era5_empty, tmp_path):
    """
    The number coordinate value does not cause issues per definition. However, when it
    is present in the first and last files of a set of files, but not in the
    intermediate ones, we get "ValueError: 'number' not present in all datasets and 
    coords='different'. [...]". For more details, see
    https://github.com/Deltares/dfm_tools/issues/1156.
    coords='minimal was added to '
    """
    ds1 = ds_era5_empty.isel(time=slice(None,3))
    ds2 = ds_era5_empty.isel(time=slice(3,6))
    ds3 = ds_era5_empty.isel(time=slice(6,None))
    ds1['number'] = 0
    ds1 = ds1.set_coords('number')
    ds3['number'] = 0
    ds3 = ds3.set_coords('number')
    ds1.to_netcdf(os.path.join(tmp_path, "file1.nc"))
    ds2.to_netcdf(os.path.join(tmp_path, "file2.nc"))
    ds3.to_netcdf(os.path.join(tmp_path, "file3.nc"))
    
    time_slice = slice('2010-01-31','2010-02-01')
    file_nc = os.path.join(tmp_path, "*.nc")
    kwargs = dict(preprocess=dfmt.preprocess_ERA5)
    _ = dfmt.merge_meteofiles(
        file_nc=file_nc,
        time_slice=time_slice, 
        **kwargs,
        )
    

@pytest.mark.unittest
def test_file_to_list_pathlib_path():
    file_nc = dfmt.data.fm_curvedbend_his(return_filepath=True)
    file_nc_list = file_to_list(Path(file_nc))
    assert file_nc_list == [file_nc]


@pytest.mark.unittest
def test_file_to_list_empty_path():
    file_nc = "path/to/dummy/dir"
    with pytest.raises(FileNotFoundError) as e:
        _ = file_to_list(file_nc)
    assert "file(s) not found, empty file_nc_list" in str(e.value)


@pytest.mark.unittest
def test_file_to_list_already_list():
    file_nc = dfmt.data.fm_curvedbend_his(return_filepath=True)
    file_nc_list = file_to_list([file_nc])
    assert file_nc_list == [file_nc]


@pytest.mark.unittest
def test_preprocess_hisnc():
    """
    not too much added value, but good to check dropping of duplicated labels
    in this case it happens for source_sinks, not sure if this is useful here.
    More useful would be a hisfile with duplicated station names like
    p:\\archivedprojects\\11206813-006-kpp2021_rmm-2d\\C_Work\\31_RMM_FMmodel
    \\computations\\model_setup\\run_206\\results\\RMM_dflowfm_0000_his.nc
    """
    file_nc = dfmt.data.fm_grevelingen_his(return_filepath=True)
    ds1 = xr.open_dataset(file_nc)
    ds2 = xr.open_mfdataset(file_nc, preprocess=dfmt.preprocess_hisnc)
    assert ds1.sizes['source_sink'] == 46
    assert ds2.sizes['source_sink'] == 1


@pytest.mark.unittest
def test_interpolate_na_multidim(cmems_dataset_notime):
    """
    check if the interpolation works as expected
    """
    
    ds = cmems_dataset_notime.copy()
    for var in ds.data_vars:
        ds[var] = interpolate_na_multidim(ds[var], ["latitude", "longitude"])
        ds[var] = interpolate_na_multidim(ds[var], ["depth"])

    expected_vals = np.array(
        [[[35.819576, 35.82568 , 35.82873 ],
        [35.819576, 35.824154, 35.831783],
        [35.822628, 35.824154, 35.82873 ]],

       [[35.802788, 35.80584 , 35.815   ],
        [35.815   , 35.810417, 35.821102],
        [35.824154, 35.813473, 35.81805 ]],

       [[35.786003, 35.789055, 35.789055],
        [35.807365, 35.796684, 35.796684],
        [35.824154, 35.80584 , 35.80584 ]],

       [[35.776848, 35.776848, 35.776848],
        [35.792107, 35.792107, 35.792107],
        [35.822628, 35.822628, 35.822628]],

       [[35.776848, 35.776848, 35.776848],
        [35.792107, 35.792107, 35.792107],
        [35.822628, 35.822628, 35.822628]]])
    assert np.allclose(ds.so.to_numpy(), expected_vals)
