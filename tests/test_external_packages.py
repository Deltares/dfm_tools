# -*- coding: utf-8 -*-
"""
Created on Mon Jun 19 22:03:43 2023

@author: veenstra
"""

import pytest
import numpy as np
import xarray as xr
import xugrid as xu
import dfm_tools as dfmt


@pytest.mark.unittest
def test_import_shapely():
    """
    tests whether shapely can be imported successfully, this is a problem in some environments
    in that case 'import shapely' works, but import 'shapely.geometry' fails
    """
    import shapely
    import shapely.geometry
    

@pytest.mark.systemtest
def test_xugrid_opendataset_ugridplot():
    """
    this one fails with xarray>=2023.3.0: https://github.com/Deltares/xugrid/issues/78
    """
    file_nc = dfmt.data.fm_curvedbend_map(return_filepath=True)
    
    uds = xu.open_dataset(file_nc,chunks={'time':1})
    
    uds['mesh2d_flowelem_bl'].ugrid.plot()


@pytest.mark.systemtest
def SKIP_test_xugrid_opendataset_ugridplot_contourf():
    """
    this one fails with xarray>=2023.3.0: https://github.com/Deltares/xugrid/issues/117
    """
    file_nc = dfmt.data.fm_curvedbend_map(return_filepath=True)
    
    uds = xu.open_dataset(file_nc,chunks={'time':1})
    
    uds['mesh2d_flowelem_bl'].ugrid.plot.contourf()


@pytest.mark.unittest
def test_xarray_pandas_resample():
    """
    this one fails with xarray<2023.3.0 (required for test_opendataset_ugridplot()) and pandas>1.5.3: https://github.com/Deltares/xugrid/issues/78#issuecomment-1597723955
    pandas 1.5.3 works with xarray<2023.3.0
    pandas 2.0.* fails with xarray<2023.3.0
    pandas 2.0.* works with xarray==2023.5.0
    """
    ds = xr.tutorial.load_dataset("air_temperature")
    ds.resample(time='D')

    
@pytest.mark.unittest
def test_xarray_interp_to_newdim():
    """
    this one fails with scipy>=1.10.0: https://github.com/pydata/xarray/issues/7701
    """
    ds = xr.Dataset()
    so_np = np.array([[[35.819576, 35.82568 , 35.82873 ],
                       [35.819576, 35.824154, 35.831783],
                       [35.822628, 35.824154, 35.82873 ]],
                      
                      [[35.802788, 35.80584 , 35.815   ],
                       [35.815   , 35.810417, 35.821102],
                       [35.824154, 35.813473, 35.81805 ]],
                      
                      [[35.786003, 35.789055, np.nan],
                       [35.807365, 35.796684, np.nan],
                       [35.824154, 35.80584 , np.nan]],
                      
                      [[35.776848, np.nan,    np.nan],
                       [35.792107, np.nan,    np.nan],
                       [35.822628, np.nan,    np.nan]],
                                              
                      [[35.781425, np.nan,    np.nan],
                       [35.792107, np.nan,    np.nan],
                       [35.789055, np.nan,    np.nan]]])
    ds['so'] = xr.DataArray(so_np,dims=('depth','latitude','longitude'))
    ds['longitude'] = xr.DataArray([-9.6, -9.5, -9.4], dims=('longitude'))
    ds['latitude'] = xr.DataArray([42.9, 43.0, 43.1], dims=('latitude'))
    
    x_xr = xr.DataArray([-9.5],dims=('plipoints'))
    y_xr = xr.DataArray([43],dims=('plipoints'))
    
    interp_with_floats = ds.interp(longitude=x_xr[0], latitude=y_xr[0], method='linear').so #selecting one value from the da drops the new plipoints dimension
    interp_with_da_existing = ds.interp(longitude=x_xr.values, latitude=y_xr.values, method='linear').so.isel(longitude=0,latitude=0) #using the DataArray values keeps lat/lon dimenions, gives the same interp result
    interp_with_da_newdim = ds.interp(longitude=x_xr, latitude=y_xr, method='linear').so.isel(plipoints=0) #using the DataArray introduces a plipoints dimension, which gives different interp result
    print(interp_with_floats.to_numpy())
    print(interp_with_da_existing.to_numpy())
    print(interp_with_da_newdim.to_numpy())
    print(xr.__version__)
    
    assert (interp_with_floats.isnull()==interp_with_da_existing.isnull()).all() #success
    assert (interp_with_floats.isnull()==interp_with_da_newdim.isnull()).all() #fails with scipy>=1.10.0


@pytest.mark.unittest
def test_xarray_decode_default_fillvals():
    """
    This test will fail as soon as xarray handles default fillvalues: https://github.com/Deltares/dfm_tools/issues/490
    After that, the minimum xarray requirement can be updated
    However, py38 support must then be dropped: https://github.com/Deltares/dfm_tools/issues/267
    In that case, this testcase and `dfmt.decode_default_fillvals()` can be removed 
    """
    
    import dfm_tools as dfmt
    import xarray as xr
    from netCDF4 import default_fillvals
    
    file_nc = dfmt.data.fm_grevelingen_map(return_filepath=True)
    file_nc = file_nc.replace('_0*','_0002')
    
    ds = xr.open_dataset(file_nc,decode_cf=False)
    
    #convert fillvalue in fnc to default fillvalue
    varn_fnc = 'mesh2d_face_nodes'
    fnc_dtype = ds[varn_fnc].dtype.str[1:]
    fill_value = ds[varn_fnc].attrs['_FillValue']
    fill_value_default = default_fillvals[fnc_dtype]
    ds[varn_fnc].attrs.pop('_FillValue')
    ds[varn_fnc] = ds[varn_fnc].where(ds[varn_fnc]!=fill_value,fill_value_default)
    
    #write file
    file_out = 'fnc_default_fillvals_map.nc'
    ds.to_netcdf(file_out)
    ds.close()
    
    #open dataset with decode_fillvals
    try:
        uds = dfmt.open_partitioned_dataset(file_out,decode_fillvals=False)
    except Exception as e:
        # this raises "ValueError: connectivity contains negative values"
        # until xarray handles default fillvalues: https://github.com/Deltares/dfm_tools/issues/490
        assert isinstance(e,ValueError)
    if 'uds' in locals():
        raise Exception("apparently xarray now decodes default fillvalues, so "
                        "`dfmt.decode_default_fillvals()` can be removed if minimum xarray "
                        "version is set as requirement")
    
    #this should be successful
    uds = dfmt.open_partitioned_dataset(file_out,decode_fillvals=True)
    fnc_new = uds.grid.face_node_connectivity
    
    assert fill_value_default in fnc_new
