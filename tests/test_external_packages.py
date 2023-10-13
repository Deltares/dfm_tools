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
    this one used to fail with xarray>=2023.3.0: https://github.com/Deltares/xugrid/issues/78
    This will probably not happen again, but it is convenient to test anyway.
    """
    file_nc = dfmt.data.fm_curvedbend_map(return_filepath=True)
    
    uds = xu.open_dataset(file_nc,chunks={'time':1})
    
    uds['mesh2d_flowelem_bl'].ugrid.plot()


@pytest.mark.systemtest
def test_xugrid_opendataset_ugridplot_contourf():
    """
    This testcase gave a DeprecationWarning with scipy<1.10.0: https://github.com/Deltares/dfm_tools/issues/557
    Contourf/contour resulted in several other issues,
    so it is useful to test if the function can be called
    """
    file_nc = dfmt.data.fm_curvedbend_map(return_filepath=True)
    uds = xu.open_dataset(file_nc)
    uds['mesh2d_flowelem_bl'].ugrid.plot.contourf()


@pytest.mark.systemtest
def test_xugrid_opendataset_ugridplot_contour_with_colorbar():
    """
    Plotting a contour plot on uniform data (bedlevel -10 meter everywhere),
    resulted in several errors upon contour, with matplotlib 3.6.0 when adding a colorbar
    """
    file_nc = dfmt.data.fm_curvedbend_map(return_filepath=True)
    uds = xu.open_dataset(file_nc)
    uds['mesh2d_flowelem_bl'].ugrid.plot.contour(add_colorbar=True)


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
