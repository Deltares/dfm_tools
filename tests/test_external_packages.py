# -*- coding: utf-8 -*-
"""
Created on Mon Jun 19 22:03:43 2023

@author: veenstra
"""

import os
import pytest
import xarray as xr
import xugrid as xu
import numpy as np
import matplotlib.pyplot as plt
import dfm_tools as dfmt
import contextily as ctx
from requests import HTTPError


@pytest.mark.unittest
def test_import_shapely():
    """
    tests whether shapely can be imported successfully, this is a problem in some environments
    in that case 'import shapely' works, but import 'shapely.geometry' fails
    """
    import shapely
    import shapely.geometry


@pytest.mark.unittest
def test_import_numba():
    """
    for some reason, numba often fails importing if versions of numpy/scipy are a bit too new
    """
    import numba as nb

    
@pytest.mark.unittest
def test_modplot_velovect():
    """
    this test will fail with matplotlib<3.6.0
    """
    x = np.linspace(-4,4,120)
    y = np.linspace(-3,3,100)
    X,Y = np.meshgrid(x,y)
    U = -1 - X**2 + Y
    V = 1 + X - Y**2
    speed = np.sqrt(U*U + V*V)
    grains = 15
    
    # dfmt.velovect requires matplotlib>=3.4.0
    fig,ax = plt.subplots()
    dfmt.velovect(ax, X, Y, U, V, color=speed, cmap='winter', arrowstyle='fancy', 
                  linewidth=speed/5, integration_direction='forward',
                  density=5, grains=grains)


@pytest.mark.unittest
def test_matplotlib_streamplot_broken_streamlines():
    """
    this test will fail with matplotlib<3.6.0
    """
    x = np.linspace(-4,4,120)
    y = np.linspace(-3,3,100)
    X,Y = np.meshgrid(x,y)
    U = -1 - X**2 + Y
    V = 1 + X - Y**2
    speed = np.sqrt(U*U + V*V)
    
    # broken_streamlines requires matplotlib>=3.6.0
    fig,ax = plt.subplots()
    ax.streamplot(X, Y, U, V, color=speed, cmap='winter', arrowstyle='fancy', linewidth=speed/5, broken_streamlines=False)


@pytest.mark.unittest
def test_xugrid_opendataset_ugridplot():
    """
    this one used to fail with xarray>=2023.3.0: https://github.com/Deltares/xugrid/issues/78
    This will probably not happen again, but it is convenient to test anyway.
    """
    file_nc = dfmt.data.fm_curvedbend_map(return_filepath=True)
    uds = xu.open_dataset(file_nc,chunks={'time':1})
    uds['mesh2d_flowelem_bl'].ugrid.plot()


@pytest.mark.unittest
def test_xugrid_opendataset_ugridplot_contourf():
    """
    This testcase gave a DeprecationWarning with scipy<1.10.0: https://github.com/Deltares/dfm_tools/issues/557
    Contourf/contour resulted in several other issues,
    so it is useful to test if the function can be called
    this test will fail with matplotlib 3.4.0 and 3.5.0
    """
    file_nc = dfmt.data.fm_curvedbend_map(return_filepath=True)
    uds = xu.open_dataset(file_nc)
    uds['mesh2d_flowelem_bl'].ugrid.plot.contourf()


@pytest.mark.unittest
def test_xugrid_opendataset_ugridplot_contour_with_colorbar():
    """
    Plotting a contour plot on uniform data (bedlevel -10 meter everywhere),
    resulted in several errors when adding a colorbar with matplotlib<=3.6.0
    """
    file_nc = dfmt.data.fm_curvedbend_map(return_filepath=True)
    uds = xu.open_dataset(file_nc)
    uds['mesh2d_flowelem_bl'].ugrid.plot.contour(add_colorbar=True)


@pytest.mark.unittest
def test_xarray_pandas_resample():
    """
    this one fails with xarray<2023.3.0 (required for test_opendataset_ugridplot()) and pandas>1.5.3: https://github.com/Deltares/xugrid/issues/78#issuecomment-1597723955
    pandas<1.5.3 works with xarray (any version)
    pandas>=2.0.0 fails with xarray<2023.3.0
    pandas>=2.0.0 cannot be installed with xarray==2023.3.0 (latter requires pandas<2 and >=1.4)
    pandas>=2.0.0 works with xarray>=2023.4.0
    therefore, xarray>=2023.4.0 works pandas (any version)
    """
    ds = xr.tutorial.load_dataset("air_temperature")
    ds.resample(time='D')

    
@pytest.mark.unittest
def test_xarray_decode_default_fillvals(tmp_path):
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
    fill_value_default = default_fillvals[fnc_dtype]
    fill_value = ds[varn_fnc].attrs.pop('_FillValue')
    ds[varn_fnc] = ds[varn_fnc].where(ds[varn_fnc]!=fill_value,fill_value_default)
    
    #write file
    file_out = os.path.join(tmp_path, 'temp_fnc_default_fillvals_map.nc')
    ds.to_netcdf(file_out)
    
    #open dataset with decode_fillvals
    with pytest.raises(ValueError) as e:
        # until xarray handles default fillvalues: https://github.com/Deltares/dfm_tools/issues/490
        # if xarray handles default fillvalues, `dfmt.decode_default_fillvals()` can be removed if minimum xarray version is set as requirement
        dfmt.open_partitioned_dataset(file_out, decode_fillvals=False)
    assert "connectivity contains negative values" in str(e.value)

    #this should be successful
    uds2 = dfmt.open_partitioned_dataset(file_out, decode_fillvals=True)
    fnc_new = uds2.grid.face_node_connectivity
    
    assert -1 in fnc_new


@pytest.mark.unittest
def test_ctx_add_basemap():
    """
    tests wheter this often used function still works, it fails in contextily 1.3.0
    mainly relevant to update the minimal contextily version in the pyproject.toml file sometimes
    https://github.com/Deltares/dfm_tools/issues/857
    """
    _, ax = plt.subplots()
    ax.set_xlim(1.8,3.2)
    ax.set_ylim(50.8, 52.2)
    try:
        ctx.add_basemap(ax=ax, crs="EPSG:4326", attribution=False)
    except HTTPError as e:
        # servers are flaky so HTTPError is raised quite often
        # this should not result in a failing dfm_tools test since it is just temporary
        print(e)
    except AttributeError as e:
        # can be removed after fixing https://github.com/geopandas/contextily/issues/252
        assert str(e) == "'NoneType' object has no attribute 'shape'"
        print(e)
    except TypeError as e:
        # can be removed after fixing https://github.com/geopandas/contextily/issues/252
        assert "not 'NoneType'" in str(e)
        print(e)
