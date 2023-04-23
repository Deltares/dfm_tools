# -*- coding: utf-8 -*-
"""
Created on Sun Apr 23 22:42:27 2023

@author: veenstra
"""


def plot_background(ax=None, projection=None, google_style='satellite', resolution=1, features=None, nticks=6, latlon_format=False, gridlines=False, **kwargs):
    raise DeprecationWarning('dfmt.plot_background() is deprecated, use contextily and cartopy instead like in https://github.com/Deltares/dfm_tools/blob/main/tests/examples_workinprogress/workinprogress_cartopy_satellite_coastlines.py')


def get_ugrid_verts(data_xr_map):
    """
    getting ugrid verts from xugrid mapfile.
    """
    raise DeprecationWarning('dfmt.get_ugrid_verts() is deprecated, use uds.grid.face_node_coordinates instead (https://github.com/Deltares/xugrid/issues/48)')


def scatter_to_regulargrid(xcoords, ycoords, values, ncellx=None, ncelly=None, reg_x_vec=None, reg_y_vec=None, method='nearest', maskland_dist=None):
    raise DeprecationWarning('dfm_tools.regulargrid.scatter_to_regulargrid() is deprecated, use ds = dfmt.rasterize_ugrid(uds) instead')


def get_varnamefromattrs(data_xr, varname):
    raise DeprecationWarning('dfmt.get_varnamefromattrs() will be deprecated in a future version of dfm_tools, ds=dfmt.rename_waqvars(ds) is a more convenient alternative')


class Polygon:
    def __init__(self, data, name, comments):
        raise DeprecationWarning('the function dfm_tools.polygon.Polygon() is deprecated, please use the new hydrolib alternative.')
        
    def fromfile(self, file_pol, pd_output=False, tekmap_output=False):
        raise DeprecationWarning('the function dfm_tools.polygon.Polygon.fromfile() is deprecated, please use the new hydrolib alternative. Example script: https://github.com/Deltares/dfm_tools/blob/main/tests/examples/preprocess_hydrolib_readwritepol.py')


def write_bcfile(filename, datablocks, metadatas, refdate=None, tzone=0, float_format='%6.2f'):
    raise DeprecationWarning('the function dfm_tools.io.bc.write_bcfile() is deprecated, please use the new hydrolib alternative. Example script: dfm_tools/tests/examples/CMEMS_interpolate_example.py')


def read_bcfile(filename, converttime=False):
    raise DeprecationWarning('the function dfm_tools.io.bc.read_bcfile() is deprecated, please use the new hydrolib alternative. Example script: dfm_tools/tests/examples/hydrolib_readbc.py')


def write_timfile(filename, datablock, header, converttime=False, refdate=None, float_format='%6.2f'):
    raise DeprecationWarning('the function dfm_tools.write_timfile() is deprecated, please use the new hydrolib alternative: https://github.com/Deltares/dfm_tools/blob/301-convert-timmodel-to-pandasdataframe/tests/examples/preprocess_hydrolib_readtim.py.')


def read_timfile(filename, converttime=False, refdate=None):
    raise DeprecationWarning('the function dfm_tools.read_timfile() is deprecated, please use the new hydrolib alternative: https://github.com/Deltares/dfm_tools/blob/301-convert-timmodel-to-pandasdataframe/tests/examples/preprocess_hydrolib_readtim.py.')

