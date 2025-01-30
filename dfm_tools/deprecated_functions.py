def get_ncmodeldata(file_nc=None, varname=None, timestep=None, layer=None, depth=None, station=None, multipart=None, get_linkedgridinfo=False, silent=False):
    raise DeprecationWarning('dfmt.get_ncmodeldata() is deprecated, use `uds = dfmt.open_partitioned_dataset(file_nc)` instead: https://github.com/Deltares/dfm_tools/blob/main/notebooks/postprocessing_example.ipynb')


def get_netdata(file_nc=None, multipart=None):
    raise DeprecationWarning('dfmt.get_netdata() is deprecated, use `uds = dfmt.open_partitioned_dataset(file_nc)` instead: https://github.com/Deltares/dfm_tools/blob/main/notebooks/postprocessing_example.ipynb')


def plot_netmapdata(verts=None, values=None, ax=None, **kwargs):
    raise DeprecationWarning('dfmt.plot_netmapdata() is deprecated, use `uds = dfmt.open_partitioned_dataset(file_nc); uds.mesh2d_s1.isel(time=0).ugrid.plot()` instead: https://github.com/Deltares/dfm_tools/blob/main/notebooks/postprocessing_example.ipynb')


def plot_background(ax=None, projection=None, google_style='satellite', resolution=1, features=None, nticks=6, latlon_format=False, gridlines=False, **kwargs):
    raise DeprecationWarning('dfmt.plot_background() is deprecated, use contextily instead')


def get_ugrid_verts(data_xr_map=None):
    raise DeprecationWarning('dfmt.get_ugrid_verts() is deprecated, use uds.grid.face_node_coordinates instead (https://github.com/Deltares/xugrid/issues/48)')


def scatter_to_regulargrid(xcoords=None, ycoords=None, values=None, ncellx=None, ncelly=None, reg_x_vec=None, reg_y_vec=None, method='nearest', maskland_dist=None):
    raise DeprecationWarning('dfm_tools.regulargrid.scatter_to_regulargrid() is deprecated, use ds = dfmt.rasterize_ugrid(uds) instead')


def get_varnamefromattrs(data_xr=None, varname=None):
    raise DeprecationWarning('dfmt.get_varnamefromattrs() will be deprecated in a future version of dfm_tools, ds=dfmt.rename_waqvars(ds) is a more convenient alternative')


class Polygon:
    def __init__(self, data=None, name=None, comments=None):
        raise DeprecationWarning('the function dfm_tools.polygon.Polygon() is deprecated, please use the new hydrolib alternative.')
        
    def fromfile(self, file_pol=None, pd_output=False, tekmap_output=False):
        raise DeprecationWarning('the function dfm_tools.polygon.Polygon.fromfile() is deprecated, please use the new hydrolib alternative. Example script: https://github.com/Deltares/dfm_tools/blob/main/tests/examples/preprocess_hydrolib_readwritepol.py')


def write_bcfile(filename=None, datablocks=None, metadatas=None, refdate=None, tzone=0, float_format='%6.2f'):
    raise DeprecationWarning('the function dfm_tools.io.bc.write_bcfile() is deprecated, please use the new hydrolib alternative. Example script: dfm_tools/tests/examples/CMEMS_interpolate_example.py')


def read_bcfile(filename=None, converttime=False):
    raise DeprecationWarning('the function dfm_tools.io.bc.read_bcfile() is deprecated, please use the new hydrolib alternative. Example script: dfm_tools/tests/examples/hydrolib_readbc.py')


def write_timfile(filename=None, datablock=None, header=None, converttime=False, refdate=None, float_format='%6.2f'):
    raise DeprecationWarning('the function dfm_tools.write_timfile() is deprecated, please use the new hydrolib alternative: https://github.com/Deltares/dfm_tools/blob/301-convert-timmodel-to-pandasdataframe/tests/examples/preprocess_hydrolib_readtim.py.')


def read_timfile(filename=None, converttime=False, refdate=None):
    raise DeprecationWarning('the function dfm_tools.read_timfile() is deprecated, please use the new hydrolib alternative: https://github.com/Deltares/dfm_tools/blob/301-convert-timmodel-to-pandasdataframe/tests/examples/preprocess_hydrolib_readtim.py.')


def generate_bndpli(**kwargs):
    raise DeprecationWarning('the function dfmt.generate_bndpli() is deprecated, please use dfmt.generate_bndpli_cutland() instead.')


def preprocess_hirlam(ds):
    raise DeprecationWarning('the function dfmt.preprocess_hirlam() is deprecated, xarray now supports datasets with multidimensional coordinates.')


def interp_regularnc_to_plipoints(ds):
    raise DeprecationWarning('the function dfmt.interp_regularnc_to_plipoints() is deprecated, '
                             'use dfmt.interp_regularnc_to_plipointsDataset() instead with gdf_points '
                             'as in https://github.com/Deltares/dfm_tools/issues/938')

def open_dataset_extra(**kwargs):
    raise DeprecationWarning('the function dfmt.open_dataset_extra() is deprecated, dfmt.open_prepare_dataset() is similar but does not support multiple quantities at once')
