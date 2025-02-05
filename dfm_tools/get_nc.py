import numpy as np
import datetime as dt
import re
import xugrid as xu
import xarray as xr
import matplotlib.pyplot as plt
from dfm_tools.xarray_helpers import Dataset_varswithdim
from dfm_tools.xugrid_helpers import get_vertical_dimensions, decode_default_fillvals

__all__ = ["polyline_mapslice",
           "reconstruct_zw_zcc",
           "get_Dataset_atdepths",
           "rasterize_ugrid",
           "plot_ztdata",
    ]

def calc_dist_pythagoras(x1,x2,y1,y2):
    distance = np.sqrt((x2 - x1)**2 + (y2 - y1)**2)
    return distance


def calc_dist_haversine(lon1,lon2,lat1,lat2):
    #calculates distance between lat/lon coordinates in meters
    #https://community.esri.com/t5/coordinate-reference-systems-blog/distance-on-a-sphere-the-haversine-formula/ba-p/902128
    
    # convert to radians
    lon1_rad = np.deg2rad(lon1)
    lon2_rad = np.deg2rad(lon2)
    lat1_rad = np.deg2rad(lat1)
    lat2_rad = np.deg2rad(lat2)
    
    # apply formulae
    a = np.sin((lat2_rad-lat1_rad)/2)**2 + np.cos(lat1_rad) * np.cos(lat2_rad) * np.sin((lon2_rad-lon1_rad)/2)**2
    c = 2 * np.arctan2( np.sqrt(a), np.sqrt(1-a) )
    R = 6371000
    distance = R * c
    if np.isnan(distance).any():
        raise ValueError('nan encountered in calc_dist_latlon distance')
    return distance


def intersect_edges_withsort(uds,edges): #TODO: move sorting to xugrid? https://deltares.github.io/xugrid/api/xugrid.Ugrid2d.intersect_edges.html
    
    edge_index, face_index, intersections = uds.grid.intersect_edges(edges) #TODO: is fast, but maybe speed can be increased with bounding box?
    
    #ordering of face_index is wrong (visible with cb3 with long line_array), so sort on distance from startpoint (in x/y units)
    
    #compute distance from start of line to start of each linepart
    edge_len = np.linalg.norm(edges[:,1] - edges[:,0], axis=1)
    edge_len_cum = np.cumsum(edge_len)
    edge_len_cum0 = np.concatenate([[0],edge_len_cum[:-1]])
    
    #compute distance from start to lineparts to start of line (via line)
    startcoord_linepart = edges[edge_index,0,:]
    dist_tostart_linepart = np.linalg.norm(intersections[:,0,:] - startcoord_linepart, axis=1)
    dist_tostart_line = dist_tostart_linepart + edge_len_cum0[edge_index]
    
    #sorting on distance
    id_sorted = np.argsort(dist_tostart_line)
    edge_index = edge_index[id_sorted]
    face_index = face_index[id_sorted]
    intersections = intersections[id_sorted]
    return edge_index, face_index, intersections


def get_xzcoords_onintersection(uds, face_index, crs_dist_starts, crs_dist_stops):
    #TODO: remove hardcoding of variable names
    if 'time' in uds.dims: #TODO: maybe make time dependent grid?
        raise Exception('time dimension present in uds, provide uds.isel(time=timestep) instead. This is necessary to retrieve correct waterlevel or fullgrid output')
    
    dimn_layer, dimn_interfaces = get_vertical_dimensions(uds)
    gridname = uds.grid.name
    
    #construct fullgrid info (zcc/zw) for 3D models
    if dimn_layer in uds.dims:
        uds = reconstruct_zw_zcc(uds)

    #drop all variables that do not contain a face dimension, then select only all sliced faceidx
    xu_facedim = uds.grid.face_dimension
    face_index_xr = xr.DataArray(face_index,dims=('ncrossed_faces'))
    uds = Dataset_varswithdim(uds,dimname=xu_facedim) #TODO: is there an xugrid alternative?
    uds_sel = uds.sel({xu_facedim:face_index_xr})
    
    # take zvals_interface
    if dimn_layer in uds_sel.dims: #3D model
        nlay = uds.sizes[dimn_layer]
        zvals_interface_filled = uds_sel[f'{gridname}_flowelem_zw'].bfill(dim=dimn_interfaces) #fill nan values (below bed) with equal values
        zvals_interface = zvals_interface_filled.to_numpy().T #transpose to make in line with 2D sigma dataset
    else: #2D model, no layers
        nlay = 1
        data_frommap_wl3_sel = uds_sel[f'{gridname}_s1'].to_numpy() #TODO: add escape for missing wl/bl vars
        data_frommap_bl_sel = uds_sel[f'{gridname}_flowelem_bl'].to_numpy()
        zvals_interface = np.linspace(data_frommap_bl_sel,data_frommap_wl3_sel,nlay+1)

    #derive crs_verts
    crs_dist_starts_matrix = np.repeat(crs_dist_starts[np.newaxis],nlay,axis=0)
    crs_dist_stops_matrix = np.repeat(crs_dist_stops[np.newaxis],nlay,axis=0)
    crs_verts_x_all = np.array([[crs_dist_starts_matrix.ravel(),crs_dist_stops_matrix.ravel(),crs_dist_stops_matrix.ravel(),crs_dist_starts_matrix.ravel()]]).T
    crs_verts_z_all = np.ma.array([zvals_interface[1:,:].ravel(),zvals_interface[1:,:].ravel(),zvals_interface[:-1,:].ravel(),zvals_interface[:-1,:].ravel()]).T[:,:,np.newaxis]
    crs_verts = np.ma.concatenate([crs_verts_x_all, crs_verts_z_all], axis=2)
    
    #define grid
    shape_crs_grid = crs_verts[:,:,0].shape
    shape_crs_flat = crs_verts[:,:,0].ravel().shape
    xr_crs_grid = xu.Ugrid2d(node_x=crs_verts[:,:,0].ravel(),
                             node_y=crs_verts[:,:,1].ravel(),
                             fill_value=-1,
                             face_node_connectivity=np.arange(shape_crs_flat[0]).reshape(shape_crs_grid),
                             )

    #define dataset
    if dimn_layer in uds_sel.dims:
        crs_plotdata_clean = uds_sel.stack({xr_crs_grid.face_dimension:[dimn_layer,'ncrossed_faces']},create_index=False)
    else: #2D: still make sure xr_crs_grid.face_dimension is created, using stack since .rename() gives "UserWarning: rename 'ncrossed_faces' to 'mesh2d_nFaces' does not create an index anymore."
        crs_plotdata_clean = uds_sel.stack({xr_crs_grid.face_dimension:['ncrossed_faces']},create_index=False)
                    
    #combine into xugrid
    xr_crs_ugrid = xu.UgridDataset(crs_plotdata_clean, grids=[xr_crs_grid])
    return xr_crs_ugrid


def polyline_mapslice(uds:xu.UgridDataset, line_array:np.array) -> xu.UgridDataset:
    """
    Slice trough mapdata, combine: intersect_edges_withsort, calculation of distances and conversion to ugrid dataset.

    Parameters
    ----------
    uds : xu.UgridDataset
        DESCRIPTION.
    line_array : np.array
        DESCRIPTION.
    calcdist_fromlatlon : bool, optional
        DESCRIPTION. The default is None.

    Raises
    ------
    Exception
        DESCRIPTION.

    Returns
    -------
    xr_crs_ugrid : xu.UgridDataset
        DESCRIPTION.

    """
    
    #compute intersection coordinates of crossings between edges and faces and their respective indices
    edges = np.stack([line_array[:-1],line_array[1:]],axis=1)
    edge_index, face_index, intersections = intersect_edges_withsort(uds=uds, edges=edges)
    if len(edge_index) == 0:
        raise ValueError('polyline does not cross mapdata')
    
    if uds.grids[0].is_geographic:
        calc_dist = calc_dist_haversine
    else:
        calc_dist = calc_dist_pythagoras
    
    #compute pyt/haversine start/stop distances for all intersections
    edge_len = calc_dist(edges[:,0,0], edges[:,1,0], edges[:,0,1], edges[:,1,1])
    edge_len_cum = np.cumsum(edge_len)
    edge_len_cum0 = np.concatenate([[0],edge_len_cum[:-1]])
    crs_dist_starts = calc_dist(edges[edge_index,0,0], intersections[:,0,0], edges[edge_index,0,1], intersections[:,0,1]) + edge_len_cum0[edge_index]
    crs_dist_stops  = calc_dist(edges[edge_index,0,0], intersections[:,1,0], edges[edge_index,0,1], intersections[:,1,1]) + edge_len_cum0[edge_index]
    
    #derive vertices from cross section (distance from first point)
    xr_crs_ugrid = get_xzcoords_onintersection(uds=uds, face_index=face_index, crs_dist_starts=crs_dist_starts, crs_dist_stops=crs_dist_stops)
    
    return xr_crs_ugrid


def get_formula_terms(uds, varn_contains):
    """
    get formula_terms for zw/zcc reconstruction, convert to list and then to dict. This can be done for layer/interface (via varn_contains)
    """
    osz_varnames = list(uds.filter_by_attrs(formula_terms=lambda v: v is not None).variables) #names of variables containing attribute "formula_terms"
    osz_varnames_contains = [x for x in osz_varnames if varn_contains in x] #TODO: to get the layer/interface ocean_*_coordinate. Not too pretty, but it works
    if len(osz_varnames_contains) != 1: #should be 1 exactly, none is the case in zlayer models
        raise ValueError(f'no or more than one {varn_contains} variable found with formula_terms attribute: {osz_varnames}')
    osz_varn = osz_varnames_contains[0]
    osz_formulaterms = uds[osz_varn].attrs['formula_terms']
    tokens = re.split('[:\\s]+', osz_formulaterms)
    osz_formulaterms_dict = dict(zip(tokens[::2], tokens[1::2]))
    return osz_formulaterms_dict


def reconstruct_zw_zcc_fromsigma(uds):
    """
    reconstruct full grid output (time/face-varying z-values) for sigma model, necessary for slicing sigmamodel on depth value
    based on https://cfconventions.org/cf-conventions/cf-conventions.html#_ocean_sigma_coordinate
    """
    osz_formulaterms_int_dict = get_formula_terms(uds,varn_contains='interface')
    osz_formulaterms_lay_dict = get_formula_terms(uds,varn_contains='layer')
    gridname = uds.grid.name
    
    uds_eta = uds[osz_formulaterms_int_dict['eta']] #mesh2d_s1
    uds_depth = uds[osz_formulaterms_int_dict['depth']] #mesh2d_bldepth in new output, mesh2d_waterdepth in old output (see comment below)
    if uds_depth.attrs['standard_name'] == 'sea_floor_depth_below_sea_surface': # previously the waterdepth instead of negative bedlevel was coupled via the formula_terms in sigmamodels (was fixed in OSS 140982 / 29-3-2022)
        uds_depth = -uds[f'{gridname}_flowelem_bl'] # assuming this variable is available, which is not guaranteed
    uds_sigma_int = uds[osz_formulaterms_int_dict['sigma']] #mesh2d_interface_sigma
    uds_sigma_lay = uds[osz_formulaterms_lay_dict['sigma']] #mesh2d_layer_sigma
    
    uds[f'{gridname}_flowelem_zw'] = uds_eta + uds_sigma_int*(uds_depth+uds_eta)
    uds[f'{gridname}_flowelem_zcc'] = uds_eta + uds_sigma_lay*(uds_depth+uds_eta)
    
    uds = uds.set_coords([f'{gridname}_flowelem_zw',f'{gridname}_flowelem_zcc'])
    return uds


def reconstruct_zw_zcc_fromz(uds):
    """
    reconstruct full grid output (time/face-varying z-values) for zvalue model. Necessary when extracting values with zdepth w.r.t. waterlevel/bedlevel
    """
    
    dimn_layer, dimn_interfaces = get_vertical_dimensions(uds)
    gridname = uds.grid.name
    
    uds_eta = uds[f'{gridname}_s1'] # assuming this variable is available, which is not guaranteed
    uds_z0 = xu.zeros_like(uds_eta)
    uds_bl = uds[f'{gridname}_flowelem_bl'] # assuming this variable is available, which is not guaranteed
    
    #deriving zcenter values, clipping zcc to bl/wl
    # zvals_center = uds[f'{gridname}_layer_z']
    # zcc = (uds_z0+zvals_center).clip(min=uds_bl, max=uds_eta)

    #deriving zinterface values, first expanding zint to wl.max(), then clipping zw to bl/wl
    zvals_interface = uds[f'{gridname}_interface_z']
    # make sure mesh2d_interface_z.max()>=wl.max() (is clipped to wl again in next step)
    if zvals_interface[-1] < uds_eta.max():
        zvals_interface[-1] = uds_eta.max()
    zw = (uds_z0+zvals_interface).clip(min=uds_bl, max=uds_eta)
       
    # correction: set interfaces below bed to nan (keeping the interface at the bed with shift)
    bool_belowbed = zw.shift({dimn_interfaces:-1}) <= uds_bl
    zw = zw.where(~bool_belowbed)
    # correction: set interfaces above wl to nan (keeping the interface at the wl with shift)
    bool_abovewl = zw.shift({dimn_interfaces:1}) >= uds_eta
    zw = zw.where(~bool_abovewl)
    
    # correction: set centers to values in between interfaces (and reshape+rename from int to cent dim)
    zcc = zw.rolling({dimn_interfaces:2}).mean()
    zcc = zcc.isel({dimn_interfaces:slice(1,None)}).rename({dimn_interfaces:dimn_layer})
    
    uds[f'{gridname}_flowelem_zw'] = zw
    uds[f'{gridname}_flowelem_zcc'] = zcc
    uds = uds.set_coords([f'{gridname}_flowelem_zw',f'{gridname}_flowelem_zcc'])
    return uds


def reconstruct_zw_zcc_fromzsigma(uds):
    """
    reconstruct full grid output (time/face-varying z-values) for zsigmavalue model without full grid output. Implemented in https://issuetracker.deltares.nl/browse/UNST-5477
    based on https://cfconventions.org/cf-conventions/cf-conventions.html#_ocean_sigma_over_z_coordinate
    """
    
    dimn_layer, dimn_interfaces = get_vertical_dimensions(uds)
    gridname = uds.grid.name
    
    # temporarily decode default fillvalues
    # TODO: xarray only decodes explicitly set fillvalues: https://github.com/Deltares/dfm_tools/issues/490
    uds = xu.UgridDataset(decode_default_fillvals(uds.obj),grids=uds.grids)
    
    osz_formulaterms_int_dict = get_formula_terms(uds,varn_contains='interface')
    # osz_formulaterms_lay_dict = get_formula_terms(uds,varn_contains='layer')
    
    uds_eta = uds[osz_formulaterms_int_dict['eta']] #mesh2d_s1
    uds_depth = uds[osz_formulaterms_int_dict['depth']] #mesh2d_bldepth: positive version of mesh2d_flowelem_bl, but this one is always in file
    uds_depth_c = uds[osz_formulaterms_int_dict['depth_c']] #mesh2d_sigmazdepth
    uds_zlev_int = uds[osz_formulaterms_int_dict['zlev']] #mesh2d_interface_z
    uds_sigma_int = uds[osz_formulaterms_int_dict['sigma']] #mesh2d_interface_sigma
    # uds_zlev_lay = uds[osz_formulaterms_lay_dict['zlev']] #mesh2d_layer_z
    # uds_sigma_lay = uds[osz_formulaterms_lay_dict['sigma']] #mesh2d_layer_sigma
    
    # for levels k where sigma(k) has a defined value and zlev(k) is not defined:
    # z(n,k,j,i) = eta(n,j,i) + sigma(k)*(min(depth_c,depth(j,i))+eta(n,j,i))
    zw_sigmapart = uds_eta + uds_sigma_int*(uds_depth.clip(max=uds_depth_c)+uds_eta)
    # zcc_sigmapart = uds_eta + uds_sigma_lay*(uds_depth.clip(max=uds_depth_c)+uds_eta)
    # for levels k where zlev(k) has a defined value and sigma(k) is not defined: 
    # z(n,k,j,i) = zlev(k)
    zw_zpart = uds_zlev_int.clip(min=-uds_depth) #added clipping of zvalues with bedlevel
    # zcc_zpart = uds_zlev_lay.clip(min=-uds_depth) #added clipping of zvalues with bedlevel
    zw = zw_sigmapart.fillna(zw_zpart)
    # zcc = zcc_sigmapart.fillna(zcc_zpart)
    
    # correction: set interfaces below bed to nan (keeping the interface at the bed with shift)
    bool_belowbed = zw.shift({dimn_interfaces:-1}) <= -uds_depth
    zw = zw.where(~bool_belowbed)
    
    # correction: set centers to values in between interfaces (and reshape+rename from int to cent dim)
    zcc = zw.rolling({dimn_interfaces:2}).mean()
    zcc = zcc.isel({dimn_interfaces:slice(1,None)}).rename({dimn_interfaces:dimn_layer})
    
    uds[f'{gridname}_flowelem_zw'] = zw
    uds[f'{gridname}_flowelem_zcc'] = zcc
    uds = uds.set_coords([f'{gridname}_flowelem_zw',f'{gridname}_flowelem_zcc'])
    return uds


def reconstruct_zw_zcc(uds:xu.UgridDataset):
    """
    reconstruct full grid output (time/face-varying z-values) for all layertypes,
    calls the respective reconstruction function

    Parameters
    ----------
    uds : xu.UgridDataset
        DESCRIPTION.

    Raises
    ------
    KeyError
        DESCRIPTION.

    Returns
    -------
    uds : TYPE
        DESCRIPTION.

    """
    
    gridname = uds.grid.name
    varname_zint = f'{gridname}_flowelem_zw'
    
    #reconstruct zw/zcc variables (if not in file) and treat as fullgrid mapfile from here
    if varname_zint in uds.variables: #fullgrid info already available, so continuing
        print(f'zw/zcc (fullgrid) values already present in Dataset in variable {varname_zint}')
    elif len(uds.filter_by_attrs(standard_name='ocean_sigma_z_coordinate')) != 0:
        print('zsigma-layer model, computing zw/zcc (fullgrid) values and treat as fullgrid model from here')
        uds = reconstruct_zw_zcc_fromzsigma(uds)
    elif len(uds.filter_by_attrs(standard_name='ocean_sigma_coordinate')) != 0:
        print('sigma-layer model, computing zw/zcc (fullgrid) values and treat as fullgrid model from here')
        uds = reconstruct_zw_zcc_fromsigma(uds)
    elif f'{gridname}_layer_z' in uds.variables:
        print('z-layer model, computing zw/zcc (fullgrid) values and treat as fullgrid model from here')
        uds = reconstruct_zw_zcc_fromz(uds)
    else:
        raise KeyError(f'layers present, but unknown layertype, expected one of variables: {gridname}_flowelem_zw, {gridname}_layer_sigma, {gridname}_layer_z')
    return uds

    
def get_Dataset_atdepths(data_xr:xu.UgridDataset, depths, reference:str ='z0'):    
    """
    Lazily depth-slice a dataset with layers. Performance can be increased by using a subset of variables or subsetting the dataset in any dimension.
    This can be done for instance with ds.isel(time=-1) or uds.ugrid.sel(x=slice(),y=slice()) to subset a ugrid dataset in space.
    The return dataset only contains the sliced variables.
    
    Parameters
    ----------
    data_xr : xu.UgridDataset
        has to be Dataset (not a DataArray), otherwise mesh2d_flowelem_zw etc are not available (interface z values)
        in case of zsigma/sigma layers (or fullgrid), it is advisable to .sel()/.isel() the time dimension first, because that is less computationally heavy
    depths : TYPE
        int/float or list/array of int/float. Depths w.r.t. reference level. If reference=='waterlevel', depth>0 returns only nans. If reference=='bedlevel', depth<0 returns only nans. Depths are sorted and only uniques are kept.
    reference : str, optional
        compute depth w.r.t. z0/waterlevel/bed. The default is 'z0'.
    zlayer_z0_selnearest : bool, optional
        Use xr.interp() to interpolate zlayer model to z-value. Only possible for reference='z' (not 'waterlevel' or 'bedlevel'). Only used if "mesh2d_layer_z" is present (zlayer model)
        This is faster but results in values interpolated between zcc (z cell centers), so it is different than slicing.. The default is False.

    Raises
    ------
    Exception
        DESCRIPTION.

    Returns
    -------
    xu.UgridDataset
        Dataset with the depth-sliced variables.

    """
    
    depth_vardimname = f'depth_from_{reference}'
    
    dimn_layer, dimn_interfaces = get_vertical_dimensions(data_xr)
    
    if dimn_layer is not None: #D-FlowFM mapfile
        gridname = data_xr.grid.name
        varname_zint = f'{gridname}_flowelem_zw'
        dimname_layc = dimn_layer
        dimname_layw = dimn_interfaces
        varname_wl = f'{gridname}_s1'
        varname_bl = f'{gridname}_flowelem_bl'
    elif 'laydim' in data_xr.dims: #D-FlowFM hisfile
        varname_zint = 'zcoordinate_w'
        dimname_layc = 'laydim'
        dimname_layw = 'laydimw'
        varname_wl = 'waterlevel'
        varname_bl = 'bedlevel'
    else:
        print(UserWarning('depth/layer dimension not found, probably 2D model, returning input Dataset')) #TODO: this can also be at depth, since slice will put parts of model dry (and allnan if below wl or below bl). Implement this
        return data_xr #early return
    
    #create depth xr.DataArray
    if isinstance(depths,(float,int)):
        depth_dims = ()
    else:
        depths = np.unique(depths) #array of unique+sorted floats/ints
        depth_dims = (depth_vardimname)
    depths_xr = xr.DataArray(depths,dims=depth_dims,attrs={'units':'m',
                                                           'reference':f'model_{reference}',
                                                           'positive':'up'}) #TODO: make more in line with CMEMS etc
    
    #potentially construct fullgrid info (zcc/zw)
    if dimn_layer is not None: #D-FlowFM mapfile
        data_xr = reconstruct_zw_zcc(data_xr)
    
    #correct reference level
    if reference=='z0':
        zw_reference = data_xr[varname_zint]
    elif reference=='waterlevel':
        if varname_wl not in data_xr.variables:
            raise KeyError(f'get_Dataset_atdepths() called with reference=waterlevel, but {varname_wl} variable not present')
        data_wl = data_xr[varname_wl]
        zw_reference = data_xr[varname_zint] - data_wl
    elif reference=='bedlevel':
        if varname_bl not in data_xr.variables:
            raise KeyError(f'get_Dataset_atdepths() called with reference=bedlevel, but {varname_bl} variable not present') #TODO: in case of zsigma/sigma it can also be -mesh2d_bldepth
        data_bl = data_xr[varname_bl]
        zw_reference = data_xr[varname_zint] - data_bl
    else:
        raise KeyError(f'unknown reference "{reference}" (options are z0, waterlevel and bedlevel)')
    
    print('>> subsetting data on fixed depth in fullgrid z-data: ',end='')
    dtstart = dt.datetime.now()
    
    #get layerbool via z-interface value (zw), check which celltop-interfaces are above/on depth and which which cellbottom-interfaces are below/on depth
    bool_topinterface_abovedepth = zw_reference.isel({dimname_layw:slice(1,None)}) >= depths_xr
    bool_botinterface_belowdepth = zw_reference.isel({dimname_layw:slice(None,-1)}) <= depths_xr
    bool_topbotinterface_arounddepth = bool_topinterface_abovedepth & bool_botinterface_belowdepth #this bool also automatically excludes all values below bed and above wl
    bool_topbotinterface_arounddepth = bool_topbotinterface_arounddepth.rename({dimname_layw:dimname_layc}) #correct dimname for interfaces to centers
    
    #subset variables that have no, time, face and/or layer dims, slice only variables with all three dims (and add to subset)
    bool_dims = [x for x in bool_topbotinterface_arounddepth.dims if x!=depth_vardimname] #exclude depth_vardimname (present if multiple depths supplied), since it is not present in pre-slice variables
    variables_toslice = [var for var in data_xr.data_vars if set(bool_dims).issubset(data_xr[var].dims)]
    
    #actual slicing with .where().max()
    ds_atdepths = data_xr[variables_toslice].where(bool_topbotinterface_arounddepth).max(dim=dimname_layc,keep_attrs=True) #set all layers but one to nan, followed by an arbitrary reduce (max in this case) #TODO: check if attributes should be passed/altered
    #TODO: suppress warning (upon plotting/load/etc): "C:\Users\veenstra\Anaconda3\envs\dfm_tools_env\lib\site-packages\dask\array\reductions.py:640: RuntimeWarning: All-NaN slice encountered" >> already merged in xarray: https://github.com/dask/dask/pull/9916
    ds_atdepths = ds_atdepths.drop_dims([dimname_layw,dimname_layc],errors='ignore') #dropping interface dim if it exists, since it does not correspond to new depths dim
    
    #add depth as coordinate var
    ds_atdepths[depth_vardimname] = depths_xr
    ds_atdepths = ds_atdepths.set_coords([depth_vardimname])
    
    print(f'{(dt.datetime.now()-dtstart).total_seconds():.2f} sec')
    
    return ds_atdepths


def rasterize_ugrid(uds:xu.UgridDataset, ds_like:xr.Dataset = None, resolution:float = None):
    """
    Rasterizing ugrid dataset to regular dataset. ds_like has higher priority than `resolution`. If both are not passed, a raster is generated of at least 200x200
    inspired by xugrid.plot.imshow and xugrid.ugrid.ugrid2d.rasterize/rasterize_like.
    
    Parameters
    ----------
    uds : xu.UgridDataset
        DESCRIPTION.
    ds_like : xr.Dataset, optional
        xr.Dataset with ed x and y variables to interpolate uds to. The default is None.
    resolution : float, optional
        Only used if ds_like is not supplied. The default is None.
    
    Raises
    ------
    Exception
        DESCRIPTION.

    Returns
    -------
    ds : TYPE
        DESCRIPTION.

    """
    
    if isinstance(uds,xu.core.wrap.UgridDataset):
        xu_facedim = uds.grid.face_dimension
        uds_facevars = Dataset_varswithdim(uds,xu_facedim)
        face_str = f'Dataset with {len(uds_facevars.data_vars)} face variables'
    elif isinstance(uds,xu.core.wrap.UgridDataArray):
        face_str = 'DataArray'
    else:
        raise TypeError(f'rasterize_ugrid expected xu.core.wrap.UgridDataset or xu.core.wrap.UgridDataArray, got {type(uds)} instead')
    
    if ds_like is None:
        xmin, ymin, xmax, ymax = uds.grid.bounds
        dx = xmax - xmin
        dy = ymax - ymin
        if resolution is None: # check if a rasterization resolution is passed, otherwise default to 200 raster cells otherwise for the smallest axis.
            resolution = min(dx, dy) / 200
        d = abs(resolution)
        regx = np.arange(xmin + 0.5 * d, xmax, d)
        regy = np.arange(ymin + 0.5 * d, ymax, d)
        ds_like = xr.DataArray(np.empty((len(regy), len(regx))), {"y": regy, "x": regx}, ["y", "x"])
    
    print(f'>> rasterizing ugrid {face_str} to shape=({len(ds_like.y)},{len(ds_like.x)}): ',end='')
    dtstart = dt.datetime.now()
    ds = uds.ugrid.rasterize_like(other=ds_like)
    print(f'{(dt.datetime.now()-dtstart).total_seconds():.2f} sec')
    
    return ds


def plot_ztdata(data_xr_sel, varname, ax=None, only_contour=False, **kwargs):
    """
    

    Parameters
    ----------
    data_xr : TYPE
        DESCRIPTION.
    varname : TYPE
        DESCRIPTION.
    ax : matplotlib.axes._subplots.AxesSubplot, optional
        the figure axis. The default is None.
    only_contour : bool, optional
        Wheter to plot contour lines of the dataset. The default is False.
    kwargs : TYPE
        properties to give on to the pcolormesh function.
    
    Raises
    ------
    Exception
        DESCRIPTION.

    Returns
    -------
    pc : matplotlib.collections.QuadMesh
        DESCRIPTION.
    
    """
    
    if not ax: ax=plt.gca()
    
    if len(data_xr_sel[varname].shape) != 2:
        raise ValueError(f'unexpected number of dimensions in requested squeezed variable ({data_xr_sel[varname].shape}), first use data_xr.isel(stations=int) to select a single station') #TODO: can also have a different cause, improve message/testing?
    
    #repair zvalues at wl/wl (filling nans and clipping to wl/bl). bfill replaces nan values with last valid value, this is necessary to enable pcolormesh to work. clip forces data to be within bl/wl
    #TODO: put clip in preproces_hisnc to make plotting easier?
    data_xr_sel['zcoordinate_c'] = data_xr_sel['zcoordinate_c'].bfill(dim='laydim').clip(min=data_xr_sel.bedlevel,max=data_xr_sel.waterlevel)
    data_xr_sel['zcoordinate_w'] = data_xr_sel['zcoordinate_w'].bfill(dim='laydimw').clip(min=data_xr_sel.bedlevel,max=data_xr_sel.waterlevel)
    
    # generate 2 2d grids for the x & y bounds (you can also give one 2D array as input in case of eg time varying z coordinates)
    data_fromhis_zcor = data_xr_sel['zcoordinate_w'].to_numpy() 
    data_fromhis_zcor = np.concatenate([data_fromhis_zcor,data_fromhis_zcor[[-1],:]],axis=0)
    time_np = data_xr_sel.time.to_numpy()
    time_cor = np.concatenate([time_np,time_np[[-1]]])
    time_mesh_cor = np.tile(time_cor,(data_fromhis_zcor.shape[-1],1)).T
    if only_contour:
        pc = data_xr_sel[varname].plot.contour(ax=ax, x='time', y='zcoordinate_c', **kwargs)
    else:
        #pc = data_xr_sel[varname].plot.pcolormesh(ax=ax, x='time', y='zcoordinate_w', **kwargs) #TODO: not possible to put center values on interfaces, so more difficult approach needed
        pc = ax.pcolormesh(time_mesh_cor, data_fromhis_zcor, data_xr_sel[varname], **kwargs)
   
    return pc

