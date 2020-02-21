# -*- coding: utf-8 -*-
"""
Created on Fri Feb 14 12:45:11 2020

@author: veenstra
"""


def get_ncmodeldata(file_nc, varname, timestep=None, layer=None, depth=None, station=None, multipart=None):
    """
    file_nc: path to netcdf file
    varname: string of netcdf variable name (standard_name?)
    timestep: (list/range/ndarray of) 0-based int or datetime. Can be used to select one or more specific timesteps, or 'all'
    layer: (list/range/ndarray of) 0-based int
    depth
    station
    multipart: set to False if you want only one of the map domains, can be left out otherwise
    """
    
    import numpy as np
    import datetime as dt
    from netCDF4 import Dataset
    
    from dfm_tools.get_nc_helpers import get_ncfilelist, get_ncvardims, get_timesfromnc, get_timeid_fromdatetime, get_hisstationlist, get_stationid_fromstationlist, ghostcell_filter, get_varname_mapnc
    
    #get variable and dimension info
    #nc_varkeys, nc_dimkeys, nc_values, nc_values_shape, nc_values_dims = get_ncvardims(file_nc, varname)
    dummy, dummy, dummy, nc_values_shape, nc_values_dims = get_ncvardims(file_nc, varname)
    data_nc = Dataset(file_nc)
    
    #CHECK VARNAME, IS THERE A SEPARATE DEFINITION TO RETRIEVE DATA?
    if varname == 'station_name' or varname == 'general_structure_id' or varname == 'cross_section_name':
        raise Exception('ERROR: variable "%s" should be retrieved with separate function:\nstation_names = get_hisstationlist(file_nc=file_nc, varname_stat="%s")'%(varname,varname))

    #TIMES CHECKS
    dimn_time = get_varname_mapnc(data_nc,'time')
    if dimn_time not in nc_values_dims: #dimension time is not available in variable
        if timestep is not None:
            raise Exception('ERROR: netcdf file variable (%s) does not contain times, but parameter timestep is provided'%(varname))
    else: #time is first dimension
        if timestep is None:
            raise Exception('ERROR: netcdf variable contains a time dimension, but parameter timestep not provided (can be "all")')
        #convert timestep to list of int if it is not already
        #get times
        data_nc_datetimes_pd = get_timesfromnc(file_nc)
        if timestep is str('all'):
            time_ids = range(len(data_nc_datetimes_pd))
        elif type(timestep)==list or type(timestep)==range or type(timestep)==type(np.arange(1,2,0.5)):
            if type(timestep[0])==int: #list/range/ndarray of int
                time_ids = timestep
            elif type(timestep[0])==type(dt.datetime(1,1,1)) or type(timestep[0])==type(np.datetime64(year=1900,month=1,day=1)): #list/range/ndarray of datetime
                time_ids = get_timeid_fromdatetime(data_nc_datetimes_pd, timestep)
            else:
                raise Exception('ERROR: 1timestep variable type not anticipated (%s), (list/range/ndarray of) datetime/int are accepted (or "all")'%(type(timestep)))
        elif type(timestep)==int:
            time_ids = [timestep]
        elif type(timestep)==type(dt.datetime(1,1,1)):
            time_ids = get_timeid_fromdatetime(data_nc_datetimes_pd, [timestep])
        else:
            raise Exception('ERROR: 2timestep variable type not anticipated (%s), (list/range/ndarray of) datetime/int are accepted (or "all")'%(type(timestep)))
        #check if requested times are within range of netcdf
        if np.min(time_ids) < 0:
            raise Exception('ERROR: requested start timestep (%d) is negative'%(np.min(time_ids)))
        if np.max(time_ids) > len(data_nc_datetimes_pd)-1:
            raise Exception('ERROR: requested end timestep (%d) is larger than available in netcdf file (%d)'%(np.max(time_ids),len(data_nc_datetimes_pd)-1))
    
    #LAYER CHECKS
    dimn_layer = get_varname_mapnc(data_nc,'nmesh2d_layer')
    if dimn_layer not in nc_values_dims: #no layer dimension in model and/or variable
        if layer is not None:
            raise Exception('ERROR: netcdf variable (%s) does not contain layers, but parameter layer is provided'%(varname))
    else: #layers are present in variable
        dimn_layer_id = nc_values_dims.index(dimn_layer)
        nlayers = nc_values_shape[dimn_layer_id]
        if layer is None:
            raise Exception('ERROR: netcdf variable contains a layer dimension, but parameter layer not provided (can be "all")')
        #convert layer to list of int if it is not already
        if layer is str('all'):
            layer_ids = range(nlayers)
        elif type(layer)==list or type(layer)==range or type(layer)==type(np.arange(1,2,0.5)):
            if type(layer[0])==int: #list/range/ndarray of int
                layer_ids = np.unique(layer)
            else:
                raise Exception('ERROR: timestep lay type not anticipated (%s), (list/range/ndarray of) int are accepted (or "all")'%(type(layer)))            
        elif type(layer)==int:
            layer_ids = [layer]
        else:
            raise Exception('ERROR: timestep lay type not anticipated (%s), (list/range/ndarray of) int are accepted (or "all")'%(type(layer)))
        #check if requested layers are within range of netcdf
        if np.min(layer_ids) < 0:
            raise Exception('ERROR: requested minimal layer (%d) is negative'%(np.min(layer_ids)))
        if np.max(layer_ids) > nlayers-1:
            raise Exception('ERROR: requested max layer (%d) is larger than available in netcdf file (%d)'%(np.max(layer_ids),nlayers-1))
    
    #DEPTH CHECKS
    if depth is not None:
        raise Exception('ERROR: depth argument is provided, but this is not implemented yet')
    
    #STATION/GENERAL_STRUCTURES CHECKS
    if 'stations' not in nc_values_dims and 'general_structures' not in nc_values_dims and 'cross_section' not in nc_values_dims: #no station dimension
        if station is not None:
            raise Exception('ERROR: netcdf file variable (%s) does not contain stations/general_structures, but parameter station is provided'%(varname))
    else: #stations are present
        if station is None:
            raise Exception('ERROR: netcdf variable contains a station/general_structures dimension, but parameter station not provided (can be "all")')
        if 'stations' in nc_values_dims:
            station_name_list_pd = get_hisstationlist(file_nc) #get stations
        elif 'general_structures' in nc_values_dims:
            station_name_list_pd = get_hisstationlist(file_nc,varname_stat='general_structure_id') #get stations            
        elif 'cross_section' in nc_values_dims:
            station_name_list_pd = get_hisstationlist(file_nc,varname_stat='cross_section_name') #get stations            
        #convert station to list of int if it is not already
        if station is str('all'):
            station_ids = range(len(station_name_list_pd))
        elif type(station)==list or type(station)==range or type(station)==type(np.arange(1,2,0.5)):
            if type(station[0])==int: #list/range/ndarray of int
                station_ids = station
            elif type(station[0])==str: #list/range/ndarray of str
                station_ids = get_stationid_fromstationlist(station_name_list_pd, station)
            else:
                raise Exception('ERROR1: station variable type not anticipated (%s), (list/range/ndarray of) strings or ints are accepted (or "all")'%(type(station)))
        elif type(station)==int:
            station_ids = [station]
        elif type(station)==str:
            station_ids = get_stationid_fromstationlist(station_name_list_pd, [station])
        else:
            raise Exception('ERROR2: station variable type not anticipated (%s), (list/range/ndarray of) strings or ints are accepted (or "all")'%(type(station)))
        #check if requested times are within range of netcdf
        if np.min(station_ids) < 0:
            raise Exception('ERROR: requested lowest station id (%d) is negative'%(np.min(station_ids)))
        if np.max(station_ids) > len(station_name_list_pd)-1:
            raise Exception('ERROR: requested highest station id (%d) is larger than available in netcdf file (%d)'%(np.max(station_ids),len(station_name_list_pd)-1))
     
    
    #check faces existence, variable could have ghost cells if partitioned
    dimn_faces = get_varname_mapnc(data_nc,'mesh2d_nFaces')

    file_ncs = get_ncfilelist(file_nc, multipart)
    
    
    for iF, file_nc_sel in enumerate(file_ncs):
        if len(file_ncs) > 1:
            print('processing mapdata from domain %04d of %04d'%(iF, len(file_ncs)-1))
        
        nc_varkeys, nc_dimkeys, nc_values, nc_values_shape, nc_values_dims = get_ncvardims(file_nc_sel, varname)
        
        concat_axis = 0 #default value, overwritten by faces/stations dimension
        values_selid = []
        values_dimlens = [] #list(nc_values.shape)
        for iD, nc_values_dimsel in enumerate(nc_values_dims):
            if nc_values_dimsel == dimn_faces: # domain variable is present, so there are multiple domains
                nonghost_ids = ghostcell_filter(file_nc_sel)
                values_selid.append(nonghost_ids)
                values_dimlens.append(0) #because concatenate axis
                concat_axis = iD
            elif nc_values_dimsel == 'stations' or nc_values_dims[iD] == 'general_structures' or nc_values_dims[iD] == 'cross_section':
                values_selid.append(station_ids)
                values_dimlens.append(0) #because concatenate axis
                concat_axis = iD
            elif nc_values_dimsel == dimn_time:
                values_selid.append(time_ids)
                values_dimlens.append(len(time_ids))
            elif nc_values_dimsel == dimn_layer:
                values_selid.append(layer_ids)
                values_dimlens.append(len(layer_ids))
            else:
                print('WARNING: not a predefined dimension name')
                values_selid.append(range(nc_values.shape[iD]))
                values_dimlens.append(nc_values.shape[iD])
        
        #initialize array
        if iF == 0:
            values_all = np.ma.empty(values_dimlens)
        #concatenate array
        values_all = np.ma.concatenate([values_all,nc_values[values_selid]],axis=concat_axis)
        
        #add metadata
        values_all.var_varname = varname
        values_all.var_dimensionnames = nc_values_dims
        if dimn_time in nc_values_dims: #only faces/stations dimensions, no times or layers
            values_all.var_times = data_nc_datetimes_pd.iloc[time_ids]
        else:
            values_all.var_times = None
        if dimn_layer in nc_values_dims: #only time and faces/stations dimensions, no layers
            values_all.var_layers = layer_ids
        else:
            values_all.var_layers = None
        if 'stations' in nc_values_dims or 'general_structures' in nc_values_dims or 'cross_section' in nc_values_dims:
            values_all.var_stations = station_name_list_pd.iloc[station_ids]
        else:
            values_all.var_stations = None
        #if :
        #    values_all.stations = ...
        #else:
        #    values_all.stations = None
    return values_all






def get_xzcoords_onintersection(file_nc, line_array=None, intersect_gridnos=None, intersect_coords=None, timestep=None, multipart=None, calcdist_fromlatlon=None):
    import numpy as np
    from netCDF4 import Dataset
    try:
        from shapely.geometry import LineString, Point
    except:
        raise Exception('ERROR: cannot execute import shapely.geometry, check known bugs on https://github.com/openearth/dfm_tools for a solution')

    from dfm_tools.get_nc_helpers import get_varname_mapnc
    
    def calc_dist(x1,x2,y1,y2):
        distance = np.sqrt((x2 - x1)**2 + (y2 - y1)**2)
        return distance

    def calc_dist_latlon(x1,x2,y1,y2):
        #https://gis.stackexchange.com/questions/80881/what-is-unit-of-shapely-length-attribute
        #calculate distance in meters between latlon coordinates
        distance = np.arccos(np.sin(np.radians(y1))*np.sin(np.radians(y2))+np.cos(np.radians(y1))*np.cos(np.radians(y2))*np.cos(np.radians(x2)-np.radians(x1)))*6371000
        return distance

    #check if all necessary arguments are provided
    if file_nc is None:
        raise Exception('ERROR: argument file_nc not provided')
    if line_array is None:
        raise Exception('ERROR: argument line_array not provided')
    if intersect_gridnos is None:
        raise Exception('ERROR: argument intersect_gridnos not provided')
    if intersect_coords is None:
        raise Exception('ERROR: argument intersect_coords not provided')
    if timestep is None:
        raise Exception('ERROR: argument timestep not provided, this is necessary to retrieve correct water level')
    
    print('calculating distance for all crossed cells, from first point of line (should not take long, but if it does, optimisation is needed)')
    nlinecoords = line_array.shape[0]
    nlinedims = len(line_array.shape)
    ncrosscells = intersect_coords.shape[0]
    if nlinecoords<2 or nlinedims != 2:
        raise Exception('ERROR: line_array should at least contain two xy points [[x,y],[x,y]]')

    crs_xstart = intersect_coords[:,0,0]
    crs_xstop = intersect_coords[:,0,1]
    crs_ystart = intersect_coords[:,1,0]
    crs_ystop = intersect_coords[:,1,1]

    #calculate distance between celledge-linepart crossing (is zero when line iL crosses cell)
    distperline_tostart = np.empty((ncrosscells,nlinecoords-1))
    #distperline_tostop = np.empty((ncrosscells,nlinecoords-1))
    linepart_length = np.zeros((nlinecoords))
    for iL in range(nlinecoords-1):
        #calculate length of lineparts
        line_section_part = LineString(line_array[iL:iL+2,:])
        if not calcdist_fromlatlon:
            linepart_length[iL+1] = line_section_part.length
        else:
            linepart_length[iL+1] = calc_dist_latlon(line_array[iL,0],line_array[iL+1,0],line_array[iL,1],line_array[iL+1,1])
        
        #get distance between all lineparts and point (later used to calculate distance from beginpoint of closest linepart)
        for iP in range(ncrosscells):
            distperline_tostart[iP,iL] = line_section_part.distance(Point(crs_xstart[iP],crs_ystart[iP]))
    linepart_lengthcum = np.cumsum(linepart_length)
    cross_points_closestlineid = np.argmin(distperline_tostart,axis=1)
    print('finished calculating distance for all crossed cells, from first point of line')
    
    data_nc = Dataset(file_nc)
    
    varn_mesh2d_s1 = get_varname_mapnc(data_nc,'mesh2d_s1')
    data_frommap_wl3 = get_ncmodeldata(file_nc, varname=varn_mesh2d_s1, timestep=timestep, multipart=multipart)
    data_frommap_wl3_sel = data_frommap_wl3[0,intersect_gridnos]
    varn_mesh2d_flowelem_bl = get_varname_mapnc(data_nc,'mesh2d_flowelem_bl')
    data_frommap_bl = get_ncmodeldata(file_nc, varname=varn_mesh2d_flowelem_bl, multipart=multipart)
    data_frommap_bl_sel = data_frommap_bl[intersect_gridnos]
    
    #nlay = data_frommap.shape[2]
    #nlay = data_nc.variables[varname].shape[2]
    dimn_layer = get_varname_mapnc(data_nc,'nmesh2d_layer')
    if dimn_layer is None: #no layers, 2D model
        nlay = 1
    else:
        nlay = data_nc.dimensions[dimn_layer].size
    
    varn_layer_z = get_varname_mapnc(data_nc,'mesh2d_layer_z')
    if varn_layer_z is None:
        laytyp = 'sigmalayer'
        #zvals_cen = np.linspace(data_frommap_bl_sel,data_frommap_wl3_sel,nlay)
        zvals_interface = np.linspace(data_frommap_bl_sel,data_frommap_wl3_sel,nlay+1)
    else:
        laytyp = 'zlayer'
        #zvals_cen = get_ncmodeldata(file_nc=file_map, varname='mesh2d_layer_z', lay='all')#, multipart=False)
        #zvals_interface = get_ncmodeldata(file_nc=file_map, varname='mesh2d_interface_z')#, multipart=False)
        zvals_interface = data_nc.variables['mesh2d_interface_z'][:]
    
    #calculate distance from points to 'previous' linepoint, add lenght of previous lineparts to it
    if not calcdist_fromlatlon:
        crs_dist_starts = calc_dist(line_array[cross_points_closestlineid,0], crs_xstart, line_array[cross_points_closestlineid,1], crs_ystart) + linepart_lengthcum[cross_points_closestlineid]
        crs_dist_stops = calc_dist(line_array[cross_points_closestlineid,0], crs_xstop, line_array[cross_points_closestlineid,1], crs_ystop) + linepart_lengthcum[cross_points_closestlineid]
    else:
        crs_dist_starts = calc_dist_latlon(line_array[cross_points_closestlineid,0], crs_xstart, line_array[cross_points_closestlineid,1], crs_ystart) + linepart_lengthcum[cross_points_closestlineid]
        crs_dist_stops = calc_dist_latlon(line_array[cross_points_closestlineid,0], crs_xstop, line_array[cross_points_closestlineid,1], crs_ystop) + linepart_lengthcum[cross_points_closestlineid]
        
    crs_verts_x_all = np.empty((0,4,1))
    crs_verts_z_all = np.empty((0,4,1))
    #data_frommap_sel_flat = np.empty((0))
    for iL in range(nlay):
        zval_lay_bot = zvals_interface[iL]
        zval_lay_top = zvals_interface[iL+1]
        crs_verts_x = np.array([[crs_dist_starts,crs_dist_stops,crs_dist_stops,crs_dist_starts]]).T
        if laytyp == 'sigmalayer':
            crs_verts_z = np.array([[zval_lay_bot,zval_lay_bot,zval_lay_top,zval_lay_top]]).T
        elif laytyp == 'zlayer':
            crs_verts_z = np.repeat(np.array([[zval_lay_bot,zval_lay_bot,zval_lay_top,zval_lay_top]]).T[np.newaxis],repeats=crs_verts_x.shape[0],axis=0)
            #top z-layer is extended to water level, if wl is higher than zval_lay_top
            if iL == nlay-1:
                crs_verts_z[:,2,0] = np.maximum(zval_lay_top,data_frommap_wl3_sel.data)
                crs_verts_z[:,3,0] = np.maximum(zval_lay_top,data_frommap_wl3_sel.data)
            # zval_lay_bot lower than bedlevel should be overwritten with bedlevel
            zvalbot_belowbl_bool = crs_verts_z[:,0,0]<data_frommap_bl_sel
            crs_verts_z[zvalbot_belowbl_bool,0,0] = data_frommap_bl_sel[zvalbot_belowbl_bool]
            crs_verts_z[zvalbot_belowbl_bool,1,0] = data_frommap_bl_sel[zvalbot_belowbl_bool]
        crs_verts_x_all = np.concatenate([crs_verts_x_all, crs_verts_x])
        crs_verts_z_all = np.concatenate([crs_verts_z_all, crs_verts_z])
        #data_frommap_sel_flat = np.concatenate([data_frommap_sel_flat,data_frommap_sel[:,iL]])
    crs_verts = np.concatenate([crs_verts_x_all, crs_verts_z_all], axis=2)
    return crs_verts





def get_netdata(file_nc, multipart=None):
    import numpy as np
    
    from dfm_tools.ugrid import UGrid
    from dfm_tools.get_nc_helpers import get_ncfilelist

    file_ncs = get_ncfilelist(file_nc, multipart)
    #get all data
    num_nodes = [0]
    verts_shape2_all = []
    for iF, file_nc_sel in enumerate(file_ncs):
        print('analyzing netdata from domain %04d of %04d'%(iF, len(file_ncs)-1))
        ugrid = UGrid.fromfile(file_nc_sel)
        verts_shape2_all.append(ugrid.verts.shape[1])
    verts_shape2_max = np.max(verts_shape2_all)
        
    for iF, file_nc_sel in enumerate(file_ncs):
        print('processing netdata from domain %04d of %04d'%(iF, len(file_ncs)-1))
        #data_nc = Dataset(file_nc_sel)
        #list(data_nc.variables.keys())
        
        ugrid = UGrid.fromfile(file_nc_sel)
        node_x = ugrid.mesh2d_node_x
        node_y = ugrid.mesh2d_node_y
        node_z = ugrid.mesh2d_node_z
        faces = ugrid.mesh2d_face_nodes
        verts = ugrid.verts
        #mesh2d_edge_x = ugrid.mesh2d_edge_x
        #mesh2d_edge_y = ugrid.mesh2d_edge_y
        edge_verts = ugrid.edge_verts
        
        #setup initial array
        if iF == 0:
            node_x_all = np.ma.empty((0,))
            node_y_all = np.ma.empty((0,))
            if node_z is not None:
                node_z_all = np.ma.empty((0,))
            else:
                node_z_all = None
            verts_all = np.ma.empty((0,verts_shape2_max,verts.shape[2]))
            faces_all = np.ma.empty((0,verts_shape2_max),dtype='int32')
            #mesh2d_edge_x_all = np.ma.empty((0,))
            #mesh2d_edge_y_all = np.ma.empty((0,))
            if edge_verts is not None:
                edge_verts_all = np.ma.empty((0,2,edge_verts.shape[2]))
            else:
                edge_verts_all = None
            
        #if necessary, add masked column(s) to increase size to max in domains
        if verts.shape[1] < verts_shape2_max:
            tofew_cols = -(verts.shape[1] - verts_shape2_max)
            vcol_extra = verts[:,[0],:]
            vcol_extra.mask = True
            fcol_extra = faces[:,[0]]
            fcol_extra.mask = True
            for iC in range(tofew_cols):
                verts = np.hstack([verts,vcol_extra])
                faces = np.hstack([faces,fcol_extra])

        #merge all
        node_x_all = np.ma.concatenate([node_x_all,node_x])
        node_y_all = np.ma.concatenate([node_y_all,node_y])
        if node_z is not None:
            node_z_all = np.ma.concatenate([node_z_all,node_z])
        verts_all = np.ma.concatenate([verts_all,verts])
        faces_all = np.ma.concatenate([faces_all,faces+np.sum(num_nodes)])
        #mesh2d_edge_x_all = np.ma.concatenate([mesh2d_edge_x_all,mesh2d_edge_x])
        #mesh2d_edge_y_all = np.ma.concatenate([mesh2d_edge_y_all,mesh2d_edge_y])
        if edge_verts is not None:
            edge_verts_all = np.ma.concatenate([edge_verts_all,edge_verts])
        num_nodes.append(node_x.shape[0])

    #set all invalid values to the same value (tends to differ between partitions)
    #faces_all.data[faces_all.mask] = -999
    #faces_all.fill_value = -999
    
    ugrid_all = UGrid(node_x_all, node_y_all, faces_all, verts_all, mesh2d_node_z=node_z_all, edge_verts=edge_verts_all)
    ugrid_all
    return ugrid_all





def plot_netmapdata(verts, values=None, ax=None, **kwargs):
    #https://stackoverflow.com/questions/52202014/how-can-i-plot-2d-fem-results-using-matplotlib
    #https://stackoverflow.com/questions/49640311/matplotlib-unstructered-quadrilaterals-instead-of-triangles
    import matplotlib.pyplot as plt
    import matplotlib.collections
    
    #check if data size is equal
    if not values is None:
        if verts.shape[0] != values.shape[0]:
            raise Exception('ERROR: size of grid and values is not equal, cannot plot')
    
    if not ax: ax=plt.gca()
    pc = matplotlib.collections.PolyCollection(verts, **kwargs)
    pc.set_array(values)
    ax.add_collection(pc)
    ax.autoscale()
    return pc






