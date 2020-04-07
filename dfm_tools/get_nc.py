# -*- coding: utf-8 -*-
"""
GNU GENERAL PUBLIC LICENSE
	      Version 3, 29 June 2007

dfm_tools are post-processing tools for Delft3D FM
Copyright (C) 2020 Deltares

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.


INFORMATION
This script is part of dfm_tools: https://github.com/openearth/dfm_tools
Check the README.rst on github for other available functions
Check the tests folder on github for example scripts (this is the dfm_tools pytest testbank)
Check the pptx and example figures in (created by the testbank): N:/Deltabox/Bulletin/veenstra/info dfm_tools

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
    
    import warnings
    import numpy as np
    import datetime as dt
    import pandas as pd
    from netCDF4 import Dataset
    
    from dfm_tools.get_nc_helpers import get_ncfilelist, get_ncvardimlist, get_ncvarobject, get_timesfromnc, get_timeid_fromdatetime, get_hisstationlist, get_stationid_fromstationlist, ghostcell_filter, get_varname_fromnc
    
    #get variable info
    nc_varobject = get_ncvarobject(file_nc, varname)
    data_nc = Dataset(file_nc)
    
    #get list of station dimnames
    vars_pd, dims_pd = get_ncvardimlist(file_nc=file_nc)
    vars_pd_stats = vars_pd[vars_pd['dtype']=='|S1']
    dimname_stat_validvals = []
    for iR, vars_pd_stat in vars_pd_stats.iterrows():
        [dimname_stat_validvals.append(x) for x in vars_pd_stat['dimensions']]
    
    
    listtype_int = [int, np.int, np.int8, np.int16, np.int32, np.int64]
    listtype_str = [str]
    listtype_range = [list, range, np.ndarray]
    listtype_datetime = [dt.datetime, np.datetime64]
    listtype_daterange = [pd.DatetimeIndex]
    
    #CHECK IF VARNAME IS STATION NAMES (STRINGS), OFFER ALTERNATIVE RETRIEVAL METHOD
    if nc_varobject.dtype == '|S1':
        warnings.warn('variable "%s" should be retrieved with separate function:\nfrom dfm_tools.get_nc_helpers import get_hisstationlist; station_names = get_hisstationlist(file_nc=file_nc, varname="%s") (or use any varname there to retrieve corresponding station list)'%(varname,varname))
    
    #TIMES CHECKS
    dimn_time = get_varname_fromnc(data_nc,'time')
    if dimn_time not in nc_varobject.dimensions: #dimension time is not available in variable
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
        elif type(timestep) in listtype_range:
            if type(timestep[0]) in listtype_int:
                time_ids = timestep
            elif type(timestep[0]) in listtype_datetime:
                time_ids = get_timeid_fromdatetime(data_nc_datetimes_pd, timestep)
            else:
                raise Exception('ERROR: timestep variable type is list/range/ndarray (%s), but type of timestep[0] not anticipated (%s), options:\n - int\n - np.int64\n - datetime\n - np.datetime64'%(type(timestep),type(timestep[0])))
        elif type(timestep) in listtype_daterange:
            time_ids = get_timeid_fromdatetime(data_nc_datetimes_pd, timestep)
        elif type(timestep) in listtype_int:
            time_ids = [timestep]
        elif type(timestep) in listtype_datetime:
            time_ids = get_timeid_fromdatetime(data_nc_datetimes_pd, [timestep])
        else:
            raise Exception('ERROR: timestep variable type not anticipated (%s), options:\n - datetime/int\n - list/range/ndarray of datetime/int\n - pandas daterange\n - "all"'%(type(timestep)))
        #check if requested times are within range of netcdf
        if np.max(time_ids) > len(data_nc_datetimes_pd)-1:
            raise Exception('ERROR: requested maximum timestep (%d) is larger than available in netcdf file (%d)'%(np.max(time_ids),len(data_nc_datetimes_pd)-1))
        if np.min(time_ids) < -len(data_nc_datetimes_pd):
            raise Exception('ERROR: requested minimum timestep (%d) is smaller than available in netcdf file (%d)'%(np.min(time_ids),-len(data_nc_datetimes_pd)))
    
    #LAYER CHECKS
    dimn_layer = get_varname_fromnc(data_nc,'nmesh2d_layer')
    if dimn_layer not in nc_varobject.dimensions: #no layer dimension in model and/or variable
        if layer is not None:
            raise Exception('ERROR: netcdf variable (%s) does not contain layers, but parameter layer is provided'%(varname))
    else: #layers are present in variable
        dimn_layer_id = nc_varobject.dimensions.index(dimn_layer)
        nlayers = nc_varobject.shape[dimn_layer_id]
        if layer is None:
            raise Exception('ERROR: netcdf variable contains a layer dimension, but parameter layer not provided (can be "all")')
        #convert layer to list of int if it is not already
        if layer is str('all'):
            layer_ids = range(nlayers)
        elif type(layer) in listtype_range:
            if type(layer[0]) in listtype_int:
                layer_ids = np.unique(layer)
            else:
                raise Exception('ERROR: timestep lay type not anticipated (%s), (list/range/ndarray of) int are accepted (or "all")'%(type(layer)))            
        elif type(layer) in listtype_int:
            layer_ids = [layer]
        else:
            raise Exception('ERROR: timestep lay type not anticipated (%s), (list/range/ndarray of) int are accepted (or "all")'%(type(layer)))
        #check if requested layers are within range of netcdf
        if np.max(layer_ids) > nlayers-1:
            raise Exception('ERROR: requested max layer (%d) is larger than available in netcdf file (%d)'%(np.max(layer_ids),nlayers-1))
        if np.min(layer_ids) < -nlayers:
            raise Exception('ERROR: requested min layer (%d) is larger than available in netcdf file (%d)'%(np.min(layer_ids),-nlayers))
    
    #DEPTH CHECKS
    if depth is not None:
        raise Exception('ERROR: depth argument is provided, but this is not implemented yet')
    
    #STATION/GENERAL_STRUCTURES CHECKS
    dimname_stat_validvals_boolpresent = [x in nc_varobject.dimensions for x in dimname_stat_validvals]
    if not any(dimname_stat_validvals_boolpresent):
        if station is not None:
            raise Exception('ERROR: netcdf file variable (%s) does not contain stations/general_structures, but parameter station is provided'%(varname))
    else: #stations are present
        #get appropriate station list
        dimname_stat_validvals_id = np.where(dimname_stat_validvals_boolpresent)[0][0]
        station_name_list_pd = get_hisstationlist(file_nc,varname=varname)
        if station is None:
            raise Exception('ERROR: netcdf variable contains a station/general_structures dimension, but parameter station not provided (can be "all"), available stations/crs/generalstructures:\n%s\nretrieve entire station list:\nfrom dfm_tools.get_nc_helpers import get_hisstationlist; get_hisstationlist(file_nc,varname="%s")'%(station_name_list_pd, varname))
        #convert station to list of int if it is not already
        if station is str('all'):
            station_ids = range(len(station_name_list_pd))
        elif type(station) in listtype_range:
            if type(station[0]) in listtype_int:
                station_ids = station
            elif type(station[0]) in listtype_str:
                station_ids = get_stationid_fromstationlist(station_name_list_pd, station, varname)
            else:
                raise Exception('ERROR1: station variable type not anticipated (%s), (list/range/ndarray of) strings or ints are accepted (or "all")'%(type(station)))
        elif type(station) in listtype_int:
            station_ids = [station]
        elif type(station) in listtype_str:
            station_ids = get_stationid_fromstationlist(station_name_list_pd, [station], varname)
        else:
            raise Exception('ERROR2: station variable type not anticipated (%s), (list/range/ndarray of) strings or ints are accepted (or "all")'%(type(station)))
        #check if requested times are within range of netcdf
        if np.max(station_ids) > len(station_name_list_pd)-1:
            raise Exception('ERROR: requested highest station id (%d) is larger than available in netcdf file (%d)'%(np.max(station_ids),len(station_name_list_pd)-1))
        if np.min(station_ids) < -len(station_name_list_pd):
            raise Exception('ERROR: requested lowest station id (%d) is smaller than available in netcdf file (%d)'%(np.min(station_ids),-len(station_name_list_pd)))
    
    
    #check faces existence, variable could have ghost cells if partitioned
    dimn_faces = get_varname_fromnc(data_nc,'mesh2d_nFaces')
    dimn_nodes = get_varname_fromnc(data_nc,'mesh2d_nNodes')
    dimn_edges = get_varname_fromnc(data_nc,'nmesh2d_edge')
    dimn_nFlowElem = get_varname_fromnc(data_nc,'nFlowElem')
    dimn_nFlowLink = get_varname_fromnc(data_nc,'nFlowLink')
    
    file_ncs = get_ncfilelist(file_nc, multipart)
    
    for iF, file_nc_sel in enumerate(file_ncs):
        if len(file_ncs) > 1:
            print('processing mapdata from domain %04d of %04d'%(iF, len(file_ncs)-1))
        
        nc_varobject_sel = get_ncvarobject(file_nc_sel, varname)

        concat_axis = 0 #default value, overwritten by faces dimension
        values_selid = []
        values_dimlens = [] #list(nc_values.shape)
        for iD, nc_values_dimsel in enumerate(nc_varobject_sel.dimensions):
            if nc_values_dimsel in [dimn_faces, dimn_nFlowElem]: # domain-like variable is present, so there are multiple domains
                nonghost_ids = ghostcell_filter(file_nc_sel)
                if nonghost_ids is None:
                    values_selid.append(range(nc_varobject_sel.shape[iD]))
                else:
                    values_selid.append(nonghost_ids)
                values_dimlens.append(0) #because concatenate axis
                concat_axis = iD
            elif nc_values_dimsel in [dimn_nodes, dimn_edges, dimn_nFlowLink]: # domain-like variable is present, so there are multiple domains
                values_selid.append(range(nc_varobject_sel.shape[iD]))
                values_dimlens.append(0) #because concatenate axis
                concat_axis = iD
            elif nc_values_dimsel in dimname_stat_validvals:
                values_selid.append(station_ids)
                values_dimlens.append(len(station_ids))
            elif nc_values_dimsel == dimn_time:
                values_selid.append(time_ids)
                values_dimlens.append(len(time_ids))
            elif nc_values_dimsel == dimn_layer:
                values_selid.append(layer_ids)
                values_dimlens.append(len(layer_ids))
            else:
                #warnings.warn('not a predefined dimension name')
                values_selid.append(range(nc_varobject_sel.shape[iD]))
                values_dimlens.append(nc_varobject_sel.shape[iD])

        if len(file_ncs) > 1:
            #initialize array
            if iF == 0:
                values_all = np.ma.empty(values_dimlens)
            #concatenate array
            values_all = np.ma.concatenate([values_all,nc_varobject_sel[values_selid]],axis=concat_axis)
        else:
            values_all = nc_varobject_sel[values_selid]
            
    #add metadata
    values_all.var_varname = varname
    values_all.var_object = nc_varobject #this is the netcdf variable, so contains properties like shape/units/dimensions
    if dimn_time in nc_varobject.dimensions:
        values_all.var_times = data_nc_datetimes_pd.iloc[time_ids]
    else:
        values_all.var_times = None
    if dimn_layer in nc_varobject.dimensions:
        values_all.var_layers = layer_ids
    else:
        values_all.var_layers = None
    if any(dimname_stat_validvals_boolpresent):
        values_all.var_stations = station_name_list_pd.iloc[station_ids]
    else:
        values_all.var_stations = None
    return values_all






def get_xzcoords_onintersection(file_nc, line_array=None, intersect_gridnos=None, intersect_coords=None, timestep=None, multipart=None, calcdist_fromlatlon=None):
    import warnings
    import numpy as np
    from netCDF4 import Dataset
    
    from dfm_tools.testutils import try_importmodule
    try_importmodule(modulename='shapely')
    from shapely.geometry import LineString, Point

    from dfm_tools.get_nc_helpers import get_varname_fromnc
    
    warnings.warn('WARNING: the function dfm_tools.get_nc.get_xzcoords_onintersection() will be improved, input variables and outputformat might change in the future')

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
    
    varn_mesh2d_s1 = get_varname_fromnc(data_nc,'mesh2d_s1')
    data_frommap_wl3 = get_ncmodeldata(file_nc, varname=varn_mesh2d_s1, timestep=timestep, multipart=multipart)
    data_frommap_wl3_sel = data_frommap_wl3[0,intersect_gridnos]
    varn_mesh2d_flowelem_bl = get_varname_fromnc(data_nc,'mesh2d_flowelem_bl')
    data_frommap_bl = get_ncmodeldata(file_nc, varname=varn_mesh2d_flowelem_bl, multipart=multipart)
    data_frommap_bl_sel = data_frommap_bl[intersect_gridnos]
    
    #nlay = data_frommap.shape[2]
    #nlay = data_nc.variables[varname].shape[2]
    dimn_layer = get_varname_fromnc(data_nc,'nmesh2d_layer')
    if dimn_layer is None: #no layers, 2D model
        nlay = 1
    else:
        nlay = data_nc.dimensions[dimn_layer].size
    
    varn_layer_z = get_varname_fromnc(data_nc,'mesh2d_layer_z')
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
        if verts.shape[:-2] != values.shape:
            raise Exception('size of first two dimensions of verts and dimensions of values is not equal, cannot plot')
    
    #convert to 3D
    if len(verts.shape) == 4 and verts.shape[-2] == 4 and verts.shape[-1] == 2: #from regular grid
        # flatten first two dimensions to one
        verts_3D = verts.reshape(-1,verts.shape[2],verts.shape[3])
        if not values is None:
            values_3D = values.reshape(-1)
        else:
            values_3D = None
    elif len(verts.shape) == 3 and verts.shape[-1] == 2: #from ugrid
        verts_3D = verts
        values_3D = values
    else:
        raise Exception('dimensions should be [m,n,4,2] or [cells,maxcorners,2], last dimension is xy')


    if not ax: ax=plt.gca()
    pc = matplotlib.collections.PolyCollection(verts_3D, **kwargs)
    pc.set_array(values_3D)
    ax.add_collection(pc)
    ax.autoscale()
    
    return pc









