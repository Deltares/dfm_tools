# -*- coding: utf-8 -*-
"""
dfm_tools are post-processing tools for Delft3D FM
Copyright (C) 2020 Deltares. All rights reserved.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

All names, logos, and references to "Deltares" are registered trademarks of
Stichting Deltares and remain full property of Stichting Deltares at all times.
All rights reserved.


INFORMATION
This script is part of dfm_tools: https://github.com/openearth/dfm_tools
Check the README.rst on github for other available functions
Check the tests folder on github for example scripts (this is the dfm_tools pytest testbank)
Check the pptx and example figures in (created by the testbank): N:/Deltabox/Bulletin/veenstra/info dfm_tools

Created on Fri Feb 14 12:45:11 2020

@author: veenstra
"""


def get_ncmodeldata(file_nc, varname=None, timestep=None, layer=None, depth=None, station=None, multipart=None, get_linkedgridinfo=False, silent=False):
    """

    Parameters
    ----------
    file_nc : str
        path to netcdf file.
    varname : str, optional
        string of netcdf variable name (key/standard_name only?).
    timestep : TYPE, optional
        (list/range/ndarray of) 0-based int or datetime. Can be used to select one or more specific timesteps, or 'all'. The default is None.
    layer : TYPE, optional
        (list/range/ndarray of) 0-based int. The default is None.
    depth : TYPE, optional
        DESCRIPTION. The default is None.
    station : TYPE, optional
        DESCRIPTION. The default is None.
    multipart : TYPE, optional
        set to False if you want only one of the map domains, can be left out otherwise. The default is None.

    Raises
    ------
    Exception
        DESCRIPTION.

    Returns
    -------
    values_all : TYPE
        DESCRIPTION.

    """

    import warnings
    import numpy as np
    import datetime as dt
    import pandas as pd
    from netCDF4 import Dataset

    from dfm_tools.get_nc_helpers import get_ncfilelist, get_ncvardimlist, get_varnamefrom_keyslongstandardname, get_variable_timevardim, get_timesfromnc, get_timeid_fromdatetime, get_hisstationlist, get_stationid_fromstationlist, ghostcell_filter, get_varname_fromnc

    #get variable info (also checks if varname exists in keys, standard name, long name)
    data_nc = Dataset(file_nc)
    varname = get_varnamefrom_keyslongstandardname(file_nc, varname) #get varname from varkeys/standardname/longname if exists
    nc_varobject = data_nc.variables[varname]

    #get list of station dimnames
    vars_pd, dims_pd = get_ncvardimlist(file_nc=file_nc)

    listtype_int = [int, np.int8, np.int16, np.int32, np.int64]
    listtype_str = [str]
    listtype_range = [list, range, np.ndarray, pd.RangeIndex]
    listtype_datetime = [dt.datetime, np.datetime64]
    listtype_daterange = [pd.DatetimeIndex]

    #CHECK IF VARNAME IS STATION NAMES (STRINGS), OFFER ALTERNATIVE RETRIEVAL METHOD
    if nc_varobject.dtype == '|S1':
        print('variable "%s" should probably be retrieved with separate function:\nfrom dfm_tools.get_nc_helpers import get_hisstationlist\nstation_names = get_hisstationlist(file_nc=file_nc, varname="%s") (or use any varname there to retrieve corresponding station list)'%(varname,varname))
    if 'time' in varname.lower():
        print('variable "%s" should probably be retrieved with separate function:\nfrom dfm_tools.get_nc_helpers import get_timesfromnc\ntimes = get_timesfromnc(file_nc=file_nc, varname="%s")'%(varname, varname))


    #TIMES CHECKS
    #dimn_time = get_varname_fromnc(data_nc,'time',vardim='dim')
    #varn_time = get_varname_fromnc(data_nc,'time',vardim='var')
    #if dimn_time is None: #dimension with a name close to 'time' is not available in variable, try to get time dimension from 'time' variable
    #    try:
    #        dimn_time = data_nc.variables[varn_time].dimensions[0]
    #    except:
    #        print('using dimn_time as variable to get dimn_time failed')
    varn_time, dimn_time = get_variable_timevardim(file_nc=file_nc, varname=varname)
    if dimn_time not in nc_varobject.dimensions: #dimension time is not available in variable
        if timestep is not None:
            raise Exception('ERROR: netcdf file variable (%s) does not contain times, but parameter timestep is provided'%(varname))
    else: #time dimension is present
        data_nc_timevar = data_nc.variables[varn_time]
        time_length = data_nc_timevar.shape[0]
        data_nc_datetimes_pd = get_timesfromnc(file_nc, varname=varname, retrieve_ids=[0,-1], silent=silent) #get selection of times
        if timestep is None:
            raise Exception('ERROR: netcdf variable contains a time dimension, but parameter timestep not provided (can be "all"), first and last timestep:\n%s\nretrieve entire times list:\nfrom dfm_tools.get_nc_helpers import get_timesfromnc\ntimes_pd = get_timesfromnc(file_nc=file_nc, varname="%s")'%(pd.DataFrame(data_nc_datetimes_pd),varname))
        #convert timestep to list of int if it is not already
        if timestep is str('all'):
            data_nc_datetimes_pd = get_timesfromnc(file_nc, varname=varname, silent=silent) #get all times
            time_ids = range(len(data_nc_datetimes_pd))
        elif type(timestep) in listtype_range:
            if len(timestep) == 0:
                raise Exception('ERROR: timestep variable type is list/range/ndarray (%s), but it has no length'%(type(timestep)))
            elif type(timestep[0]) in listtype_int:
                data_nc_datetimes_pd = get_timesfromnc(file_nc, varname=varname, retrieve_ids=timestep, silent=silent) #get selection of times
                time_ids = timestep
            elif type(timestep[0]) in listtype_datetime:
                data_nc_datetimes_pd = get_timesfromnc(file_nc, varname=varname, silent=silent) #get all times
                time_ids = get_timeid_fromdatetime(data_nc_datetimes_pd, timestep)
                data_nc_datetimes_pd = data_nc_datetimes_pd.loc[time_ids] #get selection of times
            else:
                raise Exception('ERROR: timestep variable type is list/range/ndarray (%s), but type of timestep[0] not anticipated (%s), options:\n - int\n - np.int64\n - datetime\n - np.datetime64'%(type(timestep),type(timestep[0])))
        elif type(timestep) in listtype_daterange:
            data_nc_datetimes_pd = get_timesfromnc(file_nc, varname=varname, silent=silent) #get all times
            time_ids = get_timeid_fromdatetime(data_nc_datetimes_pd, timestep)
            data_nc_datetimes_pd = data_nc_datetimes_pd.loc[time_ids] #get selection of times
        elif type(timestep) in listtype_int:
            data_nc_datetimes_pd = get_timesfromnc(file_nc, varname=varname, retrieve_ids=[timestep], silent=silent) #get selection of times
            time_ids = [timestep]
        elif type(timestep) in listtype_datetime:
            data_nc_datetimes_pd = get_timesfromnc(file_nc, varname=varname, silent=silent) #get all times
            time_ids = get_timeid_fromdatetime(data_nc_datetimes_pd, [timestep])
            data_nc_datetimes_pd = data_nc_datetimes_pd.loc[time_ids] #get selection of times
        else:
            raise Exception('ERROR: timestep variable type not anticipated (%s), options:\n - datetime/int\n - list/range/ndarray of datetime/int\n - pandas daterange\n - "all"'%(type(timestep)))
        #convert to positive index, make unique(+sort), convert to list because of indexing with np.array of len 1 errors sometimes
        time_ids = list(np.unique(np.array(range(time_length))[time_ids]))
        #check if requested times are within range of netcdf
        if np.max(time_ids) > time_length-1:
            raise Exception('ERROR: requested maximum timestep (%d) is larger than available in netcdf file (%d)'%(np.max(time_ids),time_length-1))

    #LAYER CHECKS
    dimn_layer = get_varname_fromnc(data_nc,'nmesh2d_layer',vardim='dim')
    if dimn_layer not in nc_varobject.dimensions: #no layer dimension in model and/or variable
        if layer is not None:
            raise Exception('ERROR: netcdf variable (%s) does not contain layers, but argument layer is provided'%(varname))
    else: #layers are present in variable
        dimn_layer_id = nc_varobject.dimensions.index(dimn_layer)
        nlayers = nc_varobject.shape[dimn_layer_id]
        if layer is None:
            raise Exception('ERROR: netcdf variable contains a layer dimension, but argument layer not provided (can be "all")\nnumber of layers: %d (numbered 0 to %d)'%(nlayers, nlayers-1))
        #convert layer to list of int if it is not already
        if layer is str('all') or layer is str('top') or layer is str('bottom'):
            layer_ids = range(nlayers)
        elif type(layer) in listtype_range:
            if type(layer[0]) in listtype_int:
                layer_ids = np.unique(layer)
            else:
                raise Exception('ERROR: layer variable type not anticipated (%s), (list/range/ndarray of) int are accepted (or "all")'%(type(layer)))
        elif type(layer) in listtype_int:
            layer_ids = [layer]
        else:
            raise Exception('ERROR: layer variable type not anticipated (%s), (list/range/ndarray of) int are accepted (or "all", "top" or "bottom")'%(type(layer)))
        #convert to positive index, make unique(+sort), convert to list because of indexing with np.array of len 1 errors sometimes
        layer_ids = list(np.unique(np.array(range(nlayers))[layer_ids]))
        #check if requested layers are within range of netcdf
        if np.max(layer_ids) > nlayers-1:
            raise Exception('ERROR: requested max layer (%d) is larger than available in netcdf file (%d)'%(np.max(layer_ids),nlayers-1))

    #DEPTH CHECKS
    if depth is not None:
        raise Exception('ERROR: depth argument is provided, but vertical slicing is not implemented yet, try layer argument instead')

    #STATION/GENERAL_STRUCTURES CHECKS
    vars_pd_stats = vars_pd[(vars_pd['dtype']=='|S1') & (vars_pd['dimensions'].apply(lambda x: dimn_time not in x))]
    dimname_stat_validvals = []
    for iR, vars_pd_stat in vars_pd_stats.iterrows():
        dimname_stat_validvals.append(vars_pd_stat['dimensions'][0]) #only append first dimension, the other one is often 'name_len'
    dimname_stat_validvals_boolpresent = [x in nc_varobject.dimensions for x in dimname_stat_validvals]
    if not any(dimname_stat_validvals_boolpresent):
        if station is not None:
            raise Exception('ERROR: netcdf file variable (%s) does not contain stations/general_structures, but argument station is provided'%(varname))
    else: #stations are present
        #get appropriate station list
        station_name_list_pd = get_hisstationlist(file_nc,varname=varname)
        if station is None:
            raise Exception('ERROR: netcdf variable contains a station/general_structures dimension, but argument station not provided (can be "all"), available stations/crs/generalstructures:\n%s\nretrieve entire station list:\nfrom dfm_tools.get_nc_helpers import get_hisstationlist\nstations_pd = get_hisstationlist(file_nc,varname="%s")'%(station_name_list_pd, varname))
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
        #convert to positive index, make unique(+sort), convert to list because of indexing with np.array of len 1 errors sometimes
        station_ids = list(np.unique(np.array(range(len(station_name_list_pd)))[station_ids]))
        #check if requested times are within range of netcdf
        if np.max(station_ids) > len(station_name_list_pd)-1:
            raise Exception('ERROR: requested highest station id (%d) is larger than available in netcdf file (%d)'%(np.max(station_ids),len(station_name_list_pd)-1))


    #check faces existence, variable could have ghost cells if partitioned
    dimn_faces = get_varname_fromnc(data_nc,'mesh2d_nFaces',vardim='dim')
    dimn_nodes = get_varname_fromnc(data_nc,'mesh2d_nNodes',vardim='dim')
    dimn_edges = get_varname_fromnc(data_nc,'nmesh2d_edge',vardim='dim')
    dimn_nFlowElem = get_varname_fromnc(data_nc,'nFlowElem',vardim='dim')
    dimn_nFlowLink = get_varname_fromnc(data_nc,'nFlowLink',vardim='dim')

    #revert back to single partition if non-partitioned variable is requested
    bool_varpartitioned = any([True for x in nc_varobject.dimensions if x in [dimn_faces, dimn_nodes, dimn_edges, dimn_nFlowElem, dimn_nFlowLink]])
    if not bool_varpartitioned:
        multipart = False

    #get list of partitioned files
    file_ncs = get_ncfilelist(file_nc, multipart)

    for iF, file_nc_sel in enumerate(file_ncs):
        if (len(file_ncs) > 1) and not silent:
            print('processing mapdata from domain %04d of %04d'%(iF, len(file_ncs)-1))

        data_nc_sel = Dataset(file_nc_sel)
        nc_varobject_sel = data_nc_sel.variables[varname]
        
        concat_axis = 0 #default value, overwritten by faces dimension
        ghost_removeids = [] #default value, overwritten by faces dimension

        values_selid = []
        values_dimlens = [] #list(nc_values.shape)
        values_dimlinkedgrid = [] #list(nc_values.shape)
        try:
            nc_varobject_sel_coords = nc_varobject_sel.coordinates
        except:
            nc_varobject_sel_coords = None
        if not silent:
            print('varname: %s  %s  %s, coordinates=(%s)'%(varname, nc_varobject_sel.shape, nc_varobject_sel.dimensions, nc_varobject_sel_coords))

        if len(nc_varobject_sel.dimensions) == 0:
            raise Exception('variable contains no dimensions, cannot retrieve values')

        for iD, nc_values_dimsel in enumerate(nc_varobject_sel.dimensions):
            if nc_values_dimsel in [dimn_faces, dimn_nFlowElem]: # domain-like variable is present, so there are multiple domains (with ghost cells)
                nonghost_ids = ghostcell_filter(file_nc_sel)
                if nonghost_ids is not None:
                    ghost_removeids = np.where(~nonghost_ids)[0] #remove after retrieval, since that is faster than retrieving nonghost ids
                values_selid.append(range(nc_varobject_sel.shape[iD]))
                values_dimlens.append(0) #because concatenate axis
                concat_axis = iD
            elif nc_values_dimsel in [dimn_nodes, dimn_edges, dimn_nFlowLink]: # domain-like variable is present, so there are multiple domains (no ghost cells)
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

            #get info about grid variables related to varname
            if get_linkedgridinfo and (nc_values_dimsel not in [dimn_time,dimn_layer]+dimname_stat_validvals):
                vars_pd_relevant = vars_pd[(vars_pd['ndims']<=2) & (vars_pd['dimensions'].apply(lambda x: nc_values_dimsel in x)) & -(vars_pd['dimensions'].apply(lambda x: dimn_time in x))]
                values_dimlinkedgrid.append(vars_pd_relevant)

                print('\tlinkedvars for dimension "%s":'%(nc_values_dimsel))
                #print('nc_varobject_sel.dimensions: %s'%([nc_varobject_sel.dimensions]))
                for iLV, linkedvar in vars_pd_relevant.iterrows():
                    print('\t\t%s  %s  %s'%(linkedvar['nc_varkeys'], linkedvar['shape'], linkedvar['dimensions']))
                #print('nc_values_dimsel: %s'%(nc_values_dimsel))
                #print('vars_pd_relevant:\n%s'%(vars_pd_relevant))
            else:
                values_dimlinkedgrid.append(None)

        #get selected data (including ghostcells because that is faster)
        nc_varobject_sel_selids_raw = nc_varobject_sel[values_selid]

        #remove ghost cells (cannot delete from masked array, so delete from array and mask and then couple again)
        if ghost_removeids is not []:
            nc_varobject_sel_selids = np.delete(nc_varobject_sel_selids_raw,ghost_removeids,axis=concat_axis)
            if nc_varobject_sel_selids_raw.mask.any() != False:
                nc_varobject_sel_selids_mask = np.delete(nc_varobject_sel_selids_raw.mask,ghost_removeids,axis=concat_axis)
                nc_varobject_sel_selids.mask = nc_varobject_sel_selids_mask

        #concatenate to other partitions
        if len(file_ncs) > 1:
            #initialize array
            if iF == 0:
                values_all = np.ma.empty(values_dimlens)
                values_all[:] = np.nan
            #concatenate array
            values_all = np.ma.concatenate([values_all, nc_varobject_sel_selids], axis=concat_axis)
        else:
            values_all = nc_varobject_sel_selids
        data_nc_sel.close()

    #optional extraction of top/bottom layer, convenient for z-layer models since top and/or bottom layers are often masked for part of the cells
    if layer is str('top') or layer is str('bottom'):
        warnings.warn('you are retrieving data from the %s valid layer of each cell. it is assumed that the last axis of the variable is the layer axis'%(layer))
        if not values_all.mask.any(): #if (all values in) the mask are False
            raise Exception('there is no mask present in this dataset (or all its values are False), use layer=[0,-1] to get the bottom and top layers')
        layerdim_id = nc_varobject.dimensions.index(dimn_layer)
        if layer is str('top'):
            bottomtoplay = values_all.shape[layerdim_id]-1-(~np.flip(values_all.mask,axis=layerdim_id)).argmax(axis=layerdim_id) #get index of first False value from the flipped array (over layer axis) and correct with size of that dimension. this corresponds to the top layer of each cell in case of D-Flow FM
        if layer is str('bottom'):
            bottomtoplay = (~values_all.mask).argmax(axis=layerdim_id) #get index of first False value from the original array
        values_selid_topbot = []
        for iD, dimlen in enumerate(values_all.shape):
            if iD == layerdim_id:
                values_selid_topbot.append(bottomtoplay)
            elif iD == concat_axis and not '_his.nc' in file_nc: #his files have no partitions and thus no concat_axis, this forces to 'else' and to transpose (no testcase available)
                values_selid_topbot.append(np.array(range(dimlen)))
            else:
                values_selid_topbot.append(np.array([range(dimlen)]).T)
        values_all_topbot = values_all[tuple(values_selid_topbot)] #layer dimension is removed due to advanced indexing instead of slicing
        values_all_topbot = np.expand_dims(values_all_topbot, axis=layerdim_id) #re-add layer dimension to dataset on original location
        values_all = values_all_topbot


    #add metadata
    values_all.var_filename = file_nc
    values_all.var_varname = varname
    values_all.var_dimensions = nc_varobject.dimensions
    values_all.var_shape = nc_varobject.shape
    values_all.var_dtype = nc_varobject.dtype
    values_all.var_linkedgridinfo = values_dimlinkedgrid
    #values_all.var_ncobject = data_nc #this is the netcdf object retrieved with netCDF4.Dataset() #disabled, since it becomes invalid after closing the dataset
    values_all.var_ncvarobject = f"from netCDF4 import Dataset;data_nc = Dataset('{file_nc}');nc_varobject = data_nc.variables['{varname}'];print(nc_varobject)" # nc_varobject #this is the netcdf variable, contains properties like shape/units/dimensions #disabled, since it becomes invalid after closing the dataset
    values_all.var_ncattrs = nc_varobject.__dict__ #values in nc_varobject.ncattrs() or hasattr(nc_varobject,'attributename')

    if dimn_time in nc_varobject.dimensions:
        values_all.var_times = data_nc_datetimes_pd
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
        
    data_nc.close()
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
        raise Exception('ERROR: argument timestep not provided, this is necessary to retrieve correct waterlevel or fullgrid output')

    print('calculating distance for all crossed cells, from first point of line (should not take long, but if it does, optimisation is needed)')
    nlinecoords = line_array.shape[0]
    nlinedims = len(line_array.shape)
    ncrosscellparts = intersect_coords.shape[0]
    if nlinecoords<2 or nlinedims != 2:
        raise Exception('ERROR: line_array should at least contain two xy points [[x,y],[x,y]]')

    crs_xstart = intersect_coords[:,0,0]
    crs_xstop = intersect_coords[:,0,1]
    crs_ystart = intersect_coords[:,1,0]
    crs_ystop = intersect_coords[:,1,1]

    #calculate distance between celledge-linepart crossing (is zero when line iL crosses cell)
    distperline_tostart = np.zeros((ncrosscellparts,nlinecoords-1))
    distperline_tostop = np.zeros((ncrosscellparts,nlinecoords-1))
    linepart_length = np.zeros((nlinecoords))
    for iL in range(nlinecoords-1):
        #calculate length of lineparts
        line_section_part = LineString(line_array[iL:iL+2,:])
        if not calcdist_fromlatlon:
            linepart_length[iL+1] = line_section_part.length
        else:
            linepart_length[iL+1] = calc_dist_latlon(line_array[iL,0],line_array[iL+1,0],line_array[iL,1],line_array[iL+1,1])

        #get distance between all lineparts and point (later used to calculate distance from beginpoint of closest linepart)
        for iP in range(ncrosscellparts):
            distperline_tostart[iP,iL] = line_section_part.distance(Point(crs_xstart[iP],crs_ystart[iP]))
            distperline_tostop[iP,iL] = line_section_part.distance(Point(crs_xstop[iP],crs_ystop[iP]))
    linepart_lengthcum = np.cumsum(linepart_length)
    cross_points_closestlineid = np.argmin(np.maximum(distperline_tostart,distperline_tostop),axis=1)
    print('finished calculating distance for all crossed cells, from first point of line')

    data_nc = Dataset(file_nc)

    #nlay = data_frommap.shape[2]
    #nlay = data_nc.variables[varname].shape[2]
    dimn_layer = get_varname_fromnc(data_nc,'nmesh2d_layer',vardim='dim')
    if dimn_layer is None: #no layers, 2D model
        nlay = 1
    else:
        nlay = data_nc.dimensions[dimn_layer].size

    varn_layer_z = get_varname_fromnc(data_nc,'mesh2d_layer_z',vardim='var')
    varn_layer_fullgrid = get_varname_fromnc(data_nc,'mesh2d_flowelem_zw',vardim='var')
    if varn_layer_fullgrid is not None:
        laytyp = 'fullgrid'
        zvals_interface_allfaces = get_ncmodeldata(file_nc, varname=varn_layer_fullgrid, timestep=timestep, multipart=multipart)
        zvals_interface = zvals_interface_allfaces[0,intersect_gridnos,:].T #transpose to make in line with 2D sigma dataset
    else:
        varn_mesh2d_s1 = get_varname_fromnc(data_nc,'mesh2d_s1',vardim='var')
        data_frommap_wl3 = get_ncmodeldata(file_nc, varname=varn_mesh2d_s1, timestep=timestep, multipart=multipart)
        data_frommap_wl3_sel = data_frommap_wl3[0,intersect_gridnos]
        varn_mesh2d_flowelem_bl = get_varname_fromnc(data_nc,'mesh2d_flowelem_bl',vardim='var')
        data_frommap_bl = get_ncmodeldata(file_nc, varname=varn_mesh2d_flowelem_bl, multipart=multipart)
        data_frommap_bl_sel = data_frommap_bl[intersect_gridnos]
        
        if varn_layer_z is not None:
            laytyp = 'zlayer'
            #zvals_cen = get_ncmodeldata(file_nc=file_map, varname='mesh2d_layer_z', lay='all')#, multipart=False)
            #zvals_interface = get_ncmodeldata(file_nc=file_map, varname='mesh2d_interface_z')#, multipart=False)
            zvals_interface = data_nc.variables['mesh2d_interface_z'][:]
        else: #sigma or 2D model
            laytyp = 'sigmalayer'
            #zvals_cen = np.linspace(data_frommap_bl_sel,data_frommap_wl3_sel,nlay)
            zvals_interface = np.linspace(data_frommap_bl_sel,data_frommap_wl3_sel,nlay+1)

        
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
        if laytyp in ['sigmalayer','fullgrid']:
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
    data_nc.close()
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
        print('analyzing netdata from domain %04d of %04d (counting max number of facenodes)'%(iF, len(file_ncs)-1))
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
                edge_verts_all = np.ma.empty((0,4,edge_verts.shape[2])) #create edge verts, which will contain the two edge node coordinates, as well as the two center coordinates from neighbouring faces
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
    import numpy as np
    
    if not values is None:
        #squeeze values (remove dimensions with length 1)
        values = np.squeeze(values)
        #check if data shape is equal
        if verts.shape[:-2] != values.shape:
            raise Exception('size of first dimensions of verts (%s) and dimensions of squeezed values (%s) is not equal, cannot plot. Flatten your values array or if the values are on cell edges, try providing ugrid_all.edge_verts instead'%(verts.shape[:-2],values.shape))

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








def plot_background(ax=None, projection=None, google_style='satellite', resolution=1, features=None, nticks=6, latlon_format=False, gridlines=False, **kwargs):
    """
    this definition uses cartopy to plot a geoaxis and a satellite basemap and coastlines. A faster alternative for a basemap is contextily:
    import contextily as ctx
    fig, ax = plt.subplots(1,1)
    ctx.add_basemap(ax, source=ctx.providers.Esri.WorldImagery, crs="EPSG:28992")
    More info at: https://contextily.readthedocs.io/en/latest/reference.html

    Parameters
    ----------
    ax : cartopy.mpl.geoaxes.GeoAxesSubplot, optional
        DESCRIPTION. The default is None.
    projection : integer, cartopy._crs.CRS or cartopy._epsg._EPSGProjection, optional
        DESCRIPTION. The default is None.
    google_style : Nonetype or string, optional
       The style of the Google Maps tiles. One of None, ‘street’, ‘satellite’, ‘terrain’, and ‘only_streets’. The default is 'satellite'.
    resolution : int, optional
        resolution for the Google Maps tiles. 1 works wel for global images, 12 works well for a scale of Grevelingen lake, using 12 on global scale will give you a server timeout. The default is 1.
    features : string, optional
        Features to plot, options: None, 'ocean', 'rivers', 'land', 'countries', 'countries_highres', 'coastlines', 'coastlines_highres'. The default is None.
    nticks : TYPE, optional
        DESCRIPTION. The default is 6.
    latlon_format : bool, optional
        DESCRIPTION. The default is False.
    gridlines : TYPE, optional
        DESCRIPTION. The default is False.
    **kwargs : TYPE
        additional arguments for ax.add_feature or ax.coastlines(). examples arguments and values are: alpha=0.5, facecolor='none', edgecolor='gray', linewidth=0.5, linestyle=':'

    Raises
    ------
    Exception
        DESCRIPTION.

    Returns
    -------
    ax : TYPE
        DESCRIPTION.

    """

    import matplotlib.pyplot as plt
    import numpy as np

    from dfm_tools.testutils import try_importmodule
    try_importmodule(modulename='cartopy') #check if cartopy was installed since it is an optional module, also happens in plot_cartopybasemap()

    import cartopy
    import cartopy.crs as ccrs
    import cartopy.io.img_tiles as cimgt
    import cartopy.feature as cfeature
    import cartopy.mpl.ticker as cticker

    dummy = ccrs.epsg(28992) #to make cartopy realize it has a cartopy._epsg._EPSGProjection class (maybe gets fixed with cartopy updates, see unittest test_cartopy_epsg)
    if ax is None: #provide axis projection on initialisation, cannot be edited later on
        if projection is None:
            projection=ccrs.PlateCarree() #projection of cimgt.GoogleTiles, useful default
        elif isinstance(projection, (cartopy._epsg._EPSGProjection, cartopy._crs.CRS)): #checks if argument is an EPSG projection or CRS projection (like PlateCarree, Mercator etc)
            pass
        elif type(projection) is int:
            projection = ccrs.epsg(projection)
        else:
            raise Exception('argument projection should be of type integer, cartopy._crs.CRS or cartopy._epsg._EPSGProjection')
        fig, ax = plt.subplots(subplot_kw={'projection': projection})
        #ax = plt.axes(projection=projection)
    elif type(ax) is cartopy.mpl.geoaxes.GeoAxesSubplot:
        if projection is not None:
            print('arguments ax and projection are both provided, the projection from the ax is used so the projection argument is ignored')
    else:
        raise Exception('argument ax should be of type cartopy.mpl.geoaxes.GeoAxesSubplot, leave argument empty or create correct instance with:\nimport cartopy.crs as ccrs\nfig, (ax1,ax2) = plt.subplots(1,2,figsize=(10,5), subplot_kw={"projection": ccrs.epsg(28992)})')



    if gridlines:
        ax.gridlines(draw_labels=True)
    elif nticks is not None: #only look at nticks if gridlines are not used
        extent = ax.get_extent()
        ax.set_xticks(np.linspace(extent[0],extent[1],nticks))
        ax.set_yticks(np.linspace(extent[2],extent[3],nticks))


    if google_style is not None:
        request = cimgt.GoogleTiles(style=google_style)
        ax.add_image(request,resolution)


    if features is not None:
        if type(features) is str:
            features = [features]
        elif type(features) is not list:
            raise Exception('argument features should be of type list of str')

        valid_featurelist = ['ocean','rivers','land','countries','countries_highres','coastlines','coastlines_highres']
        invalid_featurelist = [x for x in features if x not in valid_featurelist]
        if invalid_featurelist != []:
            raise Exception('invalid features %s requested, possible are: %s'%(invalid_featurelist, valid_featurelist))

        if 'ocean' in features:
            #feat = cfeature.NaturalEarthFeature(category='physical', name='ocean', facecolor=cfeature.COLORS['water'], scale='10m', edgecolor='face', alpha=alpha)
            #ax.add_feature(feat)
            ax.add_feature(cfeature.OCEAN, **kwargs)
        if 'rivers' in features:
            ax.add_feature(cfeature.RIVERS, **kwargs)
        if 'land' in features:
            #feat = cfeature.NaturalEarthFeature(category='physical', name='land', facecolor=cfeature.COLORS['land'], scale='10m', edgecolor='face', alpha=alpha)
            #ax.add_feature(feat)
            ax.add_feature(cfeature.LAND, **kwargs)
        if 'countries' in features:
            ax.add_feature(cfeature.BORDERS, **kwargs)
        if 'countries_highres' in features:
            feat = cfeature.NaturalEarthFeature(category='cultural', name='admin_0_countries', scale='10m')
            ax.add_feature(feat, **kwargs)
        if 'coastlines' in features:
            ax.add_feature(cfeature.COASTLINE, **kwargs)
        if 'coastlines_highres' in features:
            ax.coastlines(resolution='10m', **kwargs)

    if latlon_format:
        lon_formatter = cticker.LongitudeFormatter()
        lat_formatter = cticker.LatitudeFormatter()
        ax.xaxis.set_major_formatter(lon_formatter)
        ax.yaxis.set_major_formatter(lat_formatter)


    return ax








def plot_ztdata(file_nc, dfmtools_hisvar, statid_subset=0, ax=None, mask_data=True, only_contour=False, **kwargs):
    """


    Parameters
    ----------
    dfmtools_hisvar : numpy.ma.core.MaskedArray
        dfm_tools.get_nc.get_ncmodeldata output structure, which is of type numpy.ma.core.MaskedArray and has several extra properties attached which which it is possible to retrieve the correct z interface data.
    statid_subset : int, optional
        the station id for the dfmtools_hisvar, so 0-based since it retrieves from a numpy array. beware that get_ncmodeldata sorts the stations you request, so use an index like dfmtools_hisvar.var_stations['station_id'].tolist().index(stat). Avoid this issue by only retrieving data for one station. The default is 0.
    ax : matplotlib.axes._subplots.AxesSubplot, optional
        the figure axis. The default is None.
    mask_data : bool, optional
        whether to repair z_interface coordinates and mask data in inactive layers. The default is True.
    only_contour : bool, optional
        Wheter to plot contour lines of the dataset. The default is False.
    **kwargs : TYPE
        properties to give on to the pcolormesh function.

    Returns
    -------
    pc : matplotlib.collections.QuadMesh
        DESCRIPTION.

    """
    import warnings
    import numpy as np
    import matplotlib.pyplot as plt

    warnings.warn('WARNING: layers in dfowfm hisfile are currently incorrect, check your figures carefully')
    #get the original ncobject and the retrieved indices from dfmtools_hisvar, used to retrieve a corresponding z array.
    time_ids = dfmtools_hisvar.var_times.index.tolist()
    stat_ids = dfmtools_hisvar.var_stations.index.tolist()
    layer_ids = dfmtools_hisvar.var_layers
    layer_ids_interf = layer_ids+[layer_ids[-1]+1]
    #data_fromhis_zcen = nc_obj.variables['zcoordinate_c'][time_ids,stat_ids,layer_ids]
    #data_fromhis_zcor = nc_obj.variables['zcoordinate_w'][time_ids,stat_ids,layer_ids_interf]
    #data_fromhis_wl = nc_obj.variables['waterlevel'][time_ids,stat_ids]
    data_fromhis_zcen = get_ncmodeldata(file_nc=file_nc, varname='zcoordinate_c', timestep=time_ids, station=stat_ids, layer=layer_ids)
    data_fromhis_zcor = get_ncmodeldata(file_nc=file_nc, varname='zcoordinate_w', timestep=time_ids, station=stat_ids)[:,:,layer_ids_interf] #dfmtools cannot find layers in this variable
    data_fromhis_wl = get_ncmodeldata(file_nc=file_nc, varname='waterlevel', timestep=time_ids, station=stat_ids)

    #remove station dimension
    if len(data_fromhis_zcor.shape) == 3:
        data_fromhis_zcen_flat = data_fromhis_zcen[:,statid_subset,:]
        data_fromhis_zcor_flat = data_fromhis_zcor[:,statid_subset,:]
        dfmtools_hisvar_flat = dfmtools_hisvar[:,statid_subset,:]
    elif len(data_fromhis_zcor.shape) == 2:
        data_fromhis_zcen_flat = data_fromhis_zcen[:,statid_subset]
        data_fromhis_zcor_flat = data_fromhis_zcor[:,statid_subset]
        dfmtools_hisvar_flat = dfmtools_hisvar[:,statid_subset]
    else:
        raise Exception('unexpected number of dimensions')
    data_fromhis_wl_flat = data_fromhis_wl[:,statid_subset]

    """
    fig,(ax1)=plt.subplots()
    ax1.plot(np.arange(data_fromhis_zcen_flat.shape[1])+.5,data_fromhis_zcen_flat[5,:],'x',label='centers')
    ax1.plot(range(data_fromhis_zcor_flat.shape[1]),data_fromhis_zcor_flat[5,:],'o',label='corners')
    ax1.legend()
    """
    if mask_data:
        bool_zcen_equaltop = (data_fromhis_zcen_flat==data_fromhis_zcen_flat[:,-1:]).all(axis=0)
        id_zcentop = np.argmax(bool_zcen_equaltop) # id of first z_center that is equal to z_center of last layer
        if (data_fromhis_zcor_flat[:,id_zcentop] > data_fromhis_zcen_flat[:,id_zcentop]).any():
            print('correcting z interface values')
            data_fromhis_zcor_flat[:,id_zcentop+1] = data_fromhis_wl_flat
            data_fromhis_zcor_flat[:,id_zcentop] = (data_fromhis_zcen_flat[:,id_zcentop-1]+data_fromhis_zcen_flat[:,id_zcentop])/2
        bool_zcen_equaltop[id_zcentop] = False
        #bool_zcor_equaltop = (data_fromhis_zcor_flat[:,1:]==data_fromhis_zcor_flat[:,-1:]).all(axis=0)
        mask_array = np.tile(bool_zcen_equaltop,(data_fromhis_zcor_flat.shape[0],1))
        dfmtools_hisvar_flat.mask = mask_array

    if not ax: ax=plt.gca()

    # generate 2 2d grids for the x & y bounds (you can also give one 2D array as input in case of eg time varying z coordinates)
    time_mesh_cor = np.tile(dfmtools_hisvar.var_times,(data_fromhis_zcor_flat.shape[-1],1)).T
    time_mesh_cen = np.tile(dfmtools_hisvar.var_times,(data_fromhis_zcen_flat.shape[-1],1)).T
    if only_contour:
        pc = ax.contour(time_mesh_cen,data_fromhis_zcen_flat,dfmtools_hisvar_flat, **kwargs)
    else:
        #pc = ax.pcolormesh(time_mesh_cor, data_fromhis_zcor_flat, dfmtools_hisvar_flat, **kwargs)
        pc = ax.pcolor(time_mesh_cor, data_fromhis_zcor_flat, dfmtools_hisvar_flat, **kwargs) #pcolor also supports missing/masked xy data, but is slower

    return pc
