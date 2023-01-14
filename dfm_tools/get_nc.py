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
along with this program.  if not, see <http://www.gnu.org/licenses/>.

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

import warnings
import numpy as np
import datetime as dt
import pandas as pd
from netCDF4 import Dataset
import xugrid as xu
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.collections
import glob

from dfm_tools.get_nc_helpers import get_ncvarproperties, get_varnamefromattrs, get_timesfromnc, get_timeid_fromdatetime, get_hisstationlist, get_stationid_fromstationlist, ghostcell_filter, get_varname_fromnc
from dfm_tools.ugrid import UGrid
from dfm_tools.xarray_helpers import get_vertical_dimensions


def get_ncmodeldata(file_nc, varname=None, timestep=None, layer=None, station=None, multipart=None, silent=False):
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
    station : TYPE, optional
        DESCRIPTION. The default is None. Deprecated, not possible anymore (use xarray.sel instead)
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
    
    warnings.warn(DeprecationWarning('dfm_tools.get_nc.get_ncmodeldata() is deprecated, since there is an xarray alternative for multidomain FM files (xugrid).Open your file like this and use xarray sel/isel (example in postprocessing notebook):\n    data_xr_mapmerged = dfmt.open_partitioned_dataset(file_nc_map)\nFor hisfiles, use xarray with dfmt.preprocess_hisnc:\n    data_xr_his = xr.open_mfdataset(file_nc_his, preprocess=dfmt.preprocess_hisnc)'))
    
    if multipart is not None:
        warnings.warn(UserWarning('argument multipart is deprecated and will be ignored'))
    
    if isinstance(file_nc,list): #for opendap, has to support lists
        file_nc_list = file_nc
    else:
        file_nc_list = glob.glob(file_nc)
    file_nc_one = file_nc_list[0]  
    
    #get variable info (also checks if varname exists in keys, standard name, long name)
    data_nc = Dataset(file_nc_one)
    data_xr = xr.open_dataset(file_nc_one)
    varname = get_varnamefromattrs(data_xr, varname) #get varname from varkeys/standardname/longname if exists
    nc_varobject = data_nc.variables[varname]
    
    #get list of station dimnames
    vars_pd = get_ncvarproperties(file_nc=file_nc_one)
    
    listtype_int = [int, np.int8, np.int16, np.int32, np.int64]
    listtype_str = [str]
    listtype_range = [list, range, np.ndarray, pd.RangeIndex]
    listtype_datetime = [dt.datetime, np.datetime64]
    listtype_daterange = [pd.DatetimeIndex]

    #CHECK if VARNAME IS STATION NAMES (STRINGS), OFFER ALTERNATIVE RETRIEVAL METHOD
    if nc_varobject.dtype == '|S1':
        print('variable "%s" should probably be retrieved with separate function:\nfrom dfm_tools.get_nc_helpers import get_hisstationlist\nstation_names = get_hisstationlist(file_nc=file_nc, varname="%s") (or use any varname there to retrieve corresponding station list)'%(varname,varname))
    if 'time' in varname.lower():
        print('variable "%s" should probably be retrieved with separate function:\nfrom dfm_tools.get_nc_helpers import get_timesfromnc\ntimes = get_timesfromnc(file_nc=file_nc, varname="%s")'%(varname, varname))
    
    #TIMES CHECKS
    dimn_time = 'time' #hard coded, since easy to change with xarray
    if dimn_time not in nc_varobject.dimensions: #dimension time is not available in variable
        if timestep is not None:
            raise Exception('ERROR: netcdf file variable (%s) does not contain times, but parameter timestep is provided'%(varname))
    else: #time dimension is present
        data_nc_timevar = data_nc.variables['time']
        time_length = data_nc_timevar.shape[0]
        data_nc_datetimes_pd = get_timesfromnc(file_nc_one, varname=varname) #get all times
        if timestep is None:
            raise Exception('ERROR: netcdf variable contains a time dimension, but parameter timestep not provided (can be "all"), first and last timestep:\n%s\nretrieve entire times list:\nfrom dfm_tools.get_nc_helpers import get_timesfromnc\ntimes_pd = get_timesfromnc(file_nc=file_nc, varname="%s")'%(pd.DataFrame(data_nc_datetimes_pd),varname))
        #convert timestep to list of int if it is not already
        if timestep is str('all'):
            time_ids = range(len(data_nc_datetimes_pd))
        elif type(timestep) in listtype_range:
            if len(timestep) == 0:
                raise Exception('ERROR: timestep variable type is list/range/ndarray (%s), but it has no length'%(type(timestep)))
            elif type(timestep[0]) in listtype_int:
                data_nc_datetimes_pd = data_nc_datetimes_pd.iloc[timestep] #get selection of times
                time_ids = timestep
            elif type(timestep[0]) in listtype_datetime:
                time_ids = get_timeid_fromdatetime(data_nc_datetimes_pd, timestep)
                data_nc_datetimes_pd = data_nc_datetimes_pd.iloc[time_ids] #get selection of times
            else:
                raise Exception('ERROR: timestep variable type is list/range/ndarray (%s), but type of timestep[0] not anticipated (%s), options:\n - int\n - np.int64\n - datetime\n - np.datetime64'%(type(timestep),type(timestep[0])))
        elif type(timestep) in listtype_daterange:
            time_ids = get_timeid_fromdatetime(data_nc_datetimes_pd, timestep)
            data_nc_datetimes_pd = data_nc_datetimes_pd.iloc[time_ids] #get selection of times
        elif type(timestep) in listtype_int:
            time_ids = [timestep]
            data_nc_datetimes_pd = data_nc_datetimes_pd.iloc[time_ids] #get selection of times
        elif type(timestep) in listtype_datetime:
            time_ids = get_timeid_fromdatetime(data_nc_datetimes_pd, [timestep])
            data_nc_datetimes_pd = data_nc_datetimes_pd.iloc[time_ids] #get selection of times
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
    
    #STATION/GENERAL_STRUCTURES CHECKS
    vars_pd_stats = vars_pd[(vars_pd['dtype'].astype(str).str.startswith('|S') | (vars_pd['dtype']=='object')) & (vars_pd['dimensions'].apply(lambda x: dimn_time not in x))] #TODO: better check for bytes string
    dimname_stat_validvals = []
    for iR, vars_pd_stat in vars_pd_stats.iterrows():
        dimname_stat_validvals.append(vars_pd_stat['dimensions'][0]) #only append first dimension, the other one is often 'name_len'
    dimname_stat_validvals_boolpresent = [x in nc_varobject.dimensions for x in dimname_stat_validvals]
    if not any(dimname_stat_validvals_boolpresent):
        if station is not None:
            raise Exception('ERROR: netcdf file variable (%s) does not contain stations/general_structures, but argument station is provided'%(varname))
    else: #stations are present
        #get appropriate station list
        station_name_list_pd = get_hisstationlist(file_nc_one,varname=varname)
        if station is None:
            raise Exception('ERROR: netcdf variable contains a station/general_structures dimension, but argument station not provided (can be "all"), available stations/crs/generalstructures:\n%s\nretrieve entire station list:\nfrom dfm_tools.get_nc_helpers import get_hisstationlist\nstations_pd = get_hisstationlist(file_nc,varname="%s")'%(station_name_list_pd, varname))
        #convert station to list of int if it is not already
        if station is str('all'):
            station_ids = range(len(station_name_list_pd))
        elif type(station) in listtype_range:
            if type(station[0]) in listtype_int:
                station_ids = station
            elif type(station[0]) in listtype_str:
                station_ids = get_stationid_fromstationlist(station_name_list_pd, station)
            else:
                raise Exception('ERROR1: station variable type not anticipated (%s), (list/range/ndarray of) strings or ints are accepted (or "all")'%(type(station)))
        elif type(station) in listtype_int:
            station_ids = [station]
        elif type(station) in listtype_str:
            station_ids = get_stationid_fromstationlist(station_name_list_pd, [station])
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
    
    for iF, file_nc_sel in enumerate(file_nc_list):
        if (len(file_nc_list) > 1) and not silent:
            print('processing mapdata from domain %04d of %04d'%(iF, len(file_nc_list)-1))

        data_nc_sel = Dataset(file_nc_sel)
        nc_varobject_sel = data_nc_sel.variables[varname]
        
        concat_axis = 0 #default value, overwritten by faces dimension
        ghost_removeids = [] #default value, overwritten by faces/edges dimension

        values_selid = []
        values_dimlens = [] #list(nc_values.shape)
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
                nonghost_bool = ghostcell_filter(file_nc_sel)
                if nonghost_bool is not None:
                    ghost_removeids = np.where(~nonghost_bool)[0] #remove after retrieval, since that is faster than retrieving nonghost ids or using a boolean
                values_selid.append(range(nc_varobject_sel.shape[iD]))
                values_dimlens.append(0) #because concatenate axis
                concat_axis = iD
            elif nc_values_dimsel in [dimn_edges]: # domain-like variable is present, so there are multiple domains (edges from partition boundaries are removed)
                if 0:#bool_varpartitioned:
                    mesh2d_edge_faces = data_nc_sel.variables[get_varname_fromnc(data_nc_sel,'mesh2d_edge_faces',vardim='var')][:]
                    part_edges_removebool = (mesh2d_edge_faces==0).any(axis=1) #array is 1 based indexed, 0 means missing # & (np.in1d(mesh2d_edge_faces[:,0],ghost_removeids-1) | np.in1d(mesh2d_edge_faces[:,1],ghost_removeids-1))
                    part_edges_removeids = np.where(part_edges_removebool)[0]
                    ghost_removeids = part_edges_removeids #to make equal to faces varname
                values_selid.append(range(nc_varobject_sel.shape[iD]))
                values_dimlens.append(0) #because concatenate axis
                concat_axis = iD
            elif nc_values_dimsel in [dimn_nodes, dimn_nFlowLink]: # domain-like variable is present, so there are multiple domains (no ghost cells)
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

        #get selected data (including ghostcells because that is faster)
        nc_varobject_sel_selids_raw = nc_varobject_sel[values_selid]

        #remove ghost cells (cannot delete from masked array, so delete from array and mask and then couple again)
        if ghost_removeids is not []:
            nc_varobject_sel_selids = np.delete(nc_varobject_sel_selids_raw,ghost_removeids,axis=concat_axis)
            if nc_varobject_sel_selids_raw.mask.any() != False:
                nc_varobject_sel_selids_mask = np.delete(nc_varobject_sel_selids_raw.mask,ghost_removeids,axis=concat_axis)
                nc_varobject_sel_selids.mask = nc_varobject_sel_selids_mask
                
        #concatenate to other partitions
        if len(file_nc_list) > 1:
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
            bottomtoplay = values_all.shape[layerdim_id]-1-(~np.flip(values_all.mask,axis=layerdim_id)).argmax(axis=layerdim_id) #get index of first False value from the flipped array (over layer axis) and correct with size of that dimension. This corresponds to the top layer of each cell in case of D-Flow FM
        if layer is str('bottom'):
            bottomtoplay = (~values_all.mask).argmax(axis=layerdim_id) #get index of first False value from the original array. This corresponds to the top layer of each cell in case of D-Flow FM
        values_selid_topbot = []
        for iD, dimlen in enumerate(values_all.shape):
            if iD == layerdim_id:
                values_selid_topbot.append(bottomtoplay)
            elif iD == concat_axis and not '_his.nc' in file_nc_one: #his files have no partitions and thus no concat_axis, this forces to 'else' and to transpose (no testcase available)
                values_selid_topbot.append(np.array(range(dimlen)))
            else:
                values_selid_topbot.append(np.array([range(dimlen)]).T)
        values_all_topbot = values_all[tuple(values_selid_topbot)] #layer dimension is removed due to advanced indexing instead of slicing
        values_all_topbot = np.expand_dims(values_all_topbot, axis=layerdim_id) #re-add layer dimension to dataset on original location
        values_all = values_all_topbot


    #add metadata
    values_all.var_filename = file_nc_one
    values_all.var_varname = varname
    values_all.var_dimensions = nc_varobject.dimensions
    values_all.var_shape = nc_varobject.shape
    values_all.var_dtype = nc_varobject.dtype
    values_all.var_ncvarobject = f"from netCDF4 import Dataset;data_nc = Dataset('{file_nc_one}');nc_varobject = data_nc.variables['{varname}'];print(nc_varobject)" # nc_varobject #this is the netcdf variable, contains properties like shape/units/dimensions #disabled, since it becomes invalid after closing the dataset
    values_all.var_ncattrs = nc_varobject.__dict__ #values in nc_varobject.ncattrs() or hasattr(nc_varobject,'attributename')

    if dimn_time in nc_varobject.dimensions:
        values_all.var_times = data_nc_datetimes_pd
    else:
        values_all.var_times = None
    
    if dimn_layer in nc_varobject.dimensions:
        values_all.var_layers = layer_ids
    else:
        values_all.var_layers = None
    
    data_nc.close()
    return values_all


def get_ugrid_verts(data_xr_map):
    """
    getting ugrid verts from xugrid mapfile. This might be phased out in the future.
    """
    data_xr_grid = data_xr_map.ugrid.grid.to_dataset()
    face_nos = data_xr_grid.mesh2d_face_nodes.load()
    bool_nonemptyfacenode = face_nos!=-1
    facenos_nonan_min = face_nos.where(bool_nonemptyfacenode).min() #replace nans and get minval
    if facenos_nonan_min==1: #for some reason, curvedbend is 1-based indexed, grevelingen is not
        face_nos = face_nos-1
    if face_nos.dtype!='int': #for some reason, curvedbend idx is float instead of int
        face_nos = face_nos.astype(int)
    
    xu_nodedim = data_xr_map.grid.node_dimension
    
    face_nnodecoords_x = data_xr_grid.mesh2d_node_x.isel({xu_nodedim:face_nos}).where(bool_nonemptyfacenode)
    face_nnodecoords_y = data_xr_grid.mesh2d_node_y.isel({xu_nodedim:face_nos}).where(bool_nonemptyfacenode)
    ugrid_all_verts = np.c_[face_nnodecoords_x.to_numpy()[...,np.newaxis],face_nnodecoords_y.to_numpy()[...,np.newaxis]]
    return ugrid_all_verts


def calc_dist_pythagoras(x1,x2,y1,y2): # only used in dfm_tools.ugrid
    distance = np.sqrt((x2 - x1)**2 + (y2 - y1)**2)
    return distance


def calc_dist_haversine(lon1,lon2,lat1,lat2): # only used in dfm_tools.ugrid
    """
    calculates distance between lat/lon coordinates in meters
    https://community.esri.com/t5/coordinate-reference-systems-blog/distance-on-a-sphere-the-haversine-formula/ba-p/902128
    """
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
        raise Exception('nan encountered in calc_dist_latlon distance, replaced by 0') #warnings.warn
        #distance[np.isnan(distance)] = 0
    return distance


def polygon_intersect(data_frommap_merged, line_array, calcdist_fromlatlon=None):
    #data_frommap_merged: xugrid dataset (contains ds and grid)
    #TODO: remove hardcoding
    """
    #TODO: maybe move to meshkernel functionality?
    Cross section functionality is implemented in MeshKernel (C++) but still needs to be exposed in MeshKernelPy (can be done in dec2022). Here is the function with documentation: 
    https://github.com/Deltares/MeshKernel/blob/067f1493e7f972ba0cdb2a1f4deb48d1c74695d5/include/MeshKernelApi/MeshKernel.hpp#L356
    """
    
    import numpy as np
    #from matplotlib.path import Path
    import shapely #separate import, since sometimes this works, while import shapely.geometry fails
    from shapely.geometry import LineString, Polygon, MultiLineString, Point
    from dfm_tools.get_nc import calc_dist_pythagoras, calc_dist_haversine

    if calcdist_fromlatlon is None:
        #auto determine if cartesian/sperical
        if hasattr(data_frommap_merged.ugrid.obj,'projected_coordinate_system'):
            calcdist_fromlatlon = False
        elif hasattr(data_frommap_merged.ugrid.obj,'wgs84'):
            calcdist_fromlatlon = True
        else:
            raise Exception('To auto determine calcdist_fromlatlon, a variable "projected_coordinate_system" or "wgs84" is required, please provide calcdist_fromlatlon=True/False yourself.')

    dtstart_all = dt.datetime.now()

    #defining celinlinebox
    line_section = LineString(line_array)
    
    ugrid_all_verts = get_ugrid_verts(data_frommap_merged)
    verts_xmax = np.nanmax(ugrid_all_verts[:,:,0].data,axis=1)
    verts_xmin = np.nanmin(ugrid_all_verts[:,:,0].data,axis=1)
    verts_ymax = np.nanmax(ugrid_all_verts[:,:,1].data,axis=1)
    verts_ymin = np.nanmin(ugrid_all_verts[:,:,1].data,axis=1)
    
    #TODO: replace this with xr.sel() once it works for xugrid (getting verts_inlinebox_nos is than still an issue)
    cellinlinebox_all_bool = (((np.min(line_array[:,0]) <= verts_xmax) &
                               (np.max(line_array[:,0]) >= verts_xmin)) &
                              ((np.min(line_array[:,1]) <= verts_ymax) & 
                               (np.max(line_array[:,1]) >= verts_ymin))
                              )
    #cellinlinebox_all_bool[:] = 1 #to force all cells to be taken into account
    
    intersect_coords = np.empty((0,4))
    intersect_gridnos = np.empty((0),dtype=int) #has to be numbers, since a boolean is differently ordered
    verts_inlinebox = ugrid_all_verts[cellinlinebox_all_bool,:,:]
    verts_inlinebox_nos = np.where(cellinlinebox_all_bool)[0]

    
    print(f'>> finding crossing flow links (can take a while, processing {cellinlinebox_all_bool.sum()} of {len(cellinlinebox_all_bool)} cells): ',end='')
    dtstart = dt.datetime.now()
    for iP, pol_data in enumerate(verts_inlinebox):
        pol_shp = Polygon(pol_data[~np.isnan(pol_data).all(axis=1)])
        intersect_result = pol_shp.intersection(line_section)
        if isinstance(intersect_result,shapely.geometry.multilinestring.MultiLineString): #in the rare case that a cell (pol_shp) is crossed by multiple parts of the line
            intersect_result_multi = intersect_result
        elif isinstance(intersect_result,shapely.geometry.linestring.LineString): #if one linepart trough cell (ex/including node), make multilinestring anyway
            if intersect_result.coords == []: #when the line does not cross this cell, intersect_results.coords is an empty linestring and this cell can be skipped (continue makes forloop continue with next in line without finishing the rest of the steps for this instance)
                continue
            elif len(intersect_result.coords.xy[0]) == 0: #for newer cartopy versions, when line does not cross this cell, intersect_result.coords.xy is (array('d'), array('d')), and both arrays in tuple have len 0.
                continue
            intersect_result_multi = MultiLineString([intersect_result])
        for iLL, intesect_result_one in enumerate(intersect_result_multi.geoms): #loop over multilinestrings, will mostly only contain one linestring. Will be two if the line crosses a cell more than once.
            intersection_line = intesect_result_one.coords
            intline_xyshape = np.array(intersection_line.xy).shape
            #print('len(intersection_line.xy): %s'%([intline_xyshape]))
            for numlinepart_incell in range(1,intline_xyshape[1]): #is mostly 1, but more if there is a linebreakpoint in this cell (then there are two or more lineparts)
                intersect_gridnos = np.append(intersect_gridnos,verts_inlinebox_nos[iP])
                #intersect_coords = np.concatenate([intersect_coords,np.array(intersection_line.xy)[np.newaxis,:,numlinepart_incell-1:numlinepart_incell+1]],axis=0)
                intersect_coords = np.concatenate([intersect_coords,np.array(intersection_line.xy).T[numlinepart_incell-1:numlinepart_incell+1].flatten()[np.newaxis]])
    
    if intersect_coords.shape[0] != len(intersect_gridnos):
        raise Exception('something went wrong, intersect_coords.shape[0] and len(intersect_gridnos) are not equal')
    
    import pandas as pd
    intersect_pd = pd.DataFrame(intersect_coords,index=intersect_gridnos,columns=['x1','y1','x2','y2'])
    intersect_pd.index.name = 'gridnumber'
    print(f'{(dt.datetime.now()-dtstart).total_seconds():.2f} sec')
            
    
    #calculating distance for all crossed cells, from first point of line
    nlinecoords = line_array.shape[0]
    nlinedims = len(line_array.shape)
    ncrosscellparts = len(intersect_pd)
    if nlinecoords<2 or nlinedims != 2:
        raise Exception('ERROR: line_array should at least contain two xy points [[x,y],[x,y]]')
    
    #calculate distance between celledge-linepart crossing (is zero when line iL crosses cell)
    distperline_tostart = np.zeros((ncrosscellparts,nlinecoords-1))
    distperline_tostop = np.zeros((ncrosscellparts,nlinecoords-1))
    linepart_length = np.zeros((nlinecoords))
    for iL in range(nlinecoords-1):
        #calculate length of lineparts
        line_section_part = LineString(line_array[iL:iL+2,:])
        if calcdist_fromlatlon:
            linepart_length[iL+1] = calc_dist_haversine(line_array[iL,0],line_array[iL+1,0],line_array[iL,1],line_array[iL+1,1])
        else:
            linepart_length[iL+1] = line_section_part.length
    
        #get distance between all lineparts and point (later used to calculate distance from beginpoint of closest linepart)
        for iP in range(ncrosscellparts):
            distperline_tostart[iP,iL] = line_section_part.distance(Point(intersect_coords[:,0][iP],intersect_coords[:,1][iP]))
            distperline_tostop[iP,iL] = line_section_part.distance(Point(intersect_coords[:,2][iP],intersect_coords[:,3][iP]))
    linepart_lengthcum = np.cumsum(linepart_length)
    cross_points_closestlineid = np.argmin(np.maximum(distperline_tostart,distperline_tostop),axis=1)
    intersect_pd['closestlineid'] = cross_points_closestlineid
    
    
    if not calcdist_fromlatlon:
        crs_dist_starts = calc_dist_pythagoras(line_array[cross_points_closestlineid,0], intersect_coords[:,0], line_array[cross_points_closestlineid,1], intersect_coords[:,1]) + linepart_lengthcum[cross_points_closestlineid]
        crs_dist_stops = calc_dist_pythagoras(line_array[cross_points_closestlineid,0], intersect_coords[:,2], line_array[cross_points_closestlineid,1], intersect_coords[:,3]) + linepart_lengthcum[cross_points_closestlineid]
    else:
        crs_dist_starts = calc_dist_haversine(line_array[cross_points_closestlineid,0], intersect_coords[:,0], line_array[cross_points_closestlineid,1], intersect_coords[:,1]) + linepart_lengthcum[cross_points_closestlineid]
        crs_dist_stops = calc_dist_haversine(line_array[cross_points_closestlineid,0], intersect_coords[:,2], line_array[cross_points_closestlineid,1], intersect_coords[:,3]) + linepart_lengthcum[cross_points_closestlineid]
    intersect_pd['crs_dist_starts'] = crs_dist_starts
    intersect_pd['crs_dist_stops'] = crs_dist_stops
    intersect_pd['linepartlen'] = crs_dist_stops-crs_dist_starts
    intersect_pd = intersect_pd.sort_values('crs_dist_starts')
    
    print(f'>> polygon_intersect() total, found intersection trough {len(intersect_gridnos)} of {len(cellinlinebox_all_bool)} cells: {(dt.datetime.now()-dtstart_all).total_seconds():.2f} sec')
    return intersect_pd


def get_xzcoords_onintersection(data_frommap_merged, intersect_pd, timestep=None):
    #TODO: remove hardcoding of variable names
    if timestep is None: #TODO: maybe make time dependent grid?
        raise Exception('ERROR: argument timestep not provided, this is necessary to retrieve correct waterlevel or fullgrid output')
    
    dimn_layer, dimn_interfaces = get_vertical_dimensions(data_frommap_merged)
    if dimn_layer is not None:
        nlay = data_frommap_merged.dims[dimn_layer]
    else: #no layers, 2D model
        nlay = 1

    xu_facedim = data_frommap_merged.grid.face_dimension
    xu_edgedim = data_frommap_merged.grid.edge_dimension
    xu_nodedim = data_frommap_merged.grid.node_dimension
        
    #potentially construct fullgrid info (zcc/zw) #TODO: this ifloop is copied from get_mapdata_atdepth(), prevent this duplicate code
    if dimn_layer not in data_frommap_merged.dims: #2D model
        print('depth dimension not found, probably 2D model')
        pass
    elif 'mesh2d_flowelem_zw' in data_frommap_merged.variables: #fullgrid info already available, so continuing
        print('zw/zcc (fullgrid) values already present in Dataset')
        pass
    elif 'mesh2d_layer_sigma' in data_frommap_merged.variables: #reconstruct_zw_zcc_fromsigma and treat as zsigma/fullgrid mapfile from here
        print('sigma-layer model, computing zw/zcc (fullgrid) values and treat as fullgrid model from here')
        data_frommap_merged = reconstruct_zw_zcc_fromsigma(data_frommap_merged)
    elif 'mesh2d_layer_z' in data_frommap_merged.variables:        
        print('z-layer model, computing zw/zcc (fullgrid) values and treat as fullgrid model from here')
        data_frommap_merged = reconstruct_zw_zcc_fromz(data_frommap_merged)
    else:
        raise Exception('layers present, but unknown layertype, expected one of variables: mesh2d_flowelem_zw, mesh2d_layer_sigma, mesh2d_layer_z')
    
    intersect_gridnos = intersect_pd.index
    data_frommap_merged_sel = data_frommap_merged.ugrid.obj.drop_dims([xu_edgedim,xu_nodedim]).isel(time=timestep).isel({xu_facedim:intersect_gridnos}) #TODO: not possible to do isel with non-sorted gridnos on ugridDataset, but it is on xrDataset. Dropping edge/nodedims first for neatness, so there is not mismatch in face/node/edge after using .isel(faces)
    if dimn_layer not in data_frommap_merged_sel.dims:
        data_frommap_wl3_sel = data_frommap_merged_sel['mesh2d_s1'].to_numpy()
        data_frommap_bl_sel = data_frommap_merged_sel['mesh2d_flowelem_bl'].to_numpy()
        zvals_interface = np.linspace(data_frommap_bl_sel,data_frommap_wl3_sel,nlay+1)
    elif 'mesh2d_flowelem_zw' in data_frommap_merged_sel.variables:
        zvals_interface_filled = data_frommap_merged_sel['mesh2d_flowelem_zw'].bfill(dim=dimn_interfaces) #fill nan values (below bed) with equal values
        zvals_interface = zvals_interface_filled.to_numpy().T # transpose to make in line with 2D sigma dataset
    
    #convert to output for plot_netmapdata
    crs_dist_starts_matrix = np.repeat(intersect_pd['crs_dist_starts'].values[np.newaxis],nlay,axis=0)
    crs_dist_stops_matrix = np.repeat(intersect_pd['crs_dist_stops'].values[np.newaxis],nlay,axis=0)
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
    crs_plotdata_clean = data_frommap_merged_sel
    if dimn_layer in data_frommap_merged_sel.dims:
        facedim_tempname = 'facedim_tempname' #temporary new name to avoid duplicate from-to dimension name in .stack()
        crs_plotdata_clean = crs_plotdata_clean.rename({xu_facedim:facedim_tempname})
        crs_plotdata_clean = crs_plotdata_clean.stack({xr_crs_grid.face_dimension:[dimn_layer,facedim_tempname]})
        #reset_index converts face-dimension from multiindex to flat
        crs_plotdata_clean = crs_plotdata_clean.reset_index([xr_crs_grid.face_dimension])
    
    #combine into xugrid
    xr_crs_ugrid = xu.UgridDataset(crs_plotdata_clean, grids=[xr_crs_grid])
    return xr_crs_ugrid


def polyline_mapslice(data_frommap_merged, line_array, timestep, calcdist_fromlatlon=None): #TODO: merge this into one function
    #intersect function, find crossed cell numbers (gridnos) and coordinates of intersection (2 per crossed cell)
    intersect_pd = polygon_intersect(data_frommap_merged, line_array, calcdist_fromlatlon=calcdist_fromlatlon)
    if len(intersect_pd) == 0:
        raise Exception('line_array does not cross mapdata') #TODO: move exception elsewhere?
    #derive vertices from cross section (distance from first point)
    xr_crs_ugrid = get_xzcoords_onintersection(data_frommap_merged, intersect_pd=intersect_pd, timestep=timestep)
    return xr_crs_ugrid


def reconstruct_zw_zcc_fromsigma(data_xr_map):
    """
    reconstruct full grid output (time/face-varying z-values) for sigma model, necessary for slicing sigmamodel on depth value
    """
    data_frommap_wl_sel = data_xr_map['mesh2d_s1']
    data_frommap_bl_sel = data_xr_map['mesh2d_flowelem_bl']
    
    zvals_cen_percentage = data_xr_map['mesh2d_layer_sigma']
    data_xr_map['mesh2d_flowelem_zcc'] = data_frommap_wl_sel+(data_frommap_wl_sel-data_frommap_bl_sel)*zvals_cen_percentage
    
    zvals_interface_percentage = data_xr_map['mesh2d_interface_sigma']
    data_xr_map['mesh2d_flowelem_zw'] = data_frommap_wl_sel+(data_frommap_wl_sel-data_frommap_bl_sel)*zvals_interface_percentage
    
    data_xr_map = data_xr_map.set_coords(['mesh2d_flowelem_zw','mesh2d_flowelem_zcc'])
    return data_xr_map


def reconstruct_zw_zcc_fromz(data_xr_map):
    """
    reconstruct full grid output (time/face-varying z-values) for zvalue model. Necessary when extracting values with zdepth w.r.t. waterlevel/bedlevel
    #TODO: gives spotty result for 0/0.1m w.r.t. bedlevel for Grevelingen zmodel
    #TODO: remove hardcoding of varnames
    """
    
    dimn_layer, dimn_interfaces = get_vertical_dimensions(data_xr_map)
    
    data_frommap_wl_sel = data_xr_map['mesh2d_s1']
    data_frommap_z0_sel = data_frommap_wl_sel*0
    data_frommap_bl_sel = data_xr_map['mesh2d_flowelem_bl']
    
    zvals_cen_zval = data_xr_map['mesh2d_layer_z'] #no clipping for zcenter values, since otherwise interp will fail
    data_xr_map['mesh2d_flowelem_zcc'] = (data_frommap_z0_sel+zvals_cen_zval)

    zvals_interface_zval = data_xr_map['mesh2d_interface_z'] #clipping for zinterface values, to make sure layer interfaces are also at water/bed level
    data_xr_map['mesh2d_flowelem_zw'] = (data_frommap_z0_sel+zvals_interface_zval).clip(min=data_frommap_bl_sel, max=data_frommap_wl_sel)
    bool_notoplayer_int = zvals_interface_zval<zvals_interface_zval.isel({dimn_interfaces:-1})
    bool_int_abovewl = zvals_interface_zval>data_frommap_wl_sel
    data_xr_map['mesh2d_flowelem_zw'] = data_xr_map['mesh2d_flowelem_zw'].where(bool_notoplayer_int | bool_int_abovewl, other=data_frommap_wl_sel) #zvalues of top layer_interfaces that are lower than wl are replaced by wl
    
    data_xr_map = data_xr_map.set_coords(['mesh2d_flowelem_zw','mesh2d_flowelem_zcc'])
    return data_xr_map


def get_Dataset_atdepths(data_xr, depths, reference='z0', zlayer_z0_selnearest=False):
    """
    data_xr_map:
        has to be Dataset (not a DataArray), otherwise mesh2d_flowelem_zw etc are not available (interface z values)
        in case of zsigma/sigma layers (or fullgrid), it is advisable to .sel()/.isel() the time dimension first, because that is less computationally heavy
    depths:
        int/float or list/array of int/float. Depths w.r.t. reference level. If reference=='waterlevel', depth>0 returns only nans. If reference=='bedlevel', depth<0 returns only nans. Depths are sorted and only uniques are kept.
    reference:
        compute depth w.r.t. z0/waterlevel/bed
        default: reference='z0'
    zlayer_z0_interp:
        Use xr.interp() to interpolate zlayer model to z-value. Only possible for reference='z' (not 'waterlevel' or 'bedlevel'). Only used if "mesh2d_layer_z" is present (zlayer model)
        This is faster but results in values interpolated between zcc (z cell centers), so it is different than slicing.
    
    #TODO: zmodel gets depth in figure title, because of .set_index() in open_partitioned_dataset(). Sigmamodel gets percentage/fraction in title. >> set_index was removed there, so check again.
    #TODO: check if attributes should be passed/altered
    #TODO: build in check whether layers/interfaces are coherent (len(int)=len(lay)+1), since one could also supply a uds.isel(layer=range(10)) where there are still 51 interfaces. (if not, raise error to remove layer selection)
    #TODO: also waterlevelvar in 3D model gets depth_fromref coordinate, would be nice to avoid.
    #TODO: clean up unneccesary variables (like pre-interp depth values and interface dim)
    """
    
    depth_varname = 'depth_fromref'
    
    dimn_layer, dimn_interfaces = get_vertical_dimensions(data_xr) #TODO: make more generic to also work with hisfiles?
    
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
        warnings.warn(UserWarning('get_Dataset_atdepths() is not tested for hisfiles yet, please check your results.'))
    else:
        print('WARNING: depth dimension not found, probably 2D model, returning input Dataset')
        return data_xr #early return
    
    if not isinstance(data_xr,(xr.Dataset,xu.UgridDataset)):
        raise Exception(f'data_xr_map should be of type xr.Dataset, but is {type(data_xr)}')
    
    #create depth xr.DataArray
    if isinstance(depths,(float,int)):
        depths = depths #float/int
        depth_dims = ()
    else:
        depths = np.unique(depths) #array of unique+sorted floats/ints
        depth_dims = (depth_varname)
    depths_xr = xr.DataArray(depths,dims=depth_dims,attrs={'units':'m',
                                                           'reference':f'model_{reference}',
                                                           'positive':'up'}) #TODO: make more in line with CMEMS etc
    
    #extract waterlevels and bedlevels from file, to correct layers with
    data_wl = data_xr[varname_wl]
    data_bl = data_xr[varname_bl]
    
    #simplified/faster method for zlayer icm z0 reference (mapfiles only)
    if 'mesh2d_layer_z' in data_xr.variables and zlayer_z0_selnearest and reference=='z0': # selects nearest z-center values (instead of slicing), should be faster #TODO: check if this is faster than fullgrid
        print('z-layer model, zlayer_z0_selnearest=True and reference=="z0" so using xr.sel(method="nearest")]')
        data_xr = data_xr.set_index({dimn_layer:'mesh2d_layer_z'})#.rename({'nmesh2d_layer':depth_varname}) #set depth as index on layers, to be able to interp to depths instead of layernumbers
        data_xr[depth_varname] = depths_xr
        data_xr_atdepths = data_xr.sel({dimn_layer:depths_xr},method='nearest')
        data_xr_atdepths = data_xr_atdepths.where((depths_xr>=data_bl) & (depths_xr<=data_wl)) #filter above wl and below bl values
        return data_xr_atdepths #early return
    
    #potentially construct fullgrid info (zcc/zw) #TODO: maybe move to separate function, like open_partitioned_dataset() (although bl/wl are needed anyway)
    if varname_zint in data_xr.variables: #fullgrid info already available, so continuing
        print(f'zw/zcc (fullgrid) values already present in Dataset in variable {varname_zint}')
        pass
    elif 'mesh2d_layer_sigma' in data_xr.variables: #reconstruct_zw_zcc_fromsigma and treat as zsigma/fullgrid mapfile from here
        print('sigma-layer model, computing zw/zcc (fullgrid) values and treat as fullgrid model from here')
        data_xr = reconstruct_zw_zcc_fromsigma(data_xr)
    elif 'mesh2d_layer_z' in data_xr.variables:
        print('z-layer model, computing zw/zcc (fullgrid) values and treat as fullgrid model from here')
        data_xr = reconstruct_zw_zcc_fromz(data_xr)
    else:
        raise Exception('layers present, but unknown layertype/var')
    
    #correct reference level
    if reference=='z0':
        zw_reference = data_xr[varname_zint]
    elif reference=='waterlevel':
        zw_reference = data_xr[varname_zint] - data_wl
    elif reference=='bedlevel':
        zw_reference = data_xr[varname_zint] - data_bl
    else:
        raise Exception(f'unknown reference "{reference}" (possible are z0, waterlevel and bedlevel')
    
    print('>> subsetting data on fixed depth in fullgrid z-data: ',end='')
    dtstart = dt.datetime.now()
        
    if 'time' in data_xr.dims: #TODO: suppress this warning for hisfiles since it does not make sense
        warnings.warn(UserWarning('get_mapdata_onfixedepth() can be very slow when supplying dataset with time dimension, you could supply ds.isel(time=timestep) instead'))
        
    #get layerno via z-interface value (zw), check which celltop-interfaces are above/on depth and which which cellbottom-interfaces are below/on depth
    bool_topinterface_abovedepth = zw_reference.isel({dimname_layw:slice(1,None)}) >= depths_xr
    bool_botinterface_belowdepth = zw_reference.isel({dimname_layw:slice(None,-1)}) <= depths_xr
    bool_topbotinterface_arounddepth = bool_topinterface_abovedepth & bool_botinterface_belowdepth #this bool also automatically excludes all values below bed and above wl
    bool_topbotinterface_arounddepth = bool_topbotinterface_arounddepth.rename({dimname_layw:dimname_layc}) #correct dimname for interfaces to centers
    #bool_topbotinterface_arounddepth.sum(dim='nmesh2d_layer').max(dim='time').load() #TODO: max is layno, should be 1 everywhere?
    data_xr_atdepths = data_xr.where(bool_topbotinterface_arounddepth).max(dim=dimname_layc,keep_attrs=True) #set all layers but one to nan, followed by an arbitrary reduce (max in this case)
    #TODO: suppress/solve warning for DCSM (does not happen when supplying data_xr_map[['mesh2d_sa1','mesh2d_s1','mesh2d_flowelem_bl','mesh2d_flowelem_zw']]): "C:\Users\veenstra\Anaconda3\envs\dfm_tools_env\lib\site-packages\dask\array\core.py:4806: PerformanceWarning: Increasing number of chunks by factor of 20"
    #TODO: suppress warning (upon plotting/load/etc): "C:\Users\veenstra\Anaconda3\envs\dfm_tools_env\lib\site-packages\dask\array\reductions.py:640: RuntimeWarning: All-NaN slice encountered"
    
    #add depth as coordinate var
    data_xr_atdepths[depth_varname] = depths_xr
    data_xr_atdepths = data_xr_atdepths.set_coords([depth_varname])
    
    print(f'{(dt.datetime.now()-dtstart).total_seconds():.2f} sec')
    
    return data_xr_atdepths


def get_netdata(file_nc, multipart=None):

    warnings.warn(DeprecationWarning('dfm_tools.get_nc.get_netdata() is deprecated, since there is an xarray alternative for multidomain FM files (xugrid). Open it like this and use xarray sel/isel (example in postprocessing notebook):\n    data_xr_mapmerged = dfmt.open_partitioned_dataset(file_nc_map)'))
    #file_ncs = get_ncfilelist(file_nc, multipart)
    if isinstance(file_nc,list): #for opendap, has to support lists
        file_ncs = file_nc
    else:
        file_ncs = glob.glob(file_nc)
    
    if multipart is not None:
        warnings.warn(UserWarning('argument multipart is deprecated and will be ignored'))
        
    #get all data
    num_nodes = [0]
    verts_shape2_all = []
    print('processing %d partitions (first getting max number of facenodes)'%(len(file_ncs)))
    for iF, file_nc_sel in enumerate(file_ncs):
        data_nc = Dataset(file_nc_sel)
        varn_mesh2d_face_nodes = get_varname_fromnc(data_nc,'mesh2d_face_nodes',vardim='var')
        if varn_mesh2d_face_nodes is not None: # node_z variable is present
            mesh2d_face_nodes = data_nc.variables[varn_mesh2d_face_nodes]
        else:
            raise Exception('ERROR: provided file does not contain a variable mesh2d_face_nodes or similar:\n%s\nPlease do one of the following:\n- plot grid from *_map.nc file\n- import and export the grid with RGFGRID\n- import and save the grid "with cellfinfo" from interacter'%(file_nc))
        verts_shape2_all.append(mesh2d_face_nodes.shape[1])
        data_nc.close()
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
    warnings.warn(DeprecationWarning('dfm_tools.get_nc.plot_netmapdata() will be deprecated, since there is an xarray alternative. Like this (example in postprocessing notebook):\n   xr_crs_ugrid = dfmt.polyline_mapslice(data_frommap_merged, line_array, timestep=timestep)\n    fig, ax = plt.subplots()\n    xr_crs_ugrid["mesh2d_sa1"].ugrid.plot.line()'))

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

    import cartopy
    import cartopy.crs as ccrs
    import cartopy.io.img_tiles as cimgt
    import cartopy.feature as cfeature
    import cartopy.mpl.ticker as cticker

    dummy = ccrs.epsg(28992) #to make cartopy realize it has a cartopy._epsg._EPSGProjection class (maybe gets fixed with cartopy updates, see unittest test_cartopy_epsg)
    if ax is None: #provide axis projection on initialisation, cannot be edited later on
        if projection is None:
            projection=ccrs.PlateCarree() #projection of cimgt.GoogleTiles, useful default
        elif isinstance(projection, (cartopy._epsg._EPSGProjection, cartopy.crs.CRS)): #checks if argument is an EPSG projection or CRS projection (like PlateCarree, Mercator etc). Note: it was cartopy._crs.CRS before instead of cartopy.crs.CRS
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


def plot_ztdata(data_xr_sel, varname, ax=None, only_contour=False, get_ds=False, **kwargs):
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
    **kwargs : TYPE
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
        raise Exception(f'ERROR: unexpected number of dimensions in requested squeezed variable ({data_xr_sel[varname].shape}), first use data_xr.isel(stations=int) to select a single station') #TODO: can also have a different cause, improve message/testing?
    
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
        #pc = data_xr_sel[varname].plot.pcolormesh(ax=ax, x='time', y='zcoordinate_w', **kwargs) #is not possible to put center values on interfaces, som more difficult approach needed
        pc = ax.pcolormesh(time_mesh_cor, data_fromhis_zcor, data_xr_sel[varname], **kwargs)
   
    return pc

