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
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.collections

from dfm_tools.get_nc_helpers import get_ncfilelist, get_ncvarproperties, get_varnamefrom_keyslongstandardname, get_timesfromnc, get_timeid_fromdatetime, get_hisstationlist, get_stationid_fromstationlist, ghostcell_filter, get_varname_fromnc
from dfm_tools.ugrid import UGrid


def get_ncmodeldata(file_nc, varname=None, timestep=None, layer=None, station=None, multipart=None, silent=False, return_xarray=False):
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
    
    warnings.warn(DeprecationWarning('dfm_tools.get_nc.get_ncmodeldata() will be deprecated, since there is an xarray alternative for multidomain FM files (xugrid). Check the example scripts for how to use it'))
    
    #get variable info (also checks if varname exists in keys, standard name, long name)
    if isinstance(file_nc,list): #for opendap, has to support lists
        file_nc_one = file_nc[0]
    else:
        file_nc_one = file_nc        
    data_nc = Dataset(file_nc_one)
    data_xr = xr.open_dataset(file_nc_one)
    varname = get_varnamefrom_keyslongstandardname(file_nc_one, varname) #get varname from varkeys/standardname/longname if exists
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
    if not bool_varpartitioned:
        multipart = False

    #get list of partitioned files
    if isinstance(file_nc,list):
        file_ncs = file_nc
    else:
        file_ncs = get_ncfilelist(file_nc_one, multipart)

    for iF, file_nc_sel in enumerate(file_ncs):
        if (len(file_ncs) > 1) and not silent:
            print('processing mapdata from domain %04d of %04d'%(iF, len(file_ncs)-1))

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
    
    if return_xarray:
        data_xr = xr.Dataset()
        var_xr = xr.DataArray(values_all,dims=nc_varobject.dimensions,attrs=nc_varobject.__dict__,name=varname)
        data_xr[varname] = var_xr
        for vardim in nc_varobject.dimensions:
            if vardim=='time':
                data_xr['time'] = xr.DataArray(data_nc_datetimes_pd, dims=('time'), attrs=data_nc_timevar.__dict__, name='time')
            if vardim=='mesh2d_nLayers':
                data_xr['mesh2d_nLayers'] = xr.DataArray(layer_ids, dims=('mesh2d_nLayers'), name='mesh2d_nLayers')
        data_nc.close()
        return data_xr[varname]
    else:
        data_nc.close()
        return values_all


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


def polygon_intersect(data_frommap_merged, line_array, optimize_dist=False, calcdist_fromlatlon=False): #TODO: copy of ugrid function, remove ugrid when xugrid works for all stuff
    #data_frommap_merged: xugrid dataset (contains ds and grid)
    #TODO: remove hardcoding
    import numpy as np
    from matplotlib.path import Path
    import shapely #separate import, since sometimes this works, while import shapely.geometry fails
    from shapely.geometry import LineString, Polygon, MultiLineString, Point
    from dfm_tools.get_nc import calc_dist_pythagoras, calc_dist_haversine

    print('defining celinlinebox')
    
    line_section = LineString(line_array)
    
    face_nos = data_frommap_merged.ugrid.grid.to_dataset().mesh2d_face_nodes.load()
    bool_nonemptyfacenode = face_nos!=-1
    facenos_nonan_min = face_nos.where(face_nos!=-1).min() #replace nans and get minval
    if facenos_nonan_min==1: #for some reason, curvedbend is 1-based indexed, grevelingen is not
        face_nos = face_nos-1
    if face_nos.dtype!='int': #for some reason, curvedbend idx is float instead of int
        face_nos = face_nos.astype(int)
           
    face_nnodecoords_x = data_frommap_merged.ugrid.grid.to_dataset().mesh2d_node_x.isel(mesh2d_nNodes=face_nos).where(bool_nonemptyfacenode)
    face_nnodecoords_y = data_frommap_merged.ugrid.grid.to_dataset().mesh2d_node_y.isel(mesh2d_nNodes=face_nos).where(bool_nonemptyfacenode)
    ugrid_all_verts = np.c_[face_nnodecoords_x.to_numpy()[...,np.newaxis],face_nnodecoords_y.to_numpy()[...,np.newaxis]]
    
    # verts_xmax = np.nanmax(ugrid_all.verts[:,:,0].data,axis=1)
    # verts_xmin = np.nanmin(ugrid_all.verts[:,:,0].data,axis=1)
    # verts_ymax = np.nanmax(ugrid_all.verts[:,:,1].data,axis=1)
    # verts_ymin = np.nanmin(ugrid_all.verts[:,:,1].data,axis=1)
    verts_xmax = np.nanmax(face_nnodecoords_x.to_numpy(),axis=1)
    verts_xmin = np.nanmin(face_nnodecoords_x.to_numpy(),axis=1)
    verts_ymax = np.nanmax(face_nnodecoords_y.to_numpy(),axis=1)
    verts_ymin = np.nanmin(face_nnodecoords_y.to_numpy(),axis=1)
    
    if not optimize_dist: #TODO: replace this with xr.sel() once it works for xugrid
        cellinlinebox_all_bool = (((np.min(line_array[:,0]) <= verts_xmax) &
                                   (np.max(line_array[:,0]) >= verts_xmin)) &
                                  ((np.min(line_array[:,1]) <= verts_ymax) & 
                                   (np.max(line_array[:,1]) >= verts_ymin))
                                  )
    elif type(optimize_dist) in [int,float]: #not properly tested and documented
        #calculate angles wrt x axis
        angles_wrtx = []
        nlinecoords = line_array.shape[0]
        for iL in range(nlinecoords-1):
            dx = line_array[iL+1,0] - line_array[iL,0]
            dy = line_array[iL+1,1] - line_array[iL,1]
            angles_wrtx.append(np.rad2deg(np.arctan2(dy,dx)))
        angles_toprev = np.concatenate([[90],np.diff(angles_wrtx),[90]])
        angles_wrtx_ext = np.concatenate([[angles_wrtx[0]-90],np.array(angles_wrtx),[angles_wrtx[-1]+90]])
        angtot_wrtx = angles_wrtx_ext[:-1] + 0.5*(180+angles_toprev)
        #distance over xy-axis from original points
        dxynewpoints = optimize_dist * np.array([np.cos(np.deg2rad(angtot_wrtx)),np.sin(np.deg2rad(angtot_wrtx))]).T
        newpoints1 = line_array+dxynewpoints
        newpoints2 = line_array-dxynewpoints
        pol_inpol = np.concatenate([newpoints1, np.flip(newpoints2,axis=0)])
        pol_inpol_path = Path(pol_inpol)
        bool_all = []
        for iC in range(ugrid_all_verts.shape[1]):
            data_arr_ic = ugrid_all_verts[:,iC,:]
            test = pol_inpol_path.contains_points(data_arr_ic)
            bool_all.append(test)
        test_all = np.array(bool_all)
        cellinlinebox_all_bool = (test_all==True).any(axis=0)
    else:
        raise Exception('ERROR: invalid type for optimize_dist argument')
    
    intersect_coords = np.empty((0,4))
    intersect_gridnos = np.empty((0),dtype=int) #has to be numbers, since a boolean is differently ordered
    verts_inlinebox = ugrid_all_verts[cellinlinebox_all_bool,:,:]
    verts_inlinebox_nos = np.where(cellinlinebox_all_bool)[0]
    print('finding crossing flow links (can take a while if linebox over xy covers a lot of cells, %i of %i cells are being processed)'%(cellinlinebox_all_bool.sum(),len(cellinlinebox_all_bool)))
    
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
            
    print('calculating distance for all crossed cells, from first point of line (should not take long, but if it does, optimisation is needed)')
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
    print('finished calculating distance for all crossed cells, from first point of line')
    
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
    
    #dimensions (gridnos, xy, firstsecond)
    print('done finding crossing flow links: %i of %i'%(len(intersect_gridnos),len(cellinlinebox_all_bool)))
    return intersect_pd


def get_xzcoords_onintersection(data_frommap_merged, intersect_pd, timestep=None, varname=None):
    #TODO: no hardcoding of variable names?
    #check if all necessary arguments are provided
    if timestep is None:
        raise Exception('ERROR: argument timestep not provided, this is necessary to retrieve correct waterlevel or fullgrid output')
    
    varkeys_list = list(data_frommap_merged.variables.keys())
    dimn_layer = 'nmesh2d_layer'#dfmt.get_varname_fromnc(data_frommap_merged,'nmesh2d_layer',vardim='dim') #TODO: hardcoding ok because of renaming? (no, since gridname is not fixed)
    if dimn_layer in data_frommap_merged.dims:
        nlay = data_frommap_merged.dims[dimn_layer]
    else: #no layers, 2D model
        nlay = 1
    
    intersect_gridnos = intersect_pd.index
    data_frommap_merged_sel = data_frommap_merged.isel(time=timestep,mesh2d_nFaces=intersect_gridnos)
    if 'mesh2d_flowelem_zw' in varkeys_list:
        print('layertype: fullgrid output')
        zvals_interface = data_frommap_merged_sel['mesh2d_flowelem_zw'].to_numpy().T # transpose to make in line with 2D sigma dataset
    else: #no full grid output, so reconstruct
        data_frommap_wl3_sel = data_frommap_merged_sel['mesh2d_s1'].to_numpy()
        data_frommap_bl_sel = data_frommap_merged_sel['mesh2d_flowelem_bl'].to_numpy()
        if 'mesh2d_layer_z' in varkeys_list or 'LayCoord_cc' in varkeys_list:
            print('layertype: zlayer')
            warnings.warn('WARNING: your model seems to contain only z-layers. if the modeloutput is generated with an older version of dflowfm, the coordinates can be incorrect. if your model contains z-sigma-layers, use the fulloutput option in the mdu and rerun (happens automatically in newer dflowfm versions).')
            zvals_interface_vec = data_frommap_merged_sel['mesh2d_interface_z'].to_numpy()[:,np.newaxis]
            zvals_interface = np.repeat(zvals_interface_vec,len(data_frommap_wl3_sel),axis=1)
            # zvalues lower than bedlevel should be overwritten with bedlevel
            for iL in range(nlay):
                zvalbot_belowbl_bool = zvals_interface[iL,:]<data_frommap_bl_sel
                zvals_interface[iL,zvalbot_belowbl_bool] = data_frommap_bl_sel[zvalbot_belowbl_bool]
            #top z-layer is extended to water level, if wl is higher than zval_lay_top
            zvals_interface[-1,:] = np.maximum(zvals_interface[-1,:],data_frommap_wl3_sel)
        elif 'mesh2d_layer_sigma' in varkeys_list:
            print('layertype: sigmalayer')
            zvals_interface_percentage = data_frommap_merged_sel['mesh2d_interface_sigma'].to_numpy()[:,np.newaxis]
            zvals_interface = data_frommap_wl3_sel+(data_frommap_wl3_sel-data_frommap_bl_sel)[np.newaxis]*zvals_interface_percentage
        else: # 2D model
            print('layertype: 2D model')
            if nlay!=1:
                raise Exception('recheck this')
            #zvals_cen = np.linspace(data_frommap_bl_sel,data_frommap_wl3_sel,nlay)
            zvals_interface = np.linspace(data_frommap_bl_sel,data_frommap_wl3_sel,nlay+1)
    
    #convert to output for plot_netmapdata
    crs_dist_starts_matrix = np.repeat(intersect_pd['crs_dist_starts'].values[np.newaxis],nlay,axis=0)
    crs_dist_stops_matrix = np.repeat(intersect_pd['crs_dist_stops'].values[np.newaxis],nlay,axis=0)
    crs_verts_x_all = np.array([[crs_dist_starts_matrix.ravel(),crs_dist_stops_matrix.ravel(),crs_dist_stops_matrix.ravel(),crs_dist_starts_matrix.ravel()]]).T
    crs_verts_z_all = np.ma.array([zvals_interface[1:,:].ravel(),zvals_interface[1:,:].ravel(),zvals_interface[:-1,:].ravel(),zvals_interface[:-1,:].ravel()]).T[:,:,np.newaxis]
    crs_verts = np.ma.concatenate([crs_verts_x_all, crs_verts_z_all], axis=2)
    
    if varname is not None: #retrieve data for varname and return
        data_frommap_selvar = data_frommap_merged_sel[varname]
        if dimn_layer in data_frommap_selvar.dims:
            crs_plotdata = data_frommap_selvar.to_numpy().T.flatten()
        else: #for 2D models, no layers
            crs_plotdata = data_frommap_selvar.to_numpy()
        return crs_verts, crs_plotdata
    else:
        return crs_verts


def get_netdata(file_nc, multipart=None):

    warnings.warn(DeprecationWarning('dfm_tools.get_nc.get_netdata() will be deprecated, since there is an xarray alternative for multidomain FM files (xugrid). Check the example scripts for how to use it'))
    file_ncs = get_ncfilelist(file_nc, multipart)
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
    warnings.warn(PendingDeprecationWarning('dfm_tools.get_nc.plot_netmapdata() will be deprecated, since there is an xarray alternative for multidomain FM files (xugrid). Check the example scripts for how to use it'))
    
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


def plot_ztdata(data_xr_sel, varname, ax=None, mask_data=True, only_contour=False, **kwargs):
    """
    

    Parameters
    ----------
    data_xr : TYPE
        DESCRIPTION.
    varname : TYPE
        DESCRIPTION.
    ax : matplotlib.axes._subplots.AxesSubplot, optional
        the figure axis. The default is None.
    mask_data : bool, optional
        whether to repair z_interface coordinates and mask data in inactive layers. The default is True.
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

   
    print('WARNING: layers in dflowfm hisfile might be incorrect, check your figures carefully')
    
    data_fromhis_var = data_xr_sel[varname].to_numpy()
    if len(data_fromhis_var.shape) != 2:
        raise Exception(f'ERROR: unexpected number of dimensions in requested squeezed variable ({data_fromhis_var.shape}), first use data_xr.isel(stations=int) to select a single station') #TODO: can also have a different cause, improve message/testing?
    data_fromhis_zcen = data_xr_sel['zcoordinate_c'].bfill(dim='laydim').to_numpy()
    data_fromhis_zcor = data_xr_sel['zcoordinate_w'].bfill(dim='laydimw').to_numpy() #bfill replaces nan values with last valid value, this is necessary to enable pcolormesh to work
    data_fromhis_zcor = np.concatenate([data_fromhis_zcor,data_fromhis_zcor[[-1],:]],axis=0)
    data_fromhis_wl = data_xr_sel['waterlevel'].to_numpy()
    
    if mask_data:
        data_fromhis_var = np.ma.array(data_fromhis_var)
        bool_zcen_equaltop = (data_fromhis_zcen==data_fromhis_zcen[:,-1:]).all(axis=0)
        id_zcentop = np.argmax(bool_zcen_equaltop) # id of first z_center that is equal to z_center of last layer
        if (data_fromhis_zcor[:-1,id_zcentop] > data_fromhis_zcen[:,id_zcentop]).any():
            print('correcting z interface values')
            data_fromhis_zcor[:-1,id_zcentop+1] = data_fromhis_wl
            data_fromhis_zcor[:-1,id_zcentop] = (data_fromhis_zcen[:,id_zcentop-1]+data_fromhis_zcen[:,id_zcentop])/2
        bool_zcen_equaltop[id_zcentop] = False
        #bool_zcor_equaltop = (data_fromhis_zcor_flat[:,1:]==data_fromhis_zcor_flat[:,-1:]).all(axis=0)
        mask_array = np.tile(bool_zcen_equaltop,(data_fromhis_zcor.shape[0],1))
        data_fromhis_var.mask = mask_array

    if not ax: ax=plt.gca()

    # generate 2 2d grids for the x & y bounds (you can also give one 2D array as input in case of eg time varying z coordinates)
    time_np = data_xr_sel.time.to_numpy()
    time_cor = np.concatenate([time_np,time_np[[-1]]])
    time_mesh_cor = np.tile(time_cor,(data_fromhis_zcor.shape[-1],1)).T
    time_mesh_cen = np.tile(time_np,(data_fromhis_zcen.shape[-1],1)).T
    if only_contour:
        pc = ax.contour(time_mesh_cen,data_fromhis_zcen,data_fromhis_var, **kwargs)
    else: #TODO: should actually supply cell edges instead of centers to pcolor/pcolormesh, but inconvenient for time dimension.
        pc = ax.pcolormesh(time_mesh_cor, data_fromhis_zcor, data_fromhis_var, **kwargs)

    return pc

