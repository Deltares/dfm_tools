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

Created on Fri Feb 14 12:43:19 2020

@author: veenstra

helper functions for functions in get_nc.py
"""


import xarray as xr
import pandas as pd
import re
import glob
import os
import warnings
import numpy as np
from netCDF4 import Dataset
from dfm_tools.xarray_helpers import preprocess_hisnc


def get_ncfilelist(file_nc, multipart=None):
    if '' in file_nc:
        file_ncs = [file_nc]
        return file_ncs
    #get list of mapfiles
    
    if not os.path.exists(file_nc):
        raise Exception('ERROR: file does not exist: %s'%(file_nc))
    
    if '_' in file_nc:
        nctype = file_nc.split('_')[-1]
        if nctype == 'rst.nc' and len(file_nc.split('_')) >= 4:
            lastpart = file_nc.split('_')[-4]
        else:
            lastpart = file_nc.split('_')[-2]
        if file_nc.endswith('_%s'%(nctype)) and multipart != False and len(lastpart) == 4 and lastpart.isdigit(): #if part before '_map.nc' is eg '0000'
            if nctype == 'rst.nc' and len(file_nc.split('_')) >= 4:
                filename_start = re.compile('(.*)_([0-9]+)_(.*)_(.*)_%s'%(nctype)).search(file_nc).group(1)
            else:
                filename_start = re.compile('(.*)_([0-9]+)_%s'%(nctype)).search(file_nc).group(1)
            #filename_number = re.compile('(.*)_([0-9]+)_map.nc').search(file_nc).group(2)
            #file_ncs = [file_nc.replace('_%s_map.nc','_%04d_map.nc'%(filename_number, domain_id)) for domain_id in range(ndomains)]
            file_ncs = glob.glob('%s_*_%s'%(filename_start,nctype))
            filename_merged = '%s_merged_%s'%(filename_start,nctype)
            if filename_merged in file_ncs:
                file_ncs.remove(filename_merged)
            
        else:
            file_ncs = [file_nc]
    else:
        file_ncs = [file_nc]
    return file_ncs


def get_varname_fromnc(data_nc,varname_requested,vardim):
    #TODO: put this translationtable in preprocess function, optionally give that to xarray. Raise exception when eg plotnetmapdata sees old variables, saying you should use the preprocess func in xr.open_dataset()
    #VARIABLE names used within different versions of Delft3D-Flexible Mesh
    varnames_list = pd.DataFrame()
    #varnames_list['time'] = ['time','nmesh2d_dlwq_time','TIME','','',''] # time, not necessary anymore
    
    if vardim == 'var':
        varnames_list['mesh2d_node_x'] = ['mesh2d_node_x','NetNode_x','mesh2d_agg_node_x','','',''] # x-coordinate of nodes
        varnames_list['mesh2d_node_y'] = ['mesh2d_node_y','NetNode_y','mesh2d_agg_node_y','','',''] # y-coordinate of nodes
        varnames_list['mesh2d_node_z'] = ['mesh2d_node_z','NetNode_z','','','',''] # z-coordinate of nodes
        
        varnames_list['mesh2d_face_x'] = ['mesh2d_face_x','FlowElem_xzw','mesh2d_agg_face_x','','',''] # x-coordinate of faces (center)
        varnames_list['mesh2d_face_y'] = ['mesh2d_face_y','FlowElem_yzw','mesh2d_agg_face_y','','',''] # y-coordinate of faces (center)
        
        varnames_list['mesh2d_edge_x'] = ['mesh2d_edge_x','','','','',''] # x-coordinate of velocity-points
        varnames_list['mesh2d_edge_y'] = ['mesh2d_edge_y','','','','',''] # y-coordinate of velocity-points
        
        varnames_list['mesh2d_edge_nodes'] = ['mesh2d_edge_nodes','NetLink','','','',''] # 'link between two netnodes' / 'Mapping from every edge to the two nodes that it connects'
        varnames_list['mesh2d_edge_faces'] = ['mesh2d_edge_faces','','','','',''] # 'Neighboring faces of mesh edges'
        varnames_list['mesh2d_face_nodes'] = ['mesh2d_face_nodes','NetElemNode','mesh2d_agg_face_nodes','','',''] # 
        
        varnames_list['mesh2d_face_x_bnd'] = ['mesh2d_face_x_bnd','FlowElemContour_x','mesh2d_agg_face_x_bnd','','',''] # x-coordinates of flow element contours
        varnames_list['mesh2d_face_y_bnd'] = ['mesh2d_face_y_bnd','FlowElemContour_y','mesh2d_agg_face_y_bnd','','',''] # y-coordinates of flow element contours
        
        varnames_list['mesh2d_flowelem_domain'] = ['mesh2d_flowelem_domain','FlowElemDomain','','','',''] # flow element domain
        varnames_list['mesh2d_flowelem_bl'] = ['mesh2d_flowelem_bl','FlowElem_bl','','','',''] # bed level
        varnames_list['mesh2d_flowelem_ba'] = ['mesh2d_flowelem_ba','FlowElem_bac','','','',''] # area (m2) of cell faces
        varnames_list['mesh2d_s1'] = ['mesh2d_s1','','','','',''] # water level
    
        #varnames_list['mesh2d_ucx'] = ['mesh2d_ucx','ucx','',''] # 
        #varnames_list['mesh2d_ucy'] = ['mesh2d_ucy','ucy','',''] # 
        #varnames_list['mesh2d_sa1'] = ['mesh2d_sa1','sa1','',''] # 
        #varnames_list['mesh2d_tem1'] = ['mesh2d_tem1','tem1','',''] # 
    
    elif vardim == 'dim':
        ### DIMENSION names used within different versions of Delft3D-Flexible Mesh
        #dimnames_list = pd.DataFrame()
        varnames_list['nmesh2d_node'] = ['nmesh2d_node','mesh2d_nNodes','nNetNode','','','',''] # number of nodes
        varnames_list['nmesh2d_face'] = ['nmesh2d_face','mesh2d_nFaces','nNetElem','','','',''] # number of faces
        varnames_list['nmesh2d_edge'] = ['nmesh2d_edge','mesh2d_nEdges','nNetLink','','','',''] # number of velocity-points
        varnames_list['nFlowElem'] = ['nFlowElem','','','','','',''] # number of flow elements
        varnames_list['nFlowLink'] = ['nFlowLink','','','','','',''] # number of flow elements
        
        varnames_list['nmesh2d_layer'] = ['nmesh2d_layer','mesh2d_nLayers','laydim','nmesh2d_layer_dlwq','LAYER','KMAXOUT_RESTR','depth'] # layer
    else:
        raise Exception('parameter vardim can be "var" or "dim"')
    
    #look for correct pd column
    pdcol_bool = varnames_list.eq(varname_requested).any()
    varname_pdcol = pdcol_bool.index[pdcol_bool].tolist()
    if len(varname_pdcol) == 0:
        raise Exception('varname %s not found in internal database'%(varname_requested))
    elif len(varname_pdcol)>1:
        raise Exception('varname %s not found but multiple equivalents found in internal database: %s'%(varname_requested,varname_pdcol))
    else:
        varname_pdcol = varname_pdcol[0]
    
    if vardim == 'var':
        data_nc_vardimnames_list = list(data_nc.variables.keys())
    elif vardim == 'dim':
        data_nc_vardimnames_list = list(data_nc.dimensions.keys())
    else:
        raise Exception('parameter vardim can be "var" or "dim"')
    
    def get_vardimname(data_nc_names_list):
        #check what is in netcdf file
        if varname_requested in data_nc_names_list:
            varname = varname_requested
        elif varname_pdcol in data_nc_names_list:
            varname = varname_pdcol
        else:
            var_options = list(varnames_list[varname_pdcol])
            varname = [var for var in var_options if var in data_nc_names_list]
            if varname == []:
                varname = None
            else:
                varname = varname[0]
        return varname
    
    varname = get_vardimname(data_nc_vardimnames_list)
    
    return varname


def get_ncvarproperties(file_nc):
    data_xr = xr.open_dataset(file_nc)
    nc_varkeys = data_xr.variables.mapping.keys()
    
    list_varattrs_pd = []
    for varkey in nc_varkeys:
        varattrs_pd = pd.DataFrame({varkey:data_xr.variables.mapping[varkey].attrs}).T
        varattrs_pd[['shape','dimensions']] = 2*[''] #set dtype as str (float will raise an error when putting tuple in there)
        varattrs_pd.at[varkey,'shape'] = data_xr[varkey].shape
        varattrs_pd.at[varkey,'dimensions'] = data_xr.variables[varkey].dims
        varattrs_pd.loc[varkey,'dtype'] = data_xr.variables[varkey].dtype
        list_varattrs_pd.append(varattrs_pd)
    
    vars_pd = pd.concat(list_varattrs_pd,axis=0)
    vars_pd[vars_pd.isnull()] = '' #avoid nan values
    
    data_xr.close()

    return vars_pd


def get_ncvardimlist(file_nc):
    raise DeprecationWarning('use dfm_tools.get_nc_helpers.get_ncvarproperties() instead') #TODO: remove this code
    vars_pd = get_ncvarproperties(file_nc)
    
    return vars_pd, None


def get_varnamefrom_keyslongstandardname(file_nc, varname):
    DeprecationWarning('use dfm_tools.get_nc_helpers.get_varnamefromattrs() instead') #TODO: raise this warning, later remove this code
    varname_matched = get_varnamefromattrs(file_nc, varname)
    return varname_matched


def get_varnamefromattrs(file_nc, varname):
    data_xr = xr.open_dataset(file_nc)
    
    # check if requested variable is in netcdf
    varlist = list(data_xr.variables.keys())
    if varname in varlist:
        return varname
    
    #check if requested varname is in standard_name attrs of ncvars
    ds_stdname = data_xr.filter_by_attrs(standard_name=varname)
    varlist_stdname = list(ds_stdname.data_vars.keys())
    if len(varlist_stdname)==1:
        varname_matched = varlist_stdname[0]
        print(f'requested varname "{varname}" found in standard_name attribute of variable {varname_matched}')
        return varname_matched
    elif len(varlist_stdname)>1:
        raise Exception(f'ERROR: requested variable {varname} is in netcdf not 1 but {len(varlist_stdname)} times: {varlist_stdname}')
    
    #check if requested varname is in long_name attrs of ncvars
    ds_longname = data_xr.filter_by_attrs(long_name=varname)
    varlist_longname = list(ds_longname.data_vars.keys())
    if len(varlist_longname)==1:
        varname_matched = varlist_longname[0]
        print(f'requested varname "{varname}" found in long_name attribute of variable {varname_matched}')
        return varname_matched
    elif len(varlist_longname)>1:
        raise Exception(f'ERROR: requested variable {varname} is in netcdf not 1 but {len(varlist_longname)} times: {varlist_longname}')
    
    #if not returned above, the varname was not found so raise exception
    raise Exception(f'ERROR: requested variable {varname} not in netcdf, available are: {varlist} and the standard_name and long_name attrs in dfm_tools.get_nc_helpers.get_ncvarproperties(file_nc=file_nc)')    
    return varname_matched


def ghostcell_filter(file_nc):
        
    data_nc = Dataset(file_nc)
    
    varn_domain = get_varname_fromnc(data_nc,'mesh2d_flowelem_domain',vardim='var')
    if varn_domain is not None: # domain variable is present, so there are multiple domains
        domain = data_nc.variables[varn_domain][:]
        domain_no = np.bincount(domain).argmax() #meest voorkomende domeinnummer
        nonghost_bool = domain==domain_no
    else:
        nonghost_bool = None
        
    data_nc.close()
    return nonghost_bool


def get_timesfromnc(file_nc, varname='time', retrieve_ids=False, keeptimezone=True):
    """
    Retrieves time array from netcdf file. Time is converted to UTC by default, but times can optionally be returned in the original timezone.

    Parameters
    ----------
    file_nc : STR
        DESCRIPTION.
    varname : STR, optional
        DESCRIPTION. The default is 'time'.
    retrieve_ids : LIST of int, optional
        DESCRIPTION. The default is False.
    keeptimezone : BOOL, optional
        DESCRIPTION. The default is True.

    Raises
    ------
    Exception
        DESCRIPTION.

    Returns
    -------
    data_nc_datetimes_pd : TYPE
        DESCRIPTION.

    """
    
    with xr.open_dataset(file_nc) as data_xr:
        #varn_time = get_variable_timevar(file_nc,varname=varname)
        #times_xr = data_xr[varn_time]
        times_xr = data_xr[varname].time
        times_pd = times_xr.to_series().dt.round(freq='S')
    
    if len(times_xr.shape) > 1:
        raise Exception('Corrupt netCDF files with a time variable with more than one dimension')
    
    #convert back to original timezone (e.g. MET)
    if keeptimezone:
        time_units = data_xr.time.encoding['units']
        time_units_tz = pd.DatetimeIndex([time_units.split('since ')[1]]).tz
        times_pd_local = times_pd.tz_localize('UTC').tz_convert(time_units_tz)
        times_pd = pd.Series(times_pd_local.index,index=times_pd_local.index)
        
    if retrieve_ids:
        times_pd = times_pd.iloc[retrieve_ids]
    
    return times_pd


def get_timeid_fromdatetime(data_nc_datetimes_pd, timestep):
    
    timestep_pd = pd.Series(timestep)

    #check if all requested times (timestep) are in netcdf file
    times_bool_reqinfile = timestep_pd.isin(data_nc_datetimes_pd)
    if not (times_bool_reqinfile == True).all():
        raise Exception('ERROR: not all requested times are in netcdf file:\n%s\navailable in netcdf file are:\n%s\nUse this command to obtain full list as variable:\nfrom dfm_tools.get_nc_helpers import get_timesfromnc\ndata_nc_datetimes_pd = get_timesfromnc(file_nc=file_nc)'%(timestep_pd[-times_bool_reqinfile], data_nc_datetimes_pd))
        
    #get ids of requested times in netcdf file
    times_bool_fileinreq = data_nc_datetimes_pd.isin(timestep_pd)
    time_ids = np.where(times_bool_fileinreq)[0]
    
    return time_ids


def get_hisstationlist(file_nc, varname='waterlevel'):
    warnings.warn(DeprecationWarning('use data_xr[\'stations\'].to_dataframe() instead, do read in your hisfile with preprocess=preprocess_hisnc as argument, like in the postprocess_gethismodeldata.py example script'))
    data_xr = xr.open_mfdataset(file_nc)#, preprocess=preprocess_hisnc)
    
    varname = get_varnamefrom_keyslongstandardname(file_nc, varname) #get varname from varkeys/standardname/longname if exists
    vardims = data_xr[varname].dims
    
    vars_pd = get_ncvarproperties(file_nc=file_nc)
    bool_vars_dtypestr = vars_pd['dtype'].astype(str).str.startswith('|S') | (vars_pd['dtype']=='object')#& (vars_pd['ndims']==1) #TODO: better check for dtype string?
    vars_pd_stats = vars_pd.loc[bool_vars_dtypestr]
    
    if len(vars_pd_stats) == 0:
        raise Exception('ERROR: no dimension in %s variable that corresponds to station-like variables (or none present):\n%s'%(varname, vars_pd_stats.index))
        
    dimkey_use = None
    for dimtuple in vars_pd_stats['dimensions']:
        dimkey = dimtuple[0]
        if dimkey in vardims:
            dimkey_use = dimkey
    
    statlist_pd = data_xr[dimkey_use].to_dataframe()
    
    for colname in statlist_pd.columns:
        if isinstance(statlist_pd[colname][0],bytes): #TODO: better check would be data_xr.variables[colname].dtype=='S256'
            print(f'variable {colname}: converting bytes to str')
            statlist_pd[colname] = data_xr[colname].to_pandas().str.decode('utf-8',errors='ignore').str.strip() #to_pandas is essential otherwise resulting bool might not be correct. .str.strip() to remove spaces left/right from station_name (necessary for sobek models)
    
    data_xr.close()
    return statlist_pd


def get_stationid_fromstationlist(data_xr, stationlist, station_varname='station_name'): #TODO: this can be removed soon since it is only used in get_nc
    warnings.warn(DeprecationWarning('Do not use dfm_tools.get_nc_helpers.get_stationid_fromstationlist() in your script, it is only used internally and will be phased out asap. use xarray ds.sel() instead, but read in your hisfile with preprocess=preprocess_hisnc like in the examplescript postprocess_gethismodeldata.py'))
    if not isinstance(stationlist,list):
        raise Exception('ERROR: provide list of stations')
    
    stationlist_pd = pd.Series(stationlist)
    data_xr_stationlist_pd = data_xr[station_varname].to_pandas().str.decode('utf-8',errors='ignore').str.strip() #to_pandas is essential otherwise resulting bool might not be correct. .str.strip() to remove spaces left/right from station_name (necessary for sobek models)
    data_xr_stationlist_pd = data_xr_stationlist_pd.reset_index(drop=True) #in case the index contains stationnames instead of numbers
    
    #check if all requested stations are in xarray Dataset
    bool_reqstations = stationlist_pd.isin(data_xr_stationlist_pd)
    if not bool_reqstations.all():
        print(data_xr_stationlist_pd)
        raise Exception(f'ERROR: not all requested stations in netcdf:\n{stationlist_pd.loc[~bool_reqstations]}')
    
    #get idx of requested stations, in original order
    idx_stations = [data_xr_stationlist_pd.loc[data_xr_stationlist_pd.str.match(stat_req)].index[0] for stat_req in stationlist]
    
    return idx_stations



