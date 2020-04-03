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

Created on Fri Feb 14 12:43:19 2020

@author: veenstra

helper functions for functions in get_nc.py
"""

    
    

def get_ncfilelist(file_nc, multipart=None):
    #get list of mapfiles
    import re
    import glob
    import os
    
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
            file_ncs = glob.glob('%s*_%s'%(filename_start,nctype))
        else:
            file_ncs = [file_nc]
    else:
        file_ncs = [file_nc]
    return file_ncs




def get_varname_fromnc(data_nc,varname_requested):
    import pandas as pd
    
    #VARIABLE names used within different versions of Delft3D-Flexible Mesh
    varnames_list = pd.DataFrame()
    varnames_list['time'] = ['time','nmesh2d_dlwq_time','TIME','','',''] # time
    
    varnames_list['mesh2d_node_x'] = ['mesh2d_node_x','NetNode_x','mesh2d_agg_node_x','','',''] # x-coordinate of nodes
    varnames_list['mesh2d_node_y'] = ['mesh2d_node_y','NetNode_y','mesh2d_agg_node_y','','',''] # y-coordinate of nodes
    varnames_list['mesh2d_node_z'] = ['mesh2d_node_z','NetNode_z','','','',''] # z-coordinate of nodes
    
    varnames_list['mesh2d_face_x'] = ['mesh2d_face_x','FlowElem_xzw','mesh2d_agg_face_x','','',''] # x-coordinate of faces
    varnames_list['mesh2d_face_y'] = ['mesh2d_face_y','FlowElem_yzw','mesh2d_agg_face_y','','',''] # y-coordinate of faces
    
    varnames_list['mesh2d_edge_x'] = ['mesh2d_edge_x','','','','',''] # x-coordinate of velocity-points
    varnames_list['mesh2d_edge_y'] = ['mesh2d_edge_y','','','','',''] # y-coordinate of velocity-points
    
    varnames_list['mesh2d_edge_nodes'] = ['mesh2d_edge_nodes','NetLink','','','',''] # 'link between two netnodes' / 'Mapping from every edge to the two nodes that it connects'
    varnames_list['mesh2d_face_nodes'] = ['mesh2d_face_nodes','NetElemNode','mesh2d_agg_face_nodes','','',''] # 
    
    varnames_list['mesh2d_face_x_bnd'] = ['mesh2d_face_x_bnd','FlowElemContour_x','mesh2d_agg_face_x_bnd','','',''] # x-coordinates of flow element contours
    varnames_list['mesh2d_face_y_bnd'] = ['mesh2d_face_y_bnd','FlowElemContour_y','mesh2d_agg_face_y_bnd','','',''] # y-coordinates of flow element contours
    
    varnames_list['mesh2d_flowelem_domain'] = ['mesh2d_flowelem_domain','FlowElemDomain','','','',''] # flow element domain
    varnames_list['mesh2d_flowelem_bl'] = ['mesh2d_flowelem_bl','FlowElem_bl','','','',''] # bed level
    varnames_list['mesh2d_flowelem_ba'] = ['mesh2d_flowelem_ba','FlowElem_bac','','','',''] # area (m2) of cell faces
    
    varnames_list['mesh2d_layer_z'] = ['mesh2d_layer_z','LayCoord_cc','','','',''] # 
    
    #non-grid variables necessary for layer calculation for intersection/cross section) funtion
    varnames_list['mesh2d_s1'] = ['mesh2d_s1','','','','',''] # water level
    varnames_list['mesh2d_flowelem_bl'] = ['mesh2d_flowelem_bl','','','','',''] # bed level

    #varnames_list['mesh2d_ucx'] = ['mesh2d_ucx','ucx','',''] # 
    #varnames_list['mesh2d_ucy'] = ['mesh2d_ucy','ucy','',''] # 
    #varnames_list['mesh2d_sa1'] = ['mesh2d_sa1','sa1','',''] # 
    #varnames_list['mesh2d_tem1'] = ['mesh2d_tem1','tem1','',''] # 
    
    
    ### DIMENSION names used within different versions of Delft3D-Flexible Mesh
    #dimnames_list = pd.DataFrame()
    varnames_list['nmesh2d_node'] = ['nmesh2d_node','mesh2d_nNodes','nNetNode','','',''] # number of nodes
    varnames_list['nmesh2d_face'] = ['nmesh2d_face','mesh2d_nFaces','nNetElem','','',''] # number of faces
    varnames_list['nmesh2d_edge'] = ['nmesh2d_edge','nNetLink','','','',''] # number of velocity-points
    varnames_list['nFlowElem'] = ['nFlowElem','','','','',''] # number of flow elements
    varnames_list['nFlowLink'] = ['nFlowLink','','','','',''] # number of flow elements
    
    varnames_list['nmesh2d_layer'] = ['nmesh2d_layer','mesh2d_nLayers','laydim','nmesh2d_layer_dlwq','LAYER','KMAXOUT_RESTR'] # layer
    
    #look for correct pd column
    pdcol_bool = varnames_list.eq(varname_requested).any()
    varname_pdcol = pdcol_bool.index[pdcol_bool].tolist()
    if len(varname_pdcol) == 0:
        raise Exception('varname %s not found in internal database'%(varname_requested))
    elif len(varname_pdcol)>1:
        raise Exception('varname %s not found but multiple equivalents found in internal database: %s'%(varname_requested,varname_pdcol))
    else:
        varname_pdcol = varname_pdcol[0]
    
    data_nc_varnames_list = list(data_nc.variables.keys())
    data_nc_dimnames_list = list(data_nc.dimensions.keys())
    
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
    
    varname = get_vardimname(data_nc_varnames_list)
    if varname is None:
        varname = get_vardimname(data_nc_dimnames_list)
    
    return varname






def get_ncvardimlist(file_nc):
    from netCDF4 import Dataset
    import pandas as pd
    import numpy as np
    
    data_nc = Dataset(file_nc)
    
    nc_varkeys = list(data_nc.variables.keys())
    nc_dimkeys = list(data_nc.dimensions.keys())
    vars_pd = pd.DataFrame({'nc_varkeys': nc_varkeys})
    dims_pd = pd.DataFrame({'nc_dimkeys': nc_dimkeys})
    var_attr_name_list = ['standard_name','long_name','coordinates','units','mesh','location']
    for iV, nc_var in enumerate(data_nc.variables):
        #get non-attribute properties of netcdf variable
        if iV==0:
            vars_pd['shape'] = np.nan
            vars_pd['dimensions'] = np.nan
        vars_pd.loc[iV,'shape'] = str(data_nc.variables[nc_var].shape)
        vars_pd.loc[iV,'dimensions'] = str(data_nc.variables[nc_var].dimensions)
        #get attributes properties of netcdf variable
        for attr_name in var_attr_name_list:
            if iV==0:
                vars_pd[attr_name] = ''
            try:
                vars_pd.loc[iV, attr_name] = data_nc.variables[nc_var].getncattr(attr_name)
            except:
                pass
    for iD, nc_dim in enumerate(data_nc.dimensions):
        #get non-attribute properties of netcdf variable
        if iD==0:
            dims_pd['name'] = np.nan
            dims_pd['size'] = np.nan
        dims_pd.loc[iD,'name'] = data_nc.dimensions[nc_dim].name
        dims_pd.loc[iD,'size'] = data_nc.dimensions[nc_dim].size
    return vars_pd, dims_pd

    


def get_ncvarobject(file_nc, varname):
    from netCDF4 import Dataset
    
    vars_pd, dims_pd = get_ncvardimlist(file_nc=file_nc)
    nc_varkeys = list(vars_pd['nc_varkeys'])
    
    # check if requested variable is in netcdf
    if varname not in nc_varkeys:
        raise Exception('ERROR: requested variable %s not in netcdf, available are:\n%s\nUse this command to obtain full list as variable:\nfrom dfm_tools.get_nc_helpers import get_ncvardimlist; vars_pd, dims_pd = get_ncvardimlist(file_nc=file_nc)'%(varname, vars_pd))
    
    data_nc = Dataset(file_nc)
    nc_varobject = data_nc.variables[varname]
    
    return nc_varobject



def ghostcell_filter(file_nc):
    import numpy as np
    from netCDF4 import Dataset
    
    #from dfm_tools.get_nc_helpers import get_varname_fromnc
    
    data_nc = Dataset(file_nc)
    
    varn_domain = get_varname_fromnc(data_nc,'mesh2d_flowelem_domain')
    if varn_domain is not None: # domain variable is present, so there are multiple domains
        domain = data_nc.variables[varn_domain][:]
        domain_no = np.bincount(domain).argmax() #meest voorkomende domeinnummer
        nonghost_ids = domain==domain_no
    else:
        nonghost_ids = None
    return nonghost_ids



def get_timesfromnc(file_nc, force_noreconstruct=False):
    """
    retrieves time array from netcdf file.
    Since long time arrays take a long time to retrieve at once, reconstruction is tried
    in dflowfm an array can start with 0 (initial), followed by a tstart and increading with intervals to tend
    therefore, the interval at the start and end of the time array is not always equal to the 'real' time interval
    reconstruction takes care of this.
    if you still feel like the resulting timeseries is incorrect, use the keyword force_noreconstruct
    """
    
    from netCDF4 import Dataset,num2date#,date2num
    import numpy as np
    import pandas as pd
    
    #from dfm_tools.get_nc_helpers import get_varname_fromnc

    data_nc = Dataset(file_nc)
    varname_time = get_varname_fromnc(data_nc,'time')
    data_nc_timevar = data_nc.variables[varname_time]
    
    if len(data_nc_timevar)<3 or force_noreconstruct==True: #shorter than 3 rarely is the case, but just to be sure
        data_nc_times = data_nc_timevar[:]
        print('reading time dimension: read entire array (len < 3 or force_noreconstruct==True)')
    else:
        time0 = data_nc_timevar[0] 
        time1 = data_nc_timevar[1] 
        time2 = data_nc_timevar[2]
        timemin3 = data_nc_timevar[-3]
        timemin2 = data_nc_timevar[-2]
        timemin1 = data_nc_timevar[-1]
        timeinc_poststart = time2-time1 # the interval between 0 and 1 is not per definition representative, so take 1 and 2
        timeinc_preend = timemin2-timemin3
        #timeinc_end = timemin1-timemin2
        if timeinc_poststart == timeinc_preend: #reconstruct time array to save time
            print('reading time dimension: reconstruct array')
            data_nc_times_from1 = np.arange(time1,timemin1,timeinc_poststart)
            data_nc_times = np.concatenate([[time0],data_nc_times_from1,[timemin1]])
        else:
            print('reading time dimension: read entire array')
            data_nc_times = data_nc_timevar[:]
            
    data_nc_datetimes = num2date(data_nc_times, units = data_nc_timevar.units)
    #data_nc_datetimes_pd = pd.Series(data_nc_datetimes).dt.round(freq='S')
    nptimes = data_nc_datetimes.astype('datetime64[ns]') #convert to numpy first, pandas does not take all cftime datasets
    data_nc_datetimes_pd = pd.Series(nptimes).dt.round(freq='S')
    
    return data_nc_datetimes_pd




def get_timeid_fromdatetime(data_nc_datetimes_pd, timestep):
    import numpy as np
    import pandas as pd
    
    timestep_pd = pd.Series(timestep)#.dt.round(freq='S')

    #check if all requested times (timestep) are in netcdf file
    times_bool_reqinfile = timestep_pd.isin(data_nc_datetimes_pd)
    if not (times_bool_reqinfile == True).all():
        raise Exception('ERROR: not all requested times are in netcdf file:\n%s\navailable in netcdf file are:\n%s\nUse this command to obtain full list as variable:\nfrom dfm_tools.get_nc_helpers import get_timesfromnc; data_nc_datetimes_pd = get_timesfromnc(file_nc=file_nc)'%(timestep_pd[-times_bool_reqinfile], data_nc_datetimes_pd))
        
    #get ids of requested times in netcdf file
    times_bool_fileinreq = data_nc_datetimes_pd.isin(timestep_pd)
    time_ids = np.where(times_bool_fileinreq)[0]
    
    return time_ids






def get_vardimname_stat_validvals():
    #variable/dim names for:   DFM stations,   DFM gs,                 DFM crs,              Sobek stations,    WAQUA_getdata_netcdf WL/CUR-stations,  Delft3D netCDF stations
    varname_stat_validvals = ['station_name', 'general_structure_id', 'cross_section_name', 'observation_id',  'NAMWL',   'NAMC',                      'NAMST']
    dimname_stat_validvals = ['stations',     'general_structures',   'cross_section',      'id',              'STATION', 'STATIONCUR',                'NOSTAT']
    return varname_stat_validvals, dimname_stat_validvals




def get_hisstationlist(file_nc,varname_stat='station_name'):
    from netCDF4 import Dataset, chartostring
    import pandas as pd
    import numpy as np
    
    from dfm_tools.get_nc_helpers import get_vardimname_stat_validvals
    
    varname_stat_validvals, dimname_stat_validvals = get_vardimname_stat_validvals()
    if varname_stat in varname_stat_validvals:
        data_nc = Dataset(file_nc)
        station_name_char = data_nc.variables[varname_stat][:]
        station_name_list_raw = chartostring(station_name_char)
        station_name_list = np.char.strip(station_name_list_raw) #necessary step for Sobek and maybe others
        
        station_name_list_pd = pd.Series(station_name_list)
    else:
        raise Exception('ERROR: invalid value provided for varname_stat argument (%s), should be one of: %s'%(varname_stat, varname_stat_validvals))
    return station_name_list_pd




def get_stationid_fromstationlist(station_name_list_pd, station, varname_stat):
    import numpy as np
    import pandas as pd
    
    station_pd = pd.Series(station)

    #check if all requested stations are in netcdf file
    stations_bool_reqinfile = station_pd.isin(station_name_list_pd)
    if not (stations_bool_reqinfile == True).all():
        raise Exception('ERROR: not all requested stations are in netcdf file:\n%s\navailable in netcdf file are:\n%s\nUse this command to obtain full list as variable:\nfrom dfm_tools.get_nc_helpers import get_hisstationlist; station_name_list_pd = get_hisstationlist(file_nc=file_nc,varname_stat="%s")'%(station_pd[-stations_bool_reqinfile], station_name_list_pd, varname_stat))
    
    #get ids of requested stations in netcdf file
    station_bool_fileinreq = station_name_list_pd.isin(station_pd)
    station_ids = list(np.where(station_bool_fileinreq)[0])
    #station_names = np.where(station_bool_fileinreq)[0]

    return station_ids





