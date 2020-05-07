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
            filename_merged = '%s_merged_%s'%(filename_start,nctype)
            if filename_merged in file_ncs:
                file_ncs.remove(filename_merged)
            
        else:
            file_ncs = [file_nc]
    else:
        file_ncs = [file_nc]
    return file_ncs




def get_varname_fromnc(data_nc,varname_requested,vardim):
    import pandas as pd
    
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
    
    elif vardim == 'dim':
        ### DIMENSION names used within different versions of Delft3D-Flexible Mesh
        #dimnames_list = pd.DataFrame()
        varnames_list['nmesh2d_node'] = ['nmesh2d_node','mesh2d_nNodes','nNetNode','','','',''] # number of nodes
        varnames_list['nmesh2d_face'] = ['nmesh2d_face','mesh2d_nFaces','nNetElem','','','',''] # number of faces
        varnames_list['nmesh2d_edge'] = ['nmesh2d_edge','nNetLink','','','','',''] # number of velocity-points
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
    #if varname is None:
    #    varname = get_vardimname(data_nc_dimnames_list)
    
    return varname






def get_ncvardimlist(file_nc):
    from netCDF4 import Dataset
    import pandas as pd
    #import numpy as np
    
    data_nc = Dataset(file_nc)
    
    nc_varkeys = list(data_nc.variables.keys())
    nc_dimkeys = list(data_nc.dimensions.keys())
    vars_pd = pd.DataFrame({'nc_varkeys': nc_varkeys, 'shape': [['']]*len(nc_varkeys), 'dimensions': [['']]*len(nc_varkeys), 'dtype': [['']]*len(nc_varkeys)})
    dims_pd = pd.DataFrame({'nc_dimkeys': nc_dimkeys, 'name': [['']]*len(nc_dimkeys), 'size': [['']]*len(nc_dimkeys)})
    var_attr_name_list = ['standard_name','long_name','coordinates','units','mesh','location']
    for iV, nc_var in enumerate(data_nc.variables):
        #get non-attribute properties of netcdf variable
        vars_pd.loc[iV,'shape'] = data_nc.variables[nc_var].shape
        vars_pd.loc[iV,'dimensions'] = data_nc.variables[nc_var].dimensions
        vars_pd.loc[iV,'dtype'] = data_nc.variables[nc_var].dtype
        #get attributes properties of netcdf variable
        for attr_name in var_attr_name_list:
            if iV==0:
                vars_pd[attr_name] = ''
            try:
                vars_pd.loc[iV, attr_name] = data_nc.variables[nc_var].getncattr(attr_name)
            except:
                pass
    vars_pd['ndims'] = [len(x) for x in vars_pd['dimensions']]
    for iD, nc_dim in enumerate(data_nc.dimensions):
        #get non-attribute properties of netcdf variable
        dims_pd.loc[iD,'name'] = data_nc.dimensions[nc_dim].name
        dims_pd.loc[iD,'size'] = data_nc.dimensions[nc_dim].size
    return vars_pd, dims_pd

    


def get_ncvarobject(file_nc, varname):
    from netCDF4 import Dataset
    
    vars_pd, dims_pd = get_ncvardimlist(file_nc=file_nc)
    
    nc_varkeys = list(vars_pd['nc_varkeys'])
    nc_varlongnames = list(vars_pd['long_name'])
    nc_varstandardnames = list(vars_pd['standard_name'])
    
    # check if requested variable is in netcdf
    if varname == '':
        varname = None
    
    if varname in nc_varkeys:
        pass
    elif varname in nc_varlongnames:
        varid = nc_varlongnames.index(varname)
        varname = vars_pd.loc[varid,'nc_varkeys']
        print('varname found in long_name attribute')
    elif varname in nc_varstandardnames:
        varid = nc_varstandardnames.index(varname)
        varname = vars_pd.loc[varid,'nc_varkeys']
        print('varname found in standard_name attribute')
    else:
        raise Exception('ERROR: requested variable %s not in netcdf, available are:\n%s\nUse this command to obtain full list as variable:\nfrom dfm_tools.get_nc_helpers import get_ncvardimlist; vars_pd, dims_pd = get_ncvardimlist(file_nc=file_nc)\nnote that you can retrieve variables by keys, standard_name or long_name attributes'%(varname, vars_pd))
    
    data_nc = Dataset(file_nc)
    nc_varobject = data_nc.variables[varname]
    
    return nc_varobject



def ghostcell_filter(file_nc):
    import numpy as np
    from netCDF4 import Dataset
    
    #from dfm_tools.get_nc_helpers import get_varname_fromnc
    
    data_nc = Dataset(file_nc)
    
    varn_domain = get_varname_fromnc(data_nc,'mesh2d_flowelem_domain',vardim='var')
    if varn_domain is not None: # domain variable is present, so there are multiple domains
        domain = data_nc.variables[varn_domain][:]
        domain_no = np.bincount(domain).argmax() #meest voorkomende domeinnummer
        nonghost_ids = domain==domain_no
    else:
        nonghost_ids = None
    return nonghost_ids





def get_variable_timevardim(file_nc, varname):
    #get corresponding time variable name
    from netCDF4 import Dataset
    
    #from dfm_tools.get_nc_helpers import get_ncvarobject
    
    data_nc = Dataset(file_nc)
    nc_varobject = get_ncvarobject(file_nc, varname)
    
    varn_time = None
    dimn_time = None
    varlist_wunits = data_nc.get_variables_by_attributes(units=lambda v: v is not None)
    for var_lookup in varlist_wunits:
        if 'since' in var_lookup.units and var_lookup.dimensions[0] in nc_varobject.dimensions:
            dimn_time = var_lookup.dimensions[0]
            varn_time = var_lookup.name
            break
    return varn_time, dimn_time



    
def get_timesfromnc(file_nc, varname='time', retrieve_ids=False, keeptimezone=True):
    """
    retrieves time array from netcdf file.
    Since long time arrays take a long time to retrieve at once, reconstruction is tried
    in dflowfm an array can start with 0 (initial), followed by a tstart and increading with intervals to tend
    therefore, the interval at the start and end of the time array is not always equal to the 'real' time interval
    reconstruction takes care of this.
    if reconstruction fails (the length of the netCDF variable is not equal of the length of the reconstructed array), all times are read

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
    
    from netCDF4 import Dataset, num2date#,date2num
    #from cftime import num2pydate, num2date
    #from cftime import num2date as cf_num2date
    import numpy as np
    import pandas as pd
    import warnings
    import datetime as dt
    
    #from dfm_tools.get_nc_helpers import get_variable_timevardim

    data_nc = Dataset(file_nc)
    varn_time, dimn_time = get_variable_timevardim(file_nc=file_nc, varname=varname)
    data_nc_timevar = data_nc.variables[varn_time]
    time_length = data_nc_timevar.shape[0]

    if retrieve_ids is not False:
        print('reading time dimension: only requested indices')
        listtype_range = [list, range, np.ndarray]
        if type(retrieve_ids) not in listtype_range:
            raise Exception('ERROR: argument retrieve_ids should be a list')
        #convert to positive index, make unique(+sort), convert to list because of indexing with np.array of len 1 errors sometimes
        retrieve_ids = list(np.unique(np.array(range(time_length))[retrieve_ids]))
        data_nc_times = data_nc_timevar[retrieve_ids]
    elif len(data_nc_timevar)<3: #check if time dimension is shorter than 3 items
        data_nc_times = data_nc_timevar[:]
        print('reading time dimension: read entire array (because length < 3)')
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
            if data_nc_timevar.shape[0] != len(data_nc_times):#test if len of reconstructed timeseries is same as len of timevar in netCDF, retrieve entire array
                print('reading time dimension: reconstruction failed, read entire array')
                data_nc_times = data_nc_timevar[:]
        else:
            print('reading time dimension: read entire array')
            data_nc_times = data_nc_timevar[:]
        
    if len(data_nc_times.shape) > 1:
        warnings.warn('This should not happen, this exception is built in for corrupt netCDF files with a time variable with more than one dimension')
        data_nc_times = data_nc_times.flatten()
    
    #convert back to original timezone (e.g. MET)
    if keeptimezone:
        """
        #with timezone info
        #tz_str_startplus = data_nc_timevar.units.rfind('+')
        #tz_str_startmin = data_nc_timevar.units.rfind('-')
        #tz_str_start = np.max([tz_str_startplus, tz_str_startmin])
        tz_str = data_nc_timevar.units.split(' ')[-1]
        try:
            tzoffset = dt.datetime.strptime(tz_str,'%z').utcoffset()
            nptimes = data_nc_datetimes + tzoffset
            print('times converted to original timezone (%s)'%(tzoffset))
        except:
            #print('retrieving original timezone failed, time is now probably UTC')
            nptimes = data_nc_datetimes
        """
        #manual conversion which deliberately ignores timezone
        time_units_list = data_nc_timevar.units.split(' ')
        if time_units_list[1] != 'since':
            raise Exception('invalid time units string (%s)'%(data_nc_timevar.units))
        try:
            refdate_str = '%s %s'%(time_units_list[2], time_units_list[3].replace('.0','')) #remove .0 to avoid conversion issue
            refdate = dt.datetime.strptime(refdate_str,'%Y-%m-%d %H:%M:%S')
            data_nc_times_pdtd = pd.to_timedelta(data_nc_times, unit=time_units_list[0])
            data_nc_datetimes = (refdate + data_nc_times_pdtd)#.to_pydatetime()
            print('retrieving original timezone succeeded, no conversion to UTC/GMT applied')
        except:
            print('retrieving original timezone failed, using num2date output instead')
            data_nc_datetimes = num2date(data_nc_times, units=data_nc_timevar.units, only_use_cftime_datetimes=False, only_use_python_datetimes=True)
    else:
        #convert to datetime (automatically converted to UTC based on timezone in units)
        data_nc_datetimes = num2date(data_nc_times, units=data_nc_timevar.units, only_use_cftime_datetimes=False, only_use_python_datetimes=True)
        #nptimes = data_nc_datetimes.astype('datetime64[ns]') #convert to numpy first, pandas does not take all cftime datasets
        
    
    
    if retrieve_ids is not False:
        data_nc_datetimes_pd = pd.Series(data_nc_datetimes,index=retrieve_ids).dt.round(freq='S')
    else:
        data_nc_datetimes_pd = pd.Series(data_nc_datetimes).dt.round(freq='S')
    
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






def get_hisstationlist(file_nc,varname):
    from netCDF4 import Dataset, chartostring
    import pandas as pd
    import numpy as np
    import warnings
    
    #from dfm_tools.get_nc_helpers import get_ncvarobject, get_ncvardimlist, get_variable_timevardim
    
    data_nc = Dataset(file_nc)
    data_nc_varname = get_ncvarobject(file_nc, varname) #check if var exists and get var_object
    varname_dims = data_nc_varname.dimensions
    varn_time, dimn_time = get_variable_timevardim(file_nc=file_nc, varname=varname)
    
    vars_pd, dims_pd = get_ncvardimlist(file_nc=file_nc)
    vars_pd_stats = vars_pd[(vars_pd['dtype']=='|S1') & (vars_pd['dimensions'].apply(lambda x: dimn_time not in x))]
    
    
    if varname in vars_pd_stats['nc_varkeys'].tolist(): 
        vars_pd_stats = vars_pd[vars_pd['nc_varkeys']==varname]
        
    #create lists of station variable names and dimensions, that correspond to dimensions of varname
    varname_stationdimname_list = []
    varname_stationvarname_list = []
    for iR, vars_pd_statrow in vars_pd_stats.iterrows():
        for iDV, varname_dim in enumerate(varname_dims):
            if varname_dim in vars_pd_statrow['dimensions']:
                varname_stationdimname_list.append(varname_dim)
                varname_stationvarname_list.append(vars_pd_statrow['nc_varkeys'])
    
    #create dataframe of station names coupled to varname
    if varname_stationdimname_list == []:
        raise Exception('ERROR: no dimension in %s variable that corresponds to station-like variables (or none present):\n%s'%(varname, vars_pd_stats['nc_varkeys']))
    else:
        var_station_names_pd = pd.DataFrame(columns = None)
        for iSV, varname_stationvarname in enumerate(varname_stationvarname_list):
            station_name = data_nc.variables[varname_stationvarname]
            if varname_stationdimname_list[iSV] in station_name.dimensions:
                station_name_char = station_name[:]
                try:
                    station_name_list_raw = chartostring(station_name_char)
                except: #for glossis netCDF file with probably invalidly stored station names
                    warnings.warn('station list could not be decoded with utf-8, now done with latin1 but the netCDF file might be corrupt and the station names sometimes unreadable')
                    """
                    station_name_list_raw_bytes = chartostring(station_name_char,encoding='bytes')
                    for iS,stat in enumerate(station_name_list_raw_bytes):
                        try:
                            stat.decode('utf-8')
                        except:
                            print('stat %d is not utf-8:\n\tbytes decoding: %s\n\tlatin-1 decoding: %s'%(iS, stat, stat.decode('latin-1')))
                    """
                    station_name_list_raw = chartostring(station_name_char,encoding='latin1')
                    
                station_name_list = np.char.strip(station_name_list_raw) #necessary step for Sobek and maybe others
                var_station_names_pd[varname_stationvarname] = station_name_list
                
        #get coordinates of stations (only works for stations, not for crs/gs since these variables have more than 1 dimension)
        #vars_pd_statlocs = vars_pd[(vars_pd['ndims']==1) & (vars_pd['dimensions'].astype(str).str.contains(varname_stationdimname_list[iSV]))] 
        vars_pd_statlocs = vars_pd[(vars_pd['ndims']==1) & (vars_pd['dimensions'].apply(lambda x: varname_stationdimname_list[iSV] in x))]
        
        coord_varnames = vars_pd_statlocs['nc_varkeys'].tolist()
        for iC, coord_varname in enumerate(coord_varnames):
            station_coordn = data_nc.variables[coord_varname]
            var_station_names_pd[coord_varname] = station_coordn[:]

    return var_station_names_pd




def get_stationid_fromstationlist(station_name_list_pd, station, varname):
    import numpy as np
    import pandas as pd
    
    station_pd = pd.Series(station)
    bool_nrows_req = station_pd.shape[0]
    bool_nrows_avai = station_name_list_pd.shape[0]
    bool_ncols = station_name_list_pd.shape[1]
    
    stations_bool_reqinfile_allcols = np.zeros((bool_nrows_req,bool_ncols),dtype=bool)
    stations_bool_fileinreq_allcols = np.zeros((bool_nrows_avai,bool_ncols),dtype=bool)
    for iCol in range(bool_ncols):
        #check if all requested stations are in netcdf file
        stations_bool_reqinfile_allcols[:,iCol] = station_pd.isin(station_name_list_pd.iloc[:,iCol])
        stations_bool_fileinreq_allcols[:,iCol] = station_name_list_pd.iloc[:,iCol].isin(station_pd)
    
    stations_bool_reqinfile = stations_bool_reqinfile_allcols.any(axis=1)
    if not stations_bool_reqinfile.all():
        raise Exception('ERROR: not all requested stations are in netcdf file:\n%s\navailable in netcdf file are:\n%s\nUse this command to obtain full list as variable:\nfrom dfm_tools.get_nc_helpers import get_hisstationlist; station_name_list_pd = get_hisstationlist(file_nc=file_nc,varname="%s")'%(station_pd[~stations_bool_reqinfile], station_name_list_pd, varname))
    #get ids of requested stations in netcdf file
    station_bool_fileinreq = stations_bool_fileinreq_allcols.any(axis=1)
    station_ids = list(np.where(station_bool_fileinreq)[0])

    return station_ids





