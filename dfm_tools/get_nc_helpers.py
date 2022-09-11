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

Created on Fri Feb 14 12:43:19 2020

@author: veenstra

helper functions for functions in get_nc.py
"""


import xarray as xr
import pandas as pd
import re
import glob
import os


#TODO: remove this, easier with xarray
"""
def ncdump_OLD(file_nc):
    from netCDF4 import Dataset
    #import pandas as pd
    #import numpy as np
    
    data_nc = Dataset(file_nc)
    
    # NetCDF global attributes
    nc_attrs = data_nc.ncattrs()
    print("NetCDF Global Attributes:")
    for nc_attr in nc_attrs:
        print('\t%s: %s'%(nc_attr, data_nc.getncattr(nc_attr)))
    
    # Dimension shape information.
    nc_dims = list(data_nc.dimensions.keys())  # list of nc dimensions
    print("NetCDF dimension information:")
    for dim in nc_dims:
        if 'unlimited' in str(data_nc.dimensions[dim]):
            print("\t%s = UNLIMITED (currently %i)"%(dim, data_nc.dimensions[dim].size))
        else:    
            print("\t%s = %i"%(dim, data_nc.dimensions[dim].size))
    
    # Variable information.
    nc_vars = list(data_nc.variables.keys())  # list of nc variables
    print("NetCDF variable information:")
    for var in nc_vars:
        #if var not in nc_dims:
        print('\t%s %s %s'%(data_nc.variables[var].dtype, var, str(data_nc.variables[var].dimensions)))
        print("\t\tshape: %s"%(str(data_nc.variables[var].shape)))
        for ncattr in data_nc.variables[var].ncattrs():
            print('\t\t%s: %s'%(ncattr,data_nc.variables[var].getncattr(ncattr)))
    data_nc.close()
    #return nc_attrs, nc_dims, nc_vars
"""



def get_ncfilelist(file_nc, multipart=None):
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
    #if varname is None:
    #    varname = get_vardimname(data_nc_dimnames_list)
    
    return varname




def get_ncvardimlist(file_nc):
    
    data_xr = xr.open_dataset(file_nc)
    nc_varkeys = data_xr.variables.mapping.keys()
    nc_dimkeys = data_xr.dims.mapping.keys()
    
    emptycol = [['']]*len(nc_varkeys)
    emptycol_str = ['']*len(nc_varkeys)
    vars_pd = pd.DataFrame({'nc_varkeys':nc_varkeys, 'shape':emptycol, 'dimensions':emptycol, 'dtype':emptycol,
                            'standard_name':emptycol_str,'long_name':emptycol_str,'coordinates':emptycol_str,'units':emptycol_str,'mesh':emptycol_str,'location':emptycol_str},
                           index=nc_varkeys)
    dims_pd = pd.DataFrame({'nc_dimkeys': nc_dimkeys, 'size': [['']]*len(nc_dimkeys)},
                           index=nc_dimkeys)
    var_attr_name_list = ['standard_name','long_name','coordinates','units','mesh','location']
    for varkey in nc_varkeys:
        #get non-attribute properties of netcdf variable
        vars_pd.loc[varkey,'nc_varkeys'] = varkey #TODO: this one is not necessary anymore, use index instead
        vars_pd.loc[varkey,'shape'] = data_xr.variables[varkey].shape
        vars_pd.loc[varkey,'dimensions'] = data_xr.variables[varkey].dims
        vars_pd.loc[varkey,'dtype'] = data_xr.variables[varkey].dtype
        #get attributes properties of netcdf variable
        for attr_name in var_attr_name_list:
            try:
                vars_pd.loc[varkey, attr_name] = data_xr.variables[varkey].attrs[attr_name]
            except:
                pass
    vars_pd['ndims'] = vars_pd['dimensions'].apply(len)
    for dimkey in nc_dimkeys:
        #get non-attribute properties of netcdf variable
        dims_pd.loc[dimkey,'nc_dimkeys'] = dimkey
        dims_pd.loc[dimkey,'size'] = data_xr.dims[dimkey]
    data_xr.close()
    return vars_pd, dims_pd


""" #TODO: remove this def
def get_ncvardimlist_OLD(file_nc):
    from netCDF4 import Dataset
    import pandas as pd
    #import numpy as np
    
    data_nc = Dataset(file_nc)
    
    nc_varkeys = list(data_nc.variables.keys())
    nc_dimkeys = list(data_nc.dimensions.keys())
    #vars_pd = pd.DataFrame({'nc_varkeys': nc_varkeys, 'shape': [['']]*len(nc_varkeys), 'dimensions': [['']]*len(nc_varkeys), 'dtype': [['']]*len(nc_varkeys)})
    emptycol = [['']]*len(nc_varkeys)
    emptycol_str = ['']*len(nc_varkeys)
    vars_pd = pd.DataFrame({'nc_varkeys':nc_varkeys, 'shape':emptycol, 'dimensions':emptycol, 'dtype':emptycol,
                            'standard_name':emptycol_str,'long_name':emptycol_str,'coordinates':emptycol_str,'units':emptycol_str,'mesh':emptycol_str,'location':emptycol_str})
    dims_pd = pd.DataFrame({'nc_dimkeys': nc_dimkeys, 'name': [['']]*len(nc_dimkeys), 'size': [['']]*len(nc_dimkeys)})
    var_attr_name_list = ['standard_name','long_name','coordinates','units','mesh','location']
    for iV, nc_var in enumerate(data_nc.variables):
        #get non-attribute properties of netcdf variable
        vars_pd.loc[iV,'shape'] = data_nc.variables[nc_var].shape
        vars_pd.loc[iV,'dimensions'] = data_nc.variables[nc_var].dimensions
        vars_pd.loc[iV,'dtype'] = data_nc.variables[nc_var].dtype
        #get attributes properties of netcdf variable
        for attr_name in var_attr_name_list:
            #if iV==0:
            #    vars_pd[attr_name] = ''
            try:
                vars_pd.loc[iV, attr_name] = data_nc.variables[nc_var].getncattr(attr_name)
            except:
                pass
    vars_pd['ndims'] = [len(x) for x in vars_pd['dimensions']]
    for iD, nc_dim in enumerate(data_nc.dimensions):
        #get non-attribute properties of netcdf variable
        dims_pd.loc[iD,'name'] = data_nc.dimensions[nc_dim].name
        dims_pd.loc[iD,'size'] = data_nc.dimensions[nc_dim].size
    data_nc.close()
    return vars_pd, dims_pd
"""


def get_varnamefrom_keyslongstandardname(file_nc, varname):
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
        varname = vars_pd.index[varid]
        print('varname found in long_name attribute')
    elif varname in nc_varstandardnames:
        varid = nc_varstandardnames.index(varname)
        varname = vars_pd.index[varid]
        print('varname found in standard_name attribute')
    else:
        raise Exception('ERROR: requested variable %s not in netcdf, available are:\n%s\nUse this command to obtain full list as variable:\nfrom dfm_tools.get_nc_helpers import get_ncvardimlist\nvars_pd, dims_pd = get_ncvardimlist(file_nc=file_nc)\nnote that you can retrieve variables by keys, standard_name or long_name attributes'%(varname, vars_pd))
    
    return varname



def ghostcell_filter(file_nc):
    import numpy as np
    from netCDF4 import Dataset
    
    #from dfm_tools.get_nc_helpers import get_varname_fromnc
    
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





def get_variable_timevardim(file_nc, varname):
    #get corresponding time variable name
    from netCDF4 import Dataset
    
    #from dfm_tools.get_nc_helpers import get_ncvarobject
    
    data_nc = Dataset(file_nc)
    varname = get_varnamefrom_keyslongstandardname(file_nc, varname) #get varname from varkeys/standardname/longname if exists
    nc_varobject = data_nc.variables[varname]
    
    varn_time = None
    dimn_time = None
    varlist_wunits = data_nc.get_variables_by_attributes(units=lambda v: v is not None)
    for var_lookup in varlist_wunits:
        if 'since' in var_lookup.units and var_lookup.dimensions[0] in nc_varobject.dimensions:
            dimn_time = var_lookup.dimensions[0]
            varn_time = var_lookup.name
            break
    
    data_nc.close()
    return varn_time, dimn_time



    
def get_timesfromnc(file_nc, varname='time', retrieve_ids=False, keeptimezone=True, silent=False):
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
        if not silent:
            print('reading time dimension: only requested indices')
        listtype_range = [list, range, np.ndarray]
        if type(retrieve_ids) not in listtype_range:
            raise Exception('ERROR: argument retrieve_ids should be a list')
        #convert to positive index, make unique(+sort), convert to list because of indexing with np.array of len 1 errors sometimes
        retrieve_ids = list(np.unique(np.array(range(time_length))[retrieve_ids]))
        data_nc_times = data_nc_timevar[retrieve_ids]
    elif len(data_nc_timevar)<3: #check if time dimension is shorter than 3 items
        data_nc_times = data_nc_timevar[:]
        if not silent:
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
            if not silent:
                print('reading time dimension: reconstruct array')
            data_nc_times_from1 = np.arange(time1,timemin1,timeinc_poststart)
            data_nc_times = np.concatenate([[time0],data_nc_times_from1,[timemin1]])
            if data_nc_timevar.shape[0] != len(data_nc_times):#test if len of reconstructed timeseries is same as len of timevar in netCDF, retrieve entire array
                if not silent:
                    print('reading time dimension: reconstruction failed, read entire array')
                data_nc_times = data_nc_timevar[:]
        else:
            if not silent:
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
            if not silent:
                print('retrieving original timezone succeeded, no conversion to UTC/GMT applied')
        except:
            if not silent:
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
    
    data_nc.close()
    return data_nc_datetimes_pd



#TODO: remove after moving to xarray for time selection
def get_timeid_fromdatetime(data_nc_datetimes_pd, timestep):
    import numpy as np
    import pandas as pd
    
    timestep_pd = pd.Series(timestep)#.dt.round(freq='S')

    #check if all requested times (timestep) are in netcdf file
    times_bool_reqinfile = timestep_pd.isin(data_nc_datetimes_pd)
    if not (times_bool_reqinfile == True).all():
        raise Exception('ERROR: not all requested times are in netcdf file:\n%s\navailable in netcdf file are:\n%s\nUse this command to obtain full list as variable:\nfrom dfm_tools.get_nc_helpers import get_timesfromnc\ndata_nc_datetimes_pd = get_timesfromnc(file_nc=file_nc)'%(timestep_pd[-times_bool_reqinfile], data_nc_datetimes_pd))
        
    #get ids of requested times in netcdf file
    times_bool_fileinreq = data_nc_datetimes_pd.isin(timestep_pd)
    time_ids = np.where(times_bool_fileinreq)[0]
    
    return time_ids




def get_hisstationlist(file_nc, varname='waterlevel'):
    #encoding = {'station_lon': {'_FillValue': None}, #TODO lon/lat are now fillvalue if nan
    #            }
    data_xr = xr.open_dataset(file_nc)
    varname = get_varnamefrom_keyslongstandardname(file_nc, varname) #get varname from varkeys/standardname/longname if exists
    vardims = data_xr[varname].dims
    
    vars_pd, dims_pd = get_ncvardimlist(file_nc=file_nc)
    bool_vars_dtypestr = vars_pd['dtype'].astype(str).str.startswith('|S') | (vars_pd['dtype']=='object')#& (vars_pd['ndims']==1) #TODO: better check for dtype string?
    vars_pd_stats = vars_pd.loc[bool_vars_dtypestr]
    
    if len(vars_pd_stats) == 0:
        raise Exception('ERROR: no dimension in %s variable that corresponds to station-like variables (or none present):\n%s'%(varname, vars_pd_stats['nc_varkeys']))
        
    dimkey_use = None
    for dimtuple in vars_pd_stats['dimensions']:
        dimkey = dimtuple[0]
        if dimkey in vardims:
            dimkey_use = dimkey
    
    statlist_pd = data_xr[dimkey_use].to_dataframe()
    
    for colname in statlist_pd.columns:
        if isinstance(statlist_pd[colname][0],bytes): #TODO: better check would be data_xr.variables[colname].dtype=='S256'
            print(f'variable {colname}: converting bytes to str')
            statlist_pd[colname] = pd.Series(data_xr[colname].astype(str)).str.strip() #replace bytes by stripped strings
    
    data_xr.close()
    return statlist_pd


""" #TODO: remove this def
def get_hisstationlist_OLD(file_nc,varname):
    from netCDF4 import Dataset, chartostring
    import pandas as pd
    import numpy as np
    import warnings
    
    #from dfm_tools.get_nc_helpers import get_ncvarobject, get_ncvardimlist, get_variable_timevardim
    
    data_nc = Dataset(file_nc)
    varname = get_varnamefrom_keyslongstandardname(file_nc, varname) #get varname from varkeys/standardname/longname if exists
    nc_varobject = data_nc.variables[varname]
    varname_dims = nc_varobject.dimensions
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
        raise Exception('ERROR OLD: no dimension in %s variable that corresponds to station-like variables (or none present):\n%s'%(varname, vars_pd_stats['nc_varkeys']))
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
                    #station_name_list_raw_bytes = chartostring(station_name_char,encoding='bytes')
                    #for iS,stat in enumerate(station_name_list_raw_bytes):
                    #    try:
                    #        stat.decode('utf-8')
                    #    except:
                    #        print('stat %d is not utf-8:\n\tbytes decoding: %s\n\tlatin-1 decoding: %s'%(iS, stat, stat.decode('latin-1')))
                    station_name_list_raw = chartostring(station_name_char,encoding='latin-1')
                station_name_list = np.char.strip(station_name_list_raw) #necessary step for Sobek and maybe others
                var_station_names_pd[varname_stationvarname] = station_name_list
                
        #get coordinates of stations (only works for stations, not for crs/gs since these variables have more than 1 dimension)
        #vars_pd_statlocs = vars_pd[(vars_pd['ndims']==1) & (vars_pd['dimensions'].astype(str).str.contains(varname_stationdimname_list[iSV]))] 
        vars_pd_statlocs = vars_pd[(vars_pd['ndims']==1) & (vars_pd['dimensions'].apply(lambda x: varname_stationdimname_list[iSV] in x))]
        
        coord_varnames = vars_pd_statlocs['nc_varkeys'].tolist()
        for iC, coord_varname in enumerate(coord_varnames):
            station_coordn = data_nc.variables[coord_varname]
            var_station_names_pd[coord_varname] = station_coordn[:]
    
    data_nc.close()
    return var_station_names_pd
"""


#TODO: might not be necessary to have when using xarray
#TODO: can this be simplyfied? maybe possible to replace bytes by stripped string in data_xr immediately, instead of in every subsequent step
def get_stationid_fromstationlist(stations_pd, stationlist):
    import numpy as np
    import pandas as pd
    
    if not isinstance(stationlist,list):
        stationlist = [stationlist]
    station_pd = pd.Series(stationlist)
    bool_nrows_req = station_pd.shape[0]
    bool_nrows_avai = stations_pd.shape[0]
    bool_ncols = stations_pd.shape[1]
    
    stations_bool_reqinfile_allcols = np.zeros((bool_nrows_req,bool_ncols),dtype=bool)
    stations_bool_fileinreq_allcols = np.zeros((bool_nrows_avai,bool_ncols),dtype=bool)
    for iCol in range(bool_ncols):
        #check if all requested stations are in netcdf file
        stations_bool_reqinfile_allcols[:,iCol] = station_pd.isin(stations_pd.iloc[:,iCol])
        stations_bool_fileinreq_allcols[:,iCol] = stations_pd.iloc[:,iCol].isin(station_pd)
    
    stations_bool_reqinfile = stations_bool_reqinfile_allcols.any(axis=1)
    if not stations_bool_reqinfile.all():
        raise Exception(f'ERROR: not all requested stations are in netcdf file:\n{station_pd[~stations_bool_reqinfile]}\navailable in netcdf file are:\n{stations_pd}')
    #get ids of requested stations in netcdf file
    station_bool_fileinreq = stations_bool_fileinreq_allcols.any(axis=1)
    station_ids = list(np.where(station_bool_fileinreq)[0])

    return station_ids





