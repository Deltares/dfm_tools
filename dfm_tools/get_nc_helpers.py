# -*- coding: utf-8 -*-
"""
Created on Fri Feb 14 12:43:19 2020

@author: veenstra
"""

    
    

def get_ncfilelist(file_nc, multipart=None):
    #get list of mapfiles
    import re
    import glob
    import os
    
    if not os.path.exists(file_nc):
        raise Exception('ERROR: file does not exist: %s'%(file_nc))
    
    lastpart = file_nc.split('_')[-2]
    nctype = file_nc.split('_')[-1]
    if file_nc.endswith('_%s'%(nctype)) and multipart != False and len(lastpart) == 4 and lastpart.isdigit(): #if part before '_map.nc' is eg '0000'
        filename_start = re.compile('(.*)_([0-9]+)_%s'%(nctype)).search(file_nc).group(1)
        #filename_number = re.compile('(.*)_([0-9]+)_map.nc').search(file_nc).group(2)
        #file_ncs = [file_nc.replace('_%s_map.nc','_%04d_map.nc'%(filename_number, domain_id)) for domain_id in range(ndomains)]
        file_ncs = glob.glob('%s*_%s'%(filename_start,nctype))
    else:
        file_ncs = [file_nc]
    return file_ncs




def get_varname_mapnc(data_nc,varname_requested):
    import pandas as pd
    
    #VARIABLE names used within different versions of Delft3D-Flexible Mesh
    varnames_list = pd.DataFrame()
    varnames_list['time'] = ['time','nmesh2d_dlwq_time','',''] # time
    
    varnames_list['mesh2d_node_x'] = ['mesh2d_node_x','NetNode_x','mesh2d_agg_node_x',''] # x-coordinate of nodes
    varnames_list['mesh2d_node_y'] = ['mesh2d_node_y','NetNode_y','mesh2d_agg_node_y',''] # y-coordinate of nodes
    varnames_list['mesh2d_node_z'] = ['mesh2d_node_z','NetNode_z','',''] # z-coordinate of nodes
    
    varnames_list['mesh2d_face_x'] = ['mesh2d_face_x','FlowElem_xzw','mesh2d_agg_face_x',''] # x-coordinate of faces
    varnames_list['mesh2d_face_y'] = ['mesh2d_face_y','FlowElem_yzw','mesh2d_agg_face_y',''] # y-coordinate of faces
    
    varnames_list['mesh2d_edge_x'] = ['mesh2d_edge_x','','',''] # x-coordinate of velocity-points
    varnames_list['mesh2d_edge_y'] = ['mesh2d_edge_y','','',''] # y-coordinate of velocity-points
    
    varnames_list['mesh2d_edge_nodes'] = ['mesh2d_edge_nodes','NetLink','',''] # 'link between two netnodes' / 'Mapping from every edge to the two nodes that it connects'
    varnames_list['mesh2d_face_nodes'] = ['mesh2d_face_nodes','NetElemNode','mesh2d_agg_face_nodes',''] # 
    
    varnames_list['mesh2d_face_x_bnd'] = ['mesh2d_face_x_bnd','FlowElemContour_x','mesh2d_agg_face_x_bnd',''] # x-coordinates of flow element contours
    varnames_list['mesh2d_face_y_bnd'] = ['mesh2d_face_y_bnd','FlowElemContour_y','mesh2d_agg_face_y_bnd',''] # y-coordinates of flow element contours
    
    varnames_list['mesh2d_flowelem_domain'] = ['mesh2d_flowelem_domain','FlowElemDomain','',''] # flow element domain
    varnames_list['mesh2d_flowelem_bl'] = ['mesh2d_flowelem_bl','FlowElem_bl','',''] # bed level
    varnames_list['mesh2d_flowelem_ba'] = ['mesh2d_flowelem_ba','FlowElem_bac','',''] # area (m2) of cell faces
    
    varnames_list['mesh2d_layer_z'] = ['mesh2d_layer_z','LayCoord_cc','',''] # 
    
    #non-grid variables necessary for layer calculation for intersection/cross section) funtion
    varnames_list['mesh2d_s1'] = ['mesh2d_s1','','',''] # water level
    varnames_list['mesh2d_flowelem_bl'] = ['mesh2d_flowelem_bl','','',''] # bed level

    #varnames_list['mesh2d_ucx'] = ['mesh2d_ucx','ucx','',''] # 
    #varnames_list['mesh2d_ucy'] = ['mesh2d_ucy','ucy','',''] # 
    #varnames_list['mesh2d_sa1'] = ['mesh2d_sa1','sa1','',''] # 
    #varnames_list['mesh2d_tem1'] = ['mesh2d_tem1','tem1','',''] # 
    
    
    ### DIMENSION names used within different versions of Delft3D-Flexible Mesh
    #dimnames_list = pd.DataFrame()
    varnames_list['nmesh2d_node'] = ['nmesh2d_node','mesh2d_nNodes','nNetNode',''] # number of nodes
    varnames_list['nmesh2d_face'] = ['nmesh2d_face','mesh2d_nFaces','nNetElem','nFlowElem'] # number of faces
    varnames_list['nmesh2d_edge'] = ['nmesh2d_edge','nNetLink','',''] # number of velocity-points
    
    varnames_list['nmesh2d_layer'] = ['nmesh2d_layer','mesh2d_nLayers','laydim','nmesh2d_layer_dlwq'] # layer
    
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
    #if varname is None:
    #    print('WARNING: var/dim name %s or equivalent not found in netCDF file with variables:\n%s \nand dimensions:\n%s'%(varname_requested, data_nc_varnames_list, data_nc_dimnames_list))
    
    return varname



def get_ncvardims(file_nc, varname):
    from netCDF4 import Dataset
    
    data_nc = Dataset(file_nc)
    # check if requested variable is in netcdf
    nc_varkeys = list(data_nc.variables.keys())
    nc_dimkeys = list(data_nc.dimensions.keys())
    nc_varlongnames = []
    for nc_var in data_nc.variables:
        try:
            nc_varlongnames.append(data_nc.variables[nc_var].long_name)
        except:
            nc_varlongnames.append('NO long_name defined')
    if varname not in nc_varkeys:
        raise Exception('ERROR: requested variable %s not in netcdf, available are:\n%s'%(varname, '\n'.join(map(str,['%-25s: %s'%(nck,ncln) for nck,ncln in zip(nc_varkeys, nc_varlongnames)]))))
    
    nc_values = data_nc.variables[varname]
    nc_values_shape = nc_values.shape
    nc_values_dims = nc_values.dimensions
    #nc_values_ndims = len(nc_values_dims)
    return nc_varkeys, nc_dimkeys, nc_values, nc_values_shape, nc_values_dims



def ghostcell_filter(file_nc):
    import numpy as np
    from netCDF4 import Dataset
    
    #from dfm_tools.get_nc_helpers import get_varname_mapnc
    
    data_nc = Dataset(file_nc)
    
    varn_domain = get_varname_mapnc(data_nc,'mesh2d_flowelem_domain')
    if varn_domain is not None: # domain variable is present, so there are multiple domains
        ghostcells_bool = True
        domain = data_nc.variables[varn_domain][:]
        domain_no = np.bincount(domain).argmax() #meest voorkomende domeinnummer
        nonghost_ids = domain==domain_no
    else:
        ghostcells_bool = False
        nonghost_ids = None
    return ghostcells_bool, nonghost_ids



def get_timesfromnc(file_nc):
    from netCDF4 import Dataset,num2date#,date2num
    import numpy as np
    import pandas as pd
    
    #from dfm_tools.get_nc_helpers import get_varname_mapnc
    
    data_nc = Dataset(file_nc)
    varname_time = get_varname_mapnc(data_nc,'time')
    data_nc_timevar = data_nc.variables[varname_time]
    
    if len(data_nc_timevar)<3: #this rarely is the case, but just to be sure
        data_nc_times = data_nc_timevar[:]
    else:
        time0 = data_nc_timevar[0] 
        time1 = data_nc_timevar[1] 
        time2 = data_nc_timevar[2]
        timeend = data_nc_timevar[-1]
        timeinc = time2-time1 # the interval between 0 and 1 is not per definition representative, so take 1 and 2
        
        data_nc_times_from1 = np.arange(time1,timeend+timeinc,timeinc)
        data_nc_times = np.concatenate([[time0],data_nc_times_from1])
    data_nc_datetimes = num2date(data_nc_times, units = data_nc_timevar.units)
    data_nc_datetimes_pd = pd.Series(data_nc_datetimes).dt.round(freq='S')
    
    return data_nc_datetimes_pd




def get_timeid_fromdatetime(data_nc_datetimes_pd, timestep):
    import numpy as np
    import pandas as pd
    
    timestep_pd = pd.Series(timestep)#.dt.round(freq='S')

    #check if all requested times (timestep) are in netcdf file
    times_bool_reqinfile = timestep_pd.isin(data_nc_datetimes_pd)
    if not (times_bool_reqinfile == True).all():
        raise Exception('ERROR: not all requested times are in netcdf file:\n%s\navailable in netcdf file are:\n\tstart: %s\n\tstop: %s\n\tinterval: %s'%(timestep_pd[-times_bool_reqinfile], data_nc_datetimes_pd.iloc[0], data_nc_datetimes_pd.iloc[-1], data_nc_datetimes_pd.iloc[1]-data_nc_datetimes_pd.iloc[0]))
        
    #get ids of requested times in netcdf file
    times_bool_fileinreq = data_nc_datetimes_pd.isin(timestep_pd)
    time_ids = np.where(times_bool_fileinreq)[0]
    
    return time_ids



def get_hisstationlist(file_nc):
    from netCDF4 import Dataset, chartostring
    import pandas as pd
    
    data_nc = Dataset(file_nc)
    #varn_station_name = get_varname_mapnc(data_nc,'station_name')
    station_name_char = data_nc.variables['station_name'][:]
    station_name_list = chartostring(station_name_char)
    
    station_name_list_pd = pd.Series(station_name_list)
    
    return station_name_list_pd




def get_stationid_fromstationlist(station_name_list_pd, station):
    import numpy as np
    import pandas as pd
    
    station_pd = pd.Series(station)

    #check if all requested stations are in netcdf file
    stations_bool_reqinfile = station_pd.isin(station_name_list_pd)
    if not (stations_bool_reqinfile == True).all():
        raise Exception('ERROR: not all requested stations are in netcdf file:\n%s\navailable in netcdf file are:\n%s'%(station_pd[-stations_bool_reqinfile], station_name_list_pd))
    
    #get ids of requested stations in netcdf file
    station_bool_fileinreq = station_name_list_pd.isin(station_pd)
    station_ids = np.where(station_bool_fileinreq)[0]
    #station_names = np.where(station_bool_fileinreq)[0]

    return station_ids





