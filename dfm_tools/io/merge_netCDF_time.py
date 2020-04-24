# -*- coding: utf-8 -*-
"""
Created on Fri Apr 24 09:50:42 2020

@author: veenstra
"""

def merge_netCDF_time(tstart, tstop, dir_data, list_subfolders=None, nc_prefix=''):
    
    from netCDF4 import Dataset
    import glob
    import os
    import re
    import datetime as dt
    
    #SETTINGS
    tstart = dt.datetime(2016,4,28)
    tstop = dt.datetime(2016,5,3)
    dir_data = r'n:\Deltabox\Postbox\Veenstra, Jelmer\vanReimer\GLBu0.08_expt_91.2'
    list_subfolders = ''
    nc_prefix = 'HYCOM_ST_GoO_'
    file_to = os.path.join(dir_data, 'merged', '%s%sto%s.nc'%(nc_prefix, merge_tstart.strftime('%Y%m%d'), merge_tstop.strftime('%Y%m%d'))
    file_src_listall = glob.glob(os.path.join(dir_data, '%s*.nc'%(nc_prefix)))
    fn_dateformat = '%Y%m%d'
    fn_match_pattern_dates = '%s(.*).nc'%(nc_prefix)
    ###############
    
    p = re.compile(fn_match_pattern_dates)
    try:
        data_to.close()
    except:
        pass
    merge_ndays = ((merge_tstop-merge_tstart).total_seconds()/3600/24)+1

    file_src_list = []
    for iF, file_src in enumerate(file_src_listall):
        fn_result = p.search(file_src)
        fn_timestamp = fn_result.group(1)
        fn_timestamp_dt = dt.datetime.strptime(fn_timestamp,fn_dateformat)
        if (fn_timestamp_dt >= merge_tstart) & (fn_timestamp_dt <= merge_tstop):
            file_src_list.append(file_src)
    if file_src_list == []:
        raise Exception('no files found between %s and %s'%(merge_tstart, merge_tstop))
    if len(file_src_list) != merge_ndays:
        raise Exception('number of found files (%d) is not equal to number of requested days (%d)'%(len(file_src_list), merge_ndays))
        
    for iF, file_src in enumerate(file_src_list):
        data_src = Dataset(file_src)
        if iF == 0: #initiate empty file
            data_to = Dataset(file_to, 'w', format="NETCDF3_CLASSIC")
            
            #copy nc attributes other than dimensions and variables 
            data_src_attrlist = data_src.ncattrs()
            for ncattrname in data_src_attrlist:
                if ncattrname not in ['variables','dimensions']:
                    data_to.setncattr(ncattrname, data_src.getncattr(ncattrname))
            
            #copy dimensions (make time unlimited)
            data_src_dimlist = list(data_src.dimensions.keys())
            for dimname in data_src_dimlist:
                data_src_dimsize = data_src.dimensions[dimname].size
                if dimname == 'time':
                    data_to.createDimension(dimname, None)
                else:
                    data_to.createDimension(dimname, data_src_dimsize)
            
            #copy variables (fill three already)
            data_src_varlist = list(data_src.variables.keys())
            for varname in data_src_varlist:
                data_src_var=data_src.variables[varname]
                data_src_dims=data_src_var.dimensions
                data_src_type=data_src_var.dtype.str
                if '_FillValue' in data_src_var.ncattrs():
                    data_to.createVariable(varname,data_src_type,data_src_dims, fill_value=data_src_var.getncattr('_FillValue'))
                else:
                    data_to.createVariable(varname,data_src_type,data_src_dims)
                for attr_name in data_src_var.ncattrs():
                    if attr_name != '_FillValue': #to avoid error with netCDF4 library version 1.5.3, fill values are not convenient for dflowfm
                        #print(attr_name,'=',data_src_var.getncattr(attr_name))
                        data_to.variables[varname].setncattr(attr_name,data_src_var.getncattr(attr_name))
                if varname in ['depth', 'lat', 'lon']:
                    data_to.variables[varname][:] = data_src.variables[varname][:]
    
        for attr_name in data_src.variables[varname].ncattrs():
            if data_src.variables[varname].getncattr(attr_name) != data_to.variables[varname].getncattr(attr_name):
                raise Exception('attribute %s not equal for variable %s'%(attr_name,varname))
    
    
        data_src_ntimes = data_src.variables['time'].shape[0]
        data_to_ntimes = data_to.variables['time'].shape[0]
        
        #append current data to netcdf files
        data_to.variables['time'][data_to_ntimes:data_to_ntimes+data_src_ntimes] = data_src.variables['time'][:]
        data_to.variables['salinity'][data_to_ntimes:data_to_ntimes+data_src_ntimes,...] = data_src.variables['salinity'][:]
        data_to.variables['water_temp'][data_to_ntimes:data_to_ntimes+data_src_ntimes,...] = data_src.variables['water_temp'][:]
        data_src.close()
    data_to.close()
        
        
    """
    from dfm_tools.get_nc import get_netdata, get_ncmodeldata, plot_netmapdata
    from dfm_tools.get_nc_helpers import get_ncvardimlist#, get_ncfilelist
    
    file_nc = file_to
    vars_pd, dims_pd = get_ncvardimlist(file_nc=file_nc)
    
    data_fromnc = get_ncmodeldata(file_nc=file_nc, varname='salinity', timestep='all', layer='all')
    
    
    """