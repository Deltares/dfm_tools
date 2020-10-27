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

Created on Fri Apr 24 09:50:42 2020

@author: veenstra
"""

def merge_netCDF_time(tstart, tstop, tstep_sec, dir_data, nc_prefix, fn_match_pattern, fn_dateformat, subfolders=[''], dir_out=None, 
                      varn_time='time', dimn_time='time', renamevars=None, dfmtoolsattr=True):
    """
    this script works well for daily HYCOM data but should be made more generic for other types of (meteo) data.
    """
    
    from netCDF4 import Dataset
    import glob
    import os
    import re
    import datetime as dt
        
    if dir_out is None:
        dir_out = dir_data
    if not os.path.exists(dir_out):
        os.makedirs(dir_out)
        
    file_src_listall = glob.glob(os.path.join(dir_data, '%s*.nc'%(nc_prefix)))
    file_to = os.path.join(dir_out, '%s%sto%s.nc'%(nc_prefix, tstart.strftime(fn_dateformat), tstop.strftime(fn_dateformat)))
    p = re.compile(fn_match_pattern)
    
    try:
        data_to.close()
    except:
        pass
    merge_ntimesteps = ((tstop-tstart).total_seconds()/tstep_sec)+1

    file_src_list = []
    for iF, file_src in enumerate(file_src_listall):
        fn_result = p.search(file_src)
        fn_timestamp = fn_result.group(1)
        fn_timestamp_dt = dt.datetime.strptime(fn_timestamp,fn_dateformat)
        if (fn_timestamp_dt >= tstart) & (fn_timestamp_dt <= tstop):
            file_src_list.append(file_src)
    if file_src_list == []:
        raise Exception('no files found between %s and %s'%(tstart, tstop))
    if len(file_src_list) != merge_ntimesteps:
        raise Exception('number of found files (%d) is not equal to number of requested days (%d)'%(len(file_src_list), merge_ntimesteps))
        
    for iF, file_src in enumerate(file_src_list):
        print('processing: %s'%(file_src))
        data_src = Dataset(file_src)
        if iF == 0: #initiate empty file
            ncformat = data_src.file_format #eg "NETCDF3_CLASSIC"
            print('creating file in format: %s'%(ncformat))
            data_to = Dataset(file_to, 'w', format=ncformat)
            
            #copy nc attributes other than dimensions and variables 
            data_src_attrlist = data_src.ncattrs()
            for ncattrname in data_src_attrlist:
                if ncattrname not in ['variables','dimensions']:
                    data_to.setncattr(ncattrname, data_src.getncattr(ncattrname))
            if dfmtoolsattr:
                data_to.setncattr('comment', 'merged with dfm_tools.io.netCDF_utils.merge_netCDF_time() from https://github.com/openearth/dfm_tools')
            #copy dimensions (make dimn_time unlimited)
            data_src_dimlist = list(data_src.dimensions.keys())
            for dimname in data_src_dimlist:
                data_src_dimsize = data_src.dimensions[dimname].size
                if dimname == dimn_time: #this one is unlimited, necessary for merging
                    data_to.createDimension(dimname, None)
                else:
                    data_to.createDimension(dimname, data_src_dimsize)
            
            #copy variables (fill three already)
            data_src_varlist = list(data_src.variables.keys())
            timedep_varlist = []
            for varname in data_src_varlist:
                data_src_var=data_src.variables[varname]
                data_src_dims=data_src_var.dimensions
                data_src_type=data_src_var.dtype.str
                if dimn_time in data_src_dims:
                    data_src_var_timedep = True
                    timedep_varlist.append(varname)
                else:
                    data_src_var_timedep = False
                if '_FillValue' in data_src_var.ncattrs():
                    data_to.createVariable(varname,data_src_type,data_src_dims, fill_value=data_src_var.getncattr('_FillValue'))
                else:
                    data_to.createVariable(varname,data_src_type,data_src_dims)
                for attr_name in data_src_var.ncattrs():
                    if attr_name != '_FillValue': #_FillValue should be and is already set upon creating variable
                        #print(attr_name,'=',data_src_var.getncattr(attr_name))
                        data_to.variables[varname].setncattr(attr_name,data_src_var.getncattr(attr_name))
                if not data_src_var_timedep:
                    data_to.variables[varname][:] = data_src.variables[varname][:]
        
        #check if all attributes in both files are equal (mainly relevant for scale_factor and add_offset)
        for attr_name in data_src.variables[varname].ncattrs():
            if data_src.variables[varname].getncattr(attr_name) != data_to.variables[varname].getncattr(attr_name):
                raise Exception('attribute %s not equal for variable %s'%(attr_name,varname))
    
    
        data_src_ntimes = data_src.variables[varn_time].shape[0]
        data_to_ntimes = data_to.variables[varn_time].shape[0]
        
        #append current data to netcdf files
        for varname in timedep_varlist:
            data_to.variables[varname][data_to_ntimes:data_to_ntimes+data_src_ntimes,...] = data_src.variables[varname][:]
        data_src.close()
    
    if renamevars is not None:
        for varn_old in list(renamevars.keys()):
            varn_new = renamevars[varn_old]
            data_to.renameVariable(varn_old,varn_new)
    data_to.close()
    return file_to



