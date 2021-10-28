# -*- coding: utf-8 -*-
"""
Created on Wed Oct 27 22:07:34 2021

@author: veenstra

this test tests merging a (eg meteo) netcdf over time

"""

import os
import datetime as dt
import matplotlib.pyplot as plt
plt.close('all')

from netCDF4 import Dataset

from dfm_tools.get_nc_helpers import get_ncvardimlist#, get_ncfilelist
from dfm_tools.io.netCDF_utils import merge_netCDF_time

dir_testinput = r'c:\DATA\dfm_tools_testdata'
dir_output = '.'


"""
if mode == 'meteo':
    dir_data = 'P:\\1204257-dcsmzuno\\2014\\data\\meteo'
    subfolders = ['HIRLAM72_2011','HIRLAM72_2012','HIRLAM72_2013','HIRLAM72_2014','HIRLAM72_2015','HIRLAM72_2016','HIRLAM72_2017','HIRLAM72_2018','HIRLAM72_2019']
    ignorelist = ['P:\\1204257-dcsmzuno\\2014\\data\\meteo\\HIRLAM72_2016\\h72_2016_merged.nc']
    vars_orig = ['air_pressure_fixed_height','northward_wind','eastward_wind'] #voor meteo folder
    vars_dest = vars_orig
    convert_vars = [0,0,0]
elif mode == 'meteo_heatflux':
    dir_data = 'P:\\1204257-dcsmzuno\\2014\\data\\meteo-heatflux'
    subfolders = ['HIRLAM72_2011','HIRLAM72_2012','HIRLAM72_2013','HIRLAM72_2014','HIRLAM72_2015','HIRLAM72_2016','HIRLAM72_2017','HIRLAM72_2018','HIRLAM72_2019']
    ignorelist = []
    vars_orig = ['dew_point_temperature','air_temperature','cloud_area_fraction'] #voor meteo-heatflux folder
    vars_dest = vars_orig
    convert_vars = [1,1,2] #0 is no conversion, 1 is K to C (incl unit), 2 is cloud cover fraction (0-1) to cloud cover percentage (0-100) (unit is al %)
else:
    raise Exception('ERROR: wrong mode %s'%(mode))
"""

#SETTINGS
tstart = dt.datetime(2016,4,28)
tstop = dt.datetime(2016,5,3)
tstep_sec = 3600*24
dir_data = os.path.join(dir_testinput,'GLBu0.08_expt_91.2')
nc_prefix = 'HYCOM_ST_GoO_'
fn_match_pattern = '%s(.*).nc'%(nc_prefix)
fn_dateformat = '%Y%m%d'
#subfolders = ''
dir_out = dir_output #os.path.join(dir_data,'merged_new')
renamevars = {'salinity':'so', 'water_temp':'thetao'}
###############
file_merged = merge_netCDF_time(tstart=tstart, tstop=tstop, tstep_sec=tstep_sec, dir_data=dir_data, nc_prefix=nc_prefix, 
                            fn_match_pattern=fn_match_pattern, fn_dateformat=fn_dateformat, dir_out=dir_out, renamevars=renamevars)

vars_pd, dims_pd = get_ncvardimlist(file_nc=file_merged)

data_nc = Dataset(file_merged)
data_fromnc = data_nc.variables[renamevars['salinity']]
data_fromnc_all = data_fromnc[:]

varlist = data_nc.variables.keys()
print(varlist)
data_nc.close()
#data_to.variables['salinity'].setncattr('standard_name','sea_water_salinity')
#data_to.variables['water_temp'].setncattr('standard_name','sea_water_potential_temperature')
    