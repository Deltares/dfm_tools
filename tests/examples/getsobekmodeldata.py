# -*- coding: utf-8 -*-
"""
Created on Wed Oct 27 14:47:34 2021

@author: veenstra

this test retrieves sobek observation data and plots it
"""

import os
import matplotlib.pyplot as plt
plt.close('all')

from dfm_tools.get_nc import get_ncmodeldata
#from dfm_tools.get_nc_helpers import get_hisstationlist

dir_testinput = r'c:\DATA\dfm_tools_testdata'
dir_output = '.'

file_nc = os.path.join(dir_testinput,'KenmerkendeWaarden','observations.nc')

#station_names = get_hisstationlist(file_nc=file_nc, varname='observation_id')
#station_names = get_hisstationlist(file_nc=file_nc, varname='water_level')
data_fromsobek = get_ncmodeldata(file_nc=file_nc, varname='water_level', station=['Maasmond','HKVHLD','MO_1035.00'], timestep='all')

fig, ax = plt.subplots()
for iL in range(data_fromsobek.shape[1]):
    ax.plot(data_fromsobek.var_times,data_fromsobek[:,iL],'-', label=data_fromsobek.var_stations.iloc[iL,0])
ax.legend()
plt.savefig(os.path.join(dir_output,'%s_waterlevel'%(os.path.basename(file_nc).replace('.',''))))




    
    
