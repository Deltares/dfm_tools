# -*- coding: utf-8 -*-
"""
Created on Tue Oct 18 15:20:06 2022

@author: veenstra
"""

import os
from dfm_tools.download import download_ERA5

basedir = ''

# domain
longitude_min = -5
longitude_max = 4
latitude_min = 41.5
latitude_max = 42

#tstart and tstop as understood by pd.date_range with freq='MS' (month start)
tstart = '2021-01'
tstop = '2021-02'

# variables [short names, long names]
variables = ['v10n']

for varkey in variables:
    dir_out = os.path.join(basedir,varkey)
    if not os.path.isdir(dir_out):
        os.mkdir(dir_out)
    
    download_ERA5(varkey, tstart=tstart, tstop=tstop,
                  longitude_min=longitude_min, longitude_max=longitude_max, latitude_min=latitude_min, latitude_max=latitude_max)

    
