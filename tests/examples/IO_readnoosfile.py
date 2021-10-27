# -*- coding: utf-8 -*-
"""
Created on Wed Oct 27 22:10:41 2021

@author: veenstra
"""
import os
import matplotlib.pyplot as plt
plt.close('all')

from dfm_tools.io.noos import read_noosfile, write_noosfile

dir_testinput = r'c:\DATA\dfm_tools_testdata'
dir_output = '.'


file_noos = os.path.join(dir_testinput,'KORNWDZBTN_waterlevel_20061201_20190101_diffnanvals.noos')
#import matplotlib.pyplot as plt

noosdata_pd, noosheader_dict = read_noosfile(file_noos=file_noos,get_header=True,na_values=['N/A',-999])
"""
fig, ax = plt.subplots()
ax.plot(noosdata_pd['datetime'], noosdata_pd['values'])
"""
#write noosfile
header_dict = {'Source': 'from noosfile', 'Analysis time': 200201012300, 'Time zone': 'GMT', 'Coordinate system': 'WGS84',
               'x-coordinate': 200, 'y-coordinate': 300, 'Variable': 'waterlevel', 'Unit': 'm', 'teststring': ''}

write_noosfile(os.path.join(dir_output,'noos_out.noos'), noosdata_pd, metadata=noosheader_dict, na_values=-999.999)
write_noosfile(os.path.join(dir_output,'noos_out2.noos'), noosdata_pd, metadata=header_dict, na_values=-999.999)
