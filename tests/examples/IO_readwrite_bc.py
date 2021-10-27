# -*- coding: utf-8 -*-
"""
Created on Wed Oct 27 22:11:54 2021

@author: veenstra
"""

import os
import datetime as dt
import pandas as pd
import matplotlib.pyplot as plt
plt.close('all')

from dfm_tools.io.polygon import Polygon
from dfm_tools.io.bc import read_bcfile, write_bcfile    

dir_testinput = r'c:\DATA\dfm_tools_testdata'
dir_output = '.'


#RMM: read tekfile and convert to bc file
file_outbc = os.path.join(dir_output,'lateralen_RMM.bc');

refdate = dt.datetime(2007,1,1)

file_tek = os.path.join(dir_testinput,'Gouda.tek');

data_tek_pd = Polygon.fromfile(file_tek,pd_output=True)[0]

datesel_bool = (data_tek_pd['datetime'] >= refdate) & (data_tek_pd['datetime'] <= dt.datetime(2014,1,1))
data_tek_pd_sel = data_tek_pd[datesel_bool]

data_pd_tobc = pd.DataFrame({('time','datetime'):data_tek_pd_sel['datetime'], ('lateral_discharge','m3/s'): data_tek_pd_sel.iloc[:,3]})

metadata_dict = {'Name':'HY_0.77_C_ONB_Gouda', 'Function':'timeseries', 'Time-interpolation':'linear'}
write_bcfile(filename=file_outbc, datablocks=data_pd_tobc, metadatas=metadata_dict, refdate=refdate, tzone=1, float_format='%6.2f')


#read and write bc file
bc_in = os.path.join(dir_testinput,'lateralen_20072008_maas.bc')
bc_out = os.path.join(dir_output,'lateralen_20072008_maas_dfmtoolscopy.bc')

datablock_list, metadata_list = read_bcfile(filename=bc_in)

write_bcfile(filename=bc_out, datablocks=datablock_list, metadatas=metadata_list, float_format='%.2f')


