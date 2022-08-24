# -*- coding: utf-8 -*-
"""
Created on Wed Aug 24 13:25:41 2022

@author: veenstra
"""

import os
from pathlib import Path
import datetime as dt
import pandas as pd
#from netCDF4 import num2date
#import cftime
import matplotlib.pyplot as plt
plt.close('all')
import numpy as np
from hydrolib.core.io.polyfile.models import (
    #Description,
    #Metadata,
    #Point,
    PolyFile,
    #PolyObject,
)
from hydrolib.core.io.polyfile.parser import read_polyfile #TODO: should be replaced with PolyFile above

dir_testinput = r'c:\DATA\dfm_tools_testdata'
dir_output = '.'

    
if 1: #read pli/pol/ldb files (tek files with 2/3 columns)
    file_pli_list = [Path(dir_testinput,'world.ldb'),
                     #Path(dir_testinput,r'GSHHS_f_L1_world_ldb_noaa_wvs.ldb'), #huge file, so takes a lot of time
                     Path(dir_testinput,'GSHHS_high_min1000km2.ldb'), #works but slow
                     #Path(dir_testinput,'DFM_3D_z_Grevelingen\\geometry\\structures\\Grevelingen-FM_BL_fxw.pli'),
                     Path(dir_testinput,'DFM_3D_z_Grevelingen\\geometry\\structures\\Grevelingen-FM_BL_fxw.pliz'), #results also in data property of Points (not only xy)
                     ]
    
    for file_pli in file_pli_list:
        #load boundary file
        polyfile_object = read_polyfile(file_pli,has_z_values=False)
        pli_PolyObjects = polyfile_object['objects']
        
        fig,ax = plt.subplots()
        for iPO, pli_PolyObject_sel in enumerate(pli_PolyObjects):
            print(f'processing PolyObject {iPO+1} of {len(pli_PolyObjects)}: name={pli_PolyObject_sel.metadata.name}')
            xvals = np.array([p.x for p in pli_PolyObject_sel.points])
            yvals = np.array([p.y for p in pli_PolyObject_sel.points])
            if 'world.ldb' in str(file_pli):
                xvals[xvals==999.999] = np.nan
                yvals[yvals==999.999] = np.nan
            ax.plot(xvals,yvals)
            for iP, pli_Point_sel in []:#enumerate(pli_PolyObject_sel.points): #looping over plipoints within component loop, append to datablock_pd_allcomp
                print(f'processing Point {iP+1} of {len(pli_PolyObject_sel.points)}: ',end='')
                lonx, laty = pli_Point_sel.x, pli_Point_sel.y
                print(f'(x={lonx}, y={laty})')
                pli_PolyObject_name_num = f'{pli_PolyObject_sel.metadata.name}_{iP+1:04d}'
        fig.savefig(os.path.join(dir_output,os.path.basename(file_pli).replace('.','')))
    

if 1: #read tek files with more than 2 columns
    #TODO: read_polyfile does not read comments
    file_pli_list = [#Path(dir_testinput,r'ballenplot\SDS-zd003b5dec2-sal.tek'), #TODO: UserWarning: Expected valid dimensions at line 14. (3D file)
                     Path(dir_testinput,r'ballenplot\SDS-zd003b5dec2-sal_2D.tek'), #solved by removing 3rd dim, but than layers are sort of lost
                     #Path(dir_testinput,r'ballenplot\0200a.tek'), #TODO: no warning/error when reading file but result is empty
                     Path(dir_testinput,r'ballenplot\0200a_2D.tek'), #works but difficult to plot properly due to xyz-sal
                     #Path(dir_testinput,r'Gouda.tek'), #works (but slow since it is a large file)
                     Path(dir_testinput,r'Maeslant.tek'), #works
                     #Path(dir_testinput,r'ballenplot\nima-1013-lo-wl.tek'), #TODO: UserWarning: Expected a valid name or description at line 3.
                     Path(dir_testinput,r'ballenplot\nima-1013-lo-wl_validname.tek'), #removed spaces/brackets >> works
                     Path(dir_testinput,r'test.tek'), # works
                     ]
    
    for file_pli in file_pli_list:
        if ('SDS-zd003b5dec2-sal' in str(file_pli)) or ('0200a' in str(file_pli)):
            convert_time = False
        else:
            convert_time = True

        #load boundary file
        polyfile_object = read_polyfile(file_pli,has_z_values=False) #still false, since all then comes in data (not z)
        pli_PolyObjects = polyfile_object['objects']
    
        fig,ax = plt.subplots()
        for iPO, pli_PolyObject_sel in enumerate(pli_PolyObjects):
            print(f'processing PolyObject {iPO+1} of {len(pli_PolyObjects)}: name={pli_PolyObject_sel.metadata.name}')
            datavals = np.array([p.data for p in pli_PolyObject_sel.points])
            if convert_time:
                datetimevals = np.array([dt.datetime.strptime(f'{int(p.x):08d} {int(p.y):06d}','%Y%m%d %H%M%S') for p in pli_PolyObject_sel.points])
                ax.plot(datetimevals,datavals)
            else: #this is only for datasets that can currently not be plotted nicely. Not really the responsability of hydrolib I presume
                xvals = np.array([p.x for p in pli_PolyObject_sel.points])
                yvals = np.array([p.y for p in pli_PolyObject_sel.points])
                ax.scatter(xvals,yvals,c=datavals[:,0]) #TODO: valuable to be able to plot this nicely again?
        fig.savefig(os.path.join(dir_output,os.path.basename(file_pli).replace('.','')))

    
    