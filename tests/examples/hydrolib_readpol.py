# -*- coding: utf-8 -*-
"""
Created on Wed Aug 24 13:25:41 2022

@author: veenstra
"""

import os
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
plt.close('all')
from hydrolib.core.io.polyfile.models import PolyFile
from dfm_tools.hydrolib_helpers import polyobject_to_dataframe

dir_testinput = r'c:\DATA\dfm_tools_testdata'
dir_output = '.'

if 0: #read pli/pol/ldb files (tek files with 2/3 columns)
    file_pli_list = [Path(dir_testinput,'world.ldb'),
                     #Path(dir_testinput,r'GSHHS_f_L1_world_ldb_noaa_wvs.ldb'), #huge file, so takes a lot of time
                     Path(dir_testinput,'GSHHS_high_min1000km2.ldb'), #works but slow
                     #Path(dir_testinput,'DFM_3D_z_Grevelingen\\geometry\\structures\\Grevelingen-FM_BL_fxw.pli'),
                     Path(dir_testinput,'DFM_3D_z_Grevelingen\\geometry\\structures\\Grevelingen-FM_BL_fxw.pliz'), #results also in data property of Points (not only xy)
                     ]
    
    for file_pli in file_pli_list:
        #load boundary file
        polyfile_object = PolyFile(file_pli)
        
        fig,ax = plt.subplots()
        for iPO, pli_PolyObject_sel in enumerate(polyfile_object.objects):
            print(f'processing PolyObject {iPO+1} of {len(polyfile_object.objects)}: name={pli_PolyObject_sel.metadata.name}')
            polyobject_pd = polyobject_to_dataframe(pli_PolyObject_sel,dummy=999.999) #dummy is for world.ldb
            ax.plot(polyobject_pd['x'],polyobject_pd['y'])
        fig.savefig(os.path.join(dir_output,os.path.basename(file_pli).replace('.','')))
        
        #get extents of all objects in polyfile
        data_pol_pd_list = [polyobject_to_dataframe(polyobj) for polyobj in polyfile_object.objects]
        data_pol_pd_all = pd.concat(data_pol_pd_list)
        xmin,ymin = data_pol_pd_all[['x','y']].min()
        xmax,ymax = data_pol_pd_all[['x','y']].max()
        print(xmin,xmax,ymin,ymax)

    

if 1: #read tek files with more than 2 columns
    file_pli_list = [Path(dir_testinput,r'ballenplot\SDS-zd003b5dec2-sal.tek'), #TODO: UserWarning: Expected valid dimensions at line 14. (3D file). Request support for 3D file?
                     Path(dir_testinput,r'ballenplot\SDS-zd003b5dec2-sal_2D.tek'), #solved by removing 3rd dim, but than layers are sort of lost
                     #Path(dir_testinput,r'ballenplot\0200a.tek'), #TODO: UserWarning: Expected valid dimensions at line 6. (3D file). Request support for 3D file?
                     Path(dir_testinput,r'ballenplot\0200a_2D.tek'), #solved by removing 3rd dim, but than layers are sort of lost
                     #Path(dir_testinput,r'Gouda.tek'), #works (but slow since it is a large file)
                     Path(dir_testinput,r'Maeslant.tek'), #works
                     #Path(dir_testinput,r'ballenplot\nima-1013-lo-wl.tek'), # UserWarning: Expected a valid name or description at line 3. (name contains spaces and brackets)
                     Path(dir_testinput,r'ballenplot\nima-1013-lo-wl_validname.tek'), # works
                     Path(dir_testinput,r'test.tek'), # works
                     ]
    
    for file_pli in file_pli_list:
        if ('SDS-zd003b5dec2-sal' in str(file_pli)) or ('0200a' in str(file_pli)):
            convert_xy_to_time = False
        else:
            convert_xy_to_time = True
        
        #load pol/tek/pli/ldb file
        polyfile_object = PolyFile(file_pli)
        
        fig,ax = plt.subplots()
        for iPO, pli_PolyObject_sel in enumerate(polyfile_object.objects):
            print(f'processing PolyObject {iPO+1} of {len(polyfile_object.objects)}: name={pli_PolyObject_sel.metadata.name}')
            polyobject_pd = polyobject_to_dataframe(pli_PolyObject_sel,convert_xy_to_time=convert_xy_to_time)
            print(pli_PolyObject_sel.metadata)
            print(pli_PolyObject_sel.description.content)
            if convert_xy_to_time:
                ax.plot(polyobject_pd)
                ax.set_title(pli_PolyObject_sel.metadata.name)
            else: #this is only for datasets that can currently not be plotted nicely. Not really the responsability of hydrolib I presume
                ax.scatter(polyobject_pd['x'],polyobject_pd['y'],c=polyobject_pd[0]) #TODO: valuable to be able to plot this nicely again?
                ax.set_title(f'inconvenient plot of {pli_PolyObject_sel.metadata.name}')
        fig.savefig(os.path.join(dir_output,os.path.basename(file_pli).replace('.','')))

    
    