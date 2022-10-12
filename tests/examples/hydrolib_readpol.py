# -*- coding: utf-8 -*-
"""
Created on Wed Aug 24 13:25:41 2022

@author: veenstra
"""

import os
from pathlib import Path
import numpy as np
import pandas as pd
import datetime as dt
import matplotlib.pyplot as plt
plt.close('all')
from hydrolib.core.io.polyfile.models import PolyFile
from dfm_tools.hydrolib_helpers import polyobject_to_dataframe, dataframe_to_polyobject

dir_testinput = r'c:\DATA\dfm_tools_testdata'
dir_output = '.'

dtstart = dt.datetime.now()


if 1: #read pli/pol/ldb files (tek files with 2/3 columns)
    file_pli_list = [Path(dir_testinput,'world.ldb'),
                     #Path(dir_testinput,r'GSHHS_f_L1_world_ldb_noaa_wvs.ldb'), #huge file, so takes a lot of time
                     Path(dir_testinput,'GSHHS_high_min1000km2.ldb'), #works but slow
                     #Path(dir_testinput,'DFM_3D_z_Grevelingen\\geometry\\structures\\Grevelingen-FM_BL_fxw.pli'),
                     Path(dir_testinput,'DFM_3D_z_Grevelingen\\geometry\\structures\\Grevelingen-FM_BL_fxw.pliz'), #results also in data property of Points (not only xy)
                     ]
    
    for file_pli in file_pli_list:
        #load boundary file
        polyfile_object = PolyFile(file_pli)
        
        #empty polyfile object to append polyobjects to for testing full read/write workflow
        polyfile_object_out = PolyFile()
        fig,ax = plt.subplots()
        for iPO, pli_PolyObject_sel in enumerate(polyfile_object.objects):
            print(f'processing PolyObject {iPO+1} of {len(polyfile_object.objects)}: name={pli_PolyObject_sel.metadata.name}')
            polyobject_pd = polyobject_to_dataframe(pli_PolyObject_sel)
            polyobject_pd[polyobject_pd==999.999] = np.nan #for world.ldb
            ax.plot(polyobject_pd['x'],polyobject_pd['y'])
            if pli_PolyObject_sel.description is None:
                content = None
            else:
                content = pli_PolyObject_sel.description.content
            #polyobject_out = dataframe_to_polyobject(polyobject_pd, name=pli_PolyObject_sel.metadata.name, content=content)
            #polyfile_object_out.objects.append(polyobject_out)
        fig.savefig(os.path.join(dir_output,os.path.basename(file_pli).replace('.','')))
        #polyfile_object_out.save(os.path.basename(file_pli).replace('.','_out.')) #TODO: better formatting of plifile
        
        #get extents of all objects in polyfile
        data_pol_pd_list = [polyobject_to_dataframe(polyobj) for polyobj in polyfile_object.objects]
        data_pol_pd_all = pd.concat(data_pol_pd_list)
        xmin,ymin = data_pol_pd_all[['x','y']].min()
        xmax,ymax = data_pol_pd_all[['x','y']].max()
        print(xmin,xmax,ymin,ymax)
    

if 0: #read tek files with more than 2 columns
    file_pli_list = [#Path(dir_testinput,r'ballenplot\SDS-zd003b5dec2-sal.tek'), #TODO: UserWarning: Expected valid dimensions at line 14. (3D file). Request support for 3D file?
                     Path(dir_testinput,r'ballenplot\SDS-zd003b5dec2-sal_2D.tek'), #solved by removing 3rd dim, but than layers are sort of lost
                     #Path(dir_testinput,r'ballenplot\0200a.tek'), #TODO: UserWarning: Expected valid dimensions at line 6. (3D file). Request support for 3D file?
                     Path(dir_testinput,r'ballenplot\0200a_2D.tek'), #solved by removing 3rd dim, but than layers are sort of lost
                     #Path(dir_testinput,r'Gouda.tek'), #works (but slow since it is a large file)
                     Path(dir_testinput,r'Maeslant.tek'), #works
                     #Path(dir_testinput,r'ballenplot\nima-1013-lo-wl.tek'), # UserWarning: Expected a valid name or description at line 3. (name contains spaces and brackets)
                     Path(dir_testinput,r'ballenplot\nima-1013-lo-wl_validname.tek'), # works
                     Path(dir_testinput,r'test_new.tek'), # works
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
            polyobject_pd = polyobject_to_dataframe(pli_PolyObject_sel,convert_xy_to_time=convert_xy_to_time) #TODO: convert_xy_to_time is currently not supported by dataframe_to_polyobject()
            print(pli_PolyObject_sel.metadata)
            print(pli_PolyObject_sel.description.content)
            
            if convert_xy_to_time:
                ax.plot(polyobject_pd)
                ax.legend(pli_PolyObject_sel.description.content.split('\n')[2:])
                ax.set_title(pli_PolyObject_sel.metadata.name)
            elif 'SDS-zd003b5dec2-sal' in str(file_pli): #this is a 3D tekfile, depends on the contents how to handle it so it is hardcoded
                plt.close()
                fig,ax = plt.subplots(4,2,figsize=(12,7),sharex=True)
                for colno in range(8):
                    ax.flatten()[colno].scatter(polyobject_pd['x'],polyobject_pd['y'],c=polyobject_pd[colno],cmap='jet')
                    colmeta = pli_PolyObject_sel.description.content.split('\n')[colno+5]
                    ax.flatten()[colno].set_title(f'{pli_PolyObject_sel.metadata.name}: {colmeta}')
            elif '0200a' in str(file_pli): #this is a 3D tekfile, depends on the contents how to handle it so it is hardcoded
                ax.scatter(polyobject_pd['x'],polyobject_pd[0],c=polyobject_pd[1])
                ax.set_xlabel('x')
                ax.set_ylabel('z')
                ax.set_title(pli_PolyObject_sel.metadata.name)
            else: #this is not used with the current examples
                ax.scatter(polyobject_pd['x'],polyobject_pd['y'],c=polyobject_pd[0])
                ax.set_title(f'inconvenient plot of {pli_PolyObject_sel.metadata.name}')
        fig.tight_layout()
        fig.savefig(os.path.join(dir_output,os.path.basename(file_pli).replace('.','')))



if 0: #write pol/pli
    poly_pd = pd.DataFrame(np.array([[57.9730136 , 24.20954011,0],
                                    [57.49360017, 24.22389666,1],
                                    [57.12859221, 24.34592732,2],
                                    [57.01418674, 25.05657643,3],
                                    [57.29747649, 25.30063775,3],
                                    [57.864056  , 25.32217257,4]]),columns=['x','y',0])
    poly_pd = pd.DataFrame(np.array([[57.9730136 , 24.20954011],
                                    [57.49360017, 24.22389666],
                                    [57.12859221, 24.34592732],
                                    [57.01418674, 25.05657643],
                                    [57.29747649, 25.30063775],
                                    [57.864056  , 25.32217257]]),columns=['x','y'])
    polyfile_object_out = PolyFile()
    polyobject = dataframe_to_polyobject(poly_pd,name='hycom_pol')
    polyfile_object_out.objects.append(polyobject)
    polyfile_object_out.save('hycom.pli') #TODO: better formatting of plifile (also more precision)


time_passed = (dt.datetime.now()-dtstart).total_seconds()
print(f'>>time passed: {time_passed:.2f} sec')

