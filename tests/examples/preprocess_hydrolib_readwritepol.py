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
import dfm_tools as dfmt
import hydrolib.core.dflowfm as hcdfm


dir_testinput = r'c:\DATA\dfm_tools_testdata'
dir_output = '.'

dtstart = dt.datetime.now()

write_outfile = False

file_pli_list = [Path(dir_testinput,'world.ldb'),
                 #Path(dir_testinput,r'GSHHS_f_L1_world_ldb_noaa_wvs.ldb'), #huge file, so takes a lot of time
                 Path(dir_testinput,'GSHHS_high_min1000km2.ldb'), #works but slow
                 #Path(dir_testinput,'DFM_3D_z_Grevelingen\\geometry\\structures\\Grevelingen-FM_BL_fxw.pli'),
                 Path(dir_testinput,'DFM_3D_z_Grevelingen\\geometry\\structures\\Grevelingen-FM_BL_fxw.pliz'), #results also in data property of Points (not only xy)
                 #Path(dir_testinput,r'ballenplot\SDS-zd003b5dec2-sal.tek'), #TODO: UserWarning: Expected valid dimensions at line 14. (3D file). Request support for 3D file?
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
    #load pol/tek/pli/ldb file
    polyfile_object = hcdfm.PolyFile(file_pli)
    
    #empty polyfile object to append polyobjects to for testing full read/write workflow
    if write_outfile:
        polyfile_object_out = hcdfm.PolyFile()
    
    fig,ax = plt.subplots()
    for iPO, pli_PolyObject_sel in enumerate(polyfile_object.objects):
        print(f'processing PolyObject {iPO+1} of {len(polyfile_object.objects)}: name={pli_PolyObject_sel.metadata.name}')
        
        #conversion to dataframe
        polyobject_pd = dfmt.pointlike_to_DataFrame(pli_PolyObject_sel)
        if 'world.ldb' in str(file_pli):
            polyobject_pd[polyobject_pd==999.999] = np.nan #for world.ldb
        
        #get content
        if pli_PolyObject_sel.description is None:
            content_str, content = '', None #content none is necessary otherwise empty comment line added
        else:
            content_str = content = pli_PolyObject_sel.description.content
        
        #collect for writing outfile
        if write_outfile:
            polyobject_out = dfmt.DataFrame_to_PolyObject(polyobject_pd, name=pli_PolyObject_sel.metadata.name, content=content)
            polyfile_object_out.objects.append(polyobject_out)

        ax.set_title(f'{len(polyfile_object.objects)} PolyObjects, name of first is {pli_PolyObject_sel.metadata.name}')
        #plotting
        if 'Date' in content_str: #conversion of xy to datetime
            polyobject_pd_timeidx = dfmt.parse_xy_to_datetime(polyobject_pd) #TODO: not suported the other way round, necessary to add?
            ax.plot(polyobject_pd_timeidx)
            ax.legend(content_str.split('\n')[2:])
        elif 'SDS-zd003b5dec2-sal' in str(file_pli): #this is a 3D tekfile, depends on the contents how to handle it so it is hardcoded
            plt.close()
            fig,ax = plt.subplots(4,2,figsize=(12,7),sharex=True)
            for colno in range(8):
                ax.flatten()[colno].scatter(polyobject_pd['x'],polyobject_pd['y'],c=polyobject_pd[colno],cmap='jet')
                colmeta = content_str.split('\n')[colno+5]
                ax.flatten()[colno].set_title(f'{pli_PolyObject_sel.metadata.name}: {colmeta}')
        elif '0200a' in str(file_pli): #this is a 3D tekfile, depends on the contents how to handle it so it is hardcoded
            ax.scatter(polyobject_pd['x'],polyobject_pd[0],c=polyobject_pd[1])
            ax.set_xlabel('x')
            ax.set_ylabel('z')
        else: #pli/pol/ldb files with two columns
            ax.plot(polyobject_pd['x'],polyobject_pd['y'])
    
    fig.tight_layout()
    fig.savefig(os.path.join(dir_output,os.path.basename(file_pli).replace('.','')))
    
    if write_outfile:
        polyfile_object_out.save(os.path.basename(file_pli).replace('.','_out.')) #TODO: better formatting of plifile (also more precision, maybe write xy/datetime as ints?)
    
    #get extents of all objects in polyfile
    data_pol_pd_list = [dfmt.pointlike_to_DataFrame(polyobj) for polyobj in polyfile_object.objects]
    data_pol_pd_all = pd.concat(data_pol_pd_list)
    xmin,ymin = data_pol_pd_all[['x','y']].min()
    xmax,ymax = data_pol_pd_all[['x','y']].max()
    print(xmin,xmax,ymin,ymax)


time_passed = (dt.datetime.now()-dtstart).total_seconds()
print(f'>>time passed: {time_passed:.2f} sec')

