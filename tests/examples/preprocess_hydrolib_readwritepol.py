# -*- coding: utf-8 -*-
"""
Created on Wed Aug 24 13:25:41 2022

@author: veenstra
"""

import os
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

file_pli_list = [r'p:\11208054-004-dcsm-fm\models\model_input\bnd_cond\pli\DCSM-FM_OB_all_20181108.pli',
                 os.path.join(dir_testinput,'GSHHS_high_min1000km2.ldb'),
                 os.path.join(dir_testinput,'DFM_3D_z_Grevelingen\\geometry\\structures\\Grevelingen-FM_BL_fxw.pliz'),
                 os.path.join(dir_testinput,r'Maeslant.tek'),
                 os.path.join(dir_testinput,r'nima-1013-lo-wl_validname.tek'),
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
        
        #get content
        if pli_PolyObject_sel.description is None:
            content_str, content = '', None #content none is necessary otherwise empty comment line added
        else:
            content_str = content = pli_PolyObject_sel.description.content
        
        #collect for writing outfile
        if write_outfile:
            polyobject_out = dfmt.DataFrame_to_PolyObject(polyobject_pd, name=pli_PolyObject_sel.metadata.name, content=content)
            polyfile_object_out.objects.append(polyobject_out)

        #plotting
        if 'Date' in content_str: #conversion of xy to datetime
            polyobject_pd_toplot = dfmt.parse_xy_to_datetime(polyobject_pd) #TODO: not suported the other way round, necessary to add?
            polyobject_pd_toplot.columns = content_str.split('\n')[2:] #to fix legend labels
            polyobject_pd_toplot.plot(ax=ax)
        else: #pli/pol/ldb files with two columns
            polyobject_pd.plot(ax=ax,x='x',y='y',legend=False)
    ax.set_title(f'{len(polyfile_object.objects)} PolyObjects, name of last is {pli_PolyObject_sel.metadata.name}')
    fig.tight_layout()
    fig.savefig(os.path.join(dir_output,os.path.basename(file_pli).replace('.','')))
    
    if write_outfile:
        polyfile_object_out.save(os.path.join(dir_output,os.path.basename(file_pli).replace('.','_out.')))
    
    #get extents of all objects in polyfile
    data_pol_pd_list = [dfmt.pointlike_to_DataFrame(polyobj) for polyobj in polyfile_object.objects]
    data_pol_pd_all = pd.concat(data_pol_pd_list)
    xmin,ymin = data_pol_pd_all[['x','y']].min()
    xmax,ymax = data_pol_pd_all[['x','y']].max()
    print(xmin,xmax,ymin,ymax)


time_passed = (dt.datetime.now()-dtstart).total_seconds()
print(f'>>time passed: {time_passed:.2f} sec')

