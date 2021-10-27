# -*- coding: utf-8 -*-
"""
Created on Wed Oct 27 15:13:32 2021

@author: veenstra
"""

import os
import matplotlib.pyplot as plt
plt.close('all')

from dfm_tools.io.polygon import Polygon
from dfm_tools.regulargrid import center2corner

dir_testinput = r'c:\DATA\dfm_tools_testdata'
dir_output = '.'

file_pol_list = [os.path.join(dir_testinput,'DFM_3D_z_Grevelingen\\geometry\\structures\\Grevelingen-FM_BL_fxw.pliz'),
                os.path.join(dir_testinput,'DFM_3D_z_Grevelingen\\geometry\\structures\\Grevelingen-FM_BL_fxw.pli'),
                os.path.join(dir_testinput,'world.ldb'),
                os.path.join(dir_testinput,'Maeslant.tek'),
                #'p:\\11205258-006-kpp2020_rmm-g6\\C_Work\\08_RMM_FMmodel\\geometry_j13_6-w3\\rmm_v1p3_fixed_weirs.pliz',
                #'p:\\11205258-006-kpp2020_rmm-g6\\C_Work\\08_RMM_FMmodel\\geometry_j13_6-w3\\structures\\rmm_v1p3_structures.pli',
                os.path.join(dir_testinput,'ballenplot\\0200a.tek'),
                os.path.join(dir_testinput,'ballenplot\\SDS-zd003b5dec2-sal.tek'),
                os.path.join(dir_testinput,'test_new.tek'),
                ]
                
for file_pol in file_pol_list:
    pol_data_list, pol_name_list, pol_comment_list = Polygon.fromfile(file_pol, pd_output=False)
    pol_data_pd_list = Polygon.fromfile(file_pol, pd_output=True)
    #pol_object_list = Polygon.fromfile(file_pol, obj_output=True)

    fig, ax = plt.subplots()
    for iP, pol_data in enumerate(pol_data_list):
        pd_collist = pol_data_pd_list[iP].columns.tolist()
        if 'datetime' in pd_collist:
            for iV in range(3,len(pd_collist)):
                ax.plot(pol_data_pd_list[iP].loc[:,'datetime'],pol_data_pd_list[iP].iloc[:,iV],'-',label=pd_collist[iV], linewidth=0.5)
        elif '*column 1 : X' in pd_collist or '*column 1 = x coordinate' in pd_collist:
            for iV in range(len(pd_collist)):
                ax.plot(pol_data_pd_list[iP].iloc[:,iV],'-',label=pd_collist[iV], linewidth=0.5)
            #retrieve again as tekal mapdata
            pol_data_tekmap, comment_list = Polygon.fromfile(file_pol, tekmap_output=True)
            fig_hnum = len(comment_list)//2
            fig_map, axs_map = plt.subplots(fig_hnum,2,figsize=(12,4+fig_hnum),sharex=True,sharey=True)
            axs_maplist = axs_map.reshape(axs_map.size)
            if '0200a' in file_pol:
                zcolnr = 2
            else:
                zcolnr = 1
            xcor = center2corner(pol_data_tekmap[:,:,0])
            zcor = center2corner(pol_data_tekmap[:,:,zcolnr])
            for iVar,comment in enumerate(comment_list):
                ax_map = axs_maplist[iVar]
                #pc = ax_map.pcolor(pol_data_tekmap[:,:,0],pol_data_tekmap[:,:,1],pol_data_tekmap[:,:,iVar],cmap='jet')
                pc = ax_map.pcolor(xcor,zcor,pol_data_tekmap[:,:,iVar],cmap='jet')
                ax_map.set_title(comment)
                fig_map.colorbar(pc, ax=ax_map)
                if iVar > len(axs_maplist)-3:
                    ax_map.set_xlabel(comment_list[0])
                if iVar%2 == 0:
                    ax_map.set_ylabel(comment_list[zcolnr])
            fig_map.tight_layout()
            fig_map.savefig(os.path.join(dir_output,'%s_mapvariables'%(os.path.basename(file_pol).replace('.',''))))
        else: #eg ldb
            ax.plot(pol_data[:,0],pol_data[:,1],'-',linewidth=0.5)
    ax.legend()
    fig.savefig(os.path.join(dir_output,os.path.basename(file_pol).replace('.','')))

