# -*- coding: utf-8 -*-
"""
Created on Wed Oct 27 14:52:44 2021

@author: veenstra
"""

import os
import numpy as np
import matplotlib.pyplot as plt
plt.close('all')

from dfm_tools.get_nc import get_ncmodeldata, plot_netmapdata#, get_xzcoords_onintersection
from dfm_tools.regulargrid import meshgridxy2verts, center2corner
from dfm_tools.get_nc_helpers import get_ncvardimlist

dir_testinput = r'c:\DATA\dfm_tools_testdata'
dir_output = '.'

file_nc_list = [r'p:\11203869-morwaqeco3d\05-Tidal_inlet\02_FM_201910\FM_MF10_Max_30s\wave\wavm-inlet.nc',
                r'p:\11200665-c3s-codec\2_Hydro\ECWMF_meteo\meteo\ERA-5\2000\ERA5_metOcean_atm_19991201_19991231.nc',
                #r'p:\1204257-dcsmzuno\2014\data\meteo\HIRLAM72_2018\h72_201803.nc', #TODO: xarray MissingDimensionsError
                r'p:\11202255-sfincs\Testbed\Original_tests\01_Implementation\08_restartfile\sfincs_map.nc', #not available anymore
                ]

for file_nc in file_nc_list:
    
    #get cell center coordinates from regular grid
    if 'ERA5_metOcean_atm' in file_nc:
        data_fromnc_x_1D = get_ncmodeldata(file_nc=file_nc, varname='longitude')
        data_fromnc_y_1D = get_ncmodeldata(file_nc=file_nc, varname='latitude')
        data_fromnc_x, data_fromnc_y = np.meshgrid(data_fromnc_x_1D, data_fromnc_y_1D)
    else:
        data_fromnc_x = get_ncmodeldata(file_nc=file_nc, varname='x')
        data_fromnc_y = get_ncmodeldata(file_nc=file_nc, varname='y')
    
    vars_pd, dims_pd = get_ncvardimlist(file_nc=file_nc)
    x_cen_withbnd = center2corner(data_fromnc_x)
    y_cen_withbnd = center2corner(data_fromnc_y)
    grid_verts = meshgridxy2verts(x_cen_withbnd, y_cen_withbnd)
        
    fig, axs = plt.subplots(2,1, figsize=(10,9))
    ax = axs[0]
    ax.set_title('xy center data converted to xy corners')
    ax.plot(data_fromnc_x,data_fromnc_y, linewidth=0.5, color='blue')
    ax.plot(data_fromnc_x.T,data_fromnc_y.T, linewidth=0.5, color='blue')
    ax.plot(x_cen_withbnd,y_cen_withbnd, linewidth=0.5, color='crimson')
    ax.plot(x_cen_withbnd.T,y_cen_withbnd.T, linewidth=0.5, color='crimson')
    ax.set_aspect('equal')
    ax = axs[1]
    ax.set_title('xy corner data converted to vertices (useful for map plotting)')
    plot_netmapdata(grid_verts, values=None, ax=ax, linewidth=0.5, color='crimson', facecolor='None')
    ax.set_aspect('equal')
    fig.tight_layout()
    plt.savefig(os.path.join(dir_output,'%s_grid'%(os.path.basename(file_nc).replace('.',''))))
    


    