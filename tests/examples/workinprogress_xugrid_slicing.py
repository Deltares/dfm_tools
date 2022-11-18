# -*- coding: utf-8 -*-
"""
Created on Fri Nov 18 15:43:34 2022

@author: veenstra
"""


import os
import matplotlib.pyplot as plt
plt.close('all')
import numpy as np
import contextily as ctx
import datetime as dt
import xarray as xr
import dfm_tools as dfmt
import warnings

dir_testinput = r'c:\DATA\dfm_tools_testdata'
dir_output = '.'

file_nc_list = [os.path.join(dir_testinput,'DFM_sigma_curved_bend\\DFM_OUTPUT_cb_3d\\cb_3d_map.nc'), #sigmalayer
                os.path.join(dir_testinput,'DFM_3D_z_Grevelingen','computations','run01','DFM_OUTPUT_Grevelingen-FM','Grevelingen-FM_0000_map.nc'), #zlayer
                r'p:\1204257-dcsmzuno\2006-2012\3D-DCSM-FM\A18b_ntsu1\DFM_OUTPUT_DCSM-FM_0_5nm\DCSM-FM_0_5nm_0000_map.nc', #fullgrid
                r'p:\11206813-006-kpp2021_rmm-2d\C_Work\31_RMM_FMmodel\computations\model_setup\run_207\results\RMM_dflowfm_0000_map.nc', #2D model
                ]

for file_nc in file_nc_list:
    print('processing %s'%(os.path.basename(file_nc)))
    basename = os.path.basename(file_nc).replace('.','')
    
    if 'cb_3d_map' in file_nc:
        timestep = 72
        layno = 5
        calcdist_fromlatlon = None
        line_array = np.array([[ 185.08667065, 2461.11775254],
                               [2934.63837418, 1134.16019127]])
        line_array = np.array([[ 104.15421399, 2042.7077107 ],
                               [2913.47878063, 2102.48057382]])
        #line_array = np.array([[2084.67741935, 3353.02419355], #with linebend in cell en with line crossing same cell twice
        #   [2255.79637097, 3307.15725806],
        #   [2222.27822581, 3206.60282258],
        #   [2128.78024194, 3266.58266129]])
        val_ylim = None
        clim_bl = None
        clim_sal = None
        crs = None
        file_nc_fou = None
    elif 'Grevelingen' in file_nc:
        timestep = 3
        layno = 33 #35 is top
        calcdist_fromlatlon = None
        line_array = np.array([[ 56267.59146475, 415644.67447155],
                               [ 64053.73427496, 419407.58239502]])
        line_array = np.array([[ 53181.96942503, 424270.83361629],
                               [ 55160.15232593, 416913.77136685]])
        #line_array = np.array([[ 53181.96942503, 424270.83361629],
        #                       [ 55160.15232593, 416913.77136685],
        #                       [ 65288.15232593, 419360.77136685]])
        val_ylim = [-25,5]
        clim_bl = None
        clim_sal = [28,30.2]
        crs = "EPSG:28992"
        file_nc_fou = None
    elif 'DCSM-FM_0_5nm' in file_nc:
        timestep = 365
        layno = 45
        calcdist_fromlatlon = True
        #provide xy order, so lonlat
        line_array = np.array([[ 0.97452229, 51.13407643],
                               [ 1.89808917, 50.75191083]])
        line_array = np.array([[10.17702481, 57.03663877], #dummy for partition 0000
                               [12.38583134, 57.61284917]])
        line_array = np.array([[ 8.92659074, 56.91538014],
                               [ 8.58447136, 58.66874192]])
        val_ylim = [-600,1]
        clim_bl = [-500,0]
        clim_sal = [25,36]
        crs = "EPSG:4326"
        file_nc_fou = None
    elif 'RMM_dflowfm' in file_nc:
        timestep = 365 #50
        layno = None
        calcdist_fromlatlon = None
        #provide xy order, so lonlat
        line_array = np.array([[ 65655.72699961, 444092.54776465],
                               [ 78880.42720631, 435019.78832052]])
        line_array = np.array([[ 52444.56849912, 434039.27970214], #HVSL
                               [ 61304.25484967, 430703.86837017],
                               [ 62164.16558369, 428619.23628769]])
        line_array = np.array([[ 61013.8966525 , 446291.69129373], #NWW
                               [ 67151.68543524, 444096.96681991],
                               [ 69011.62143001, 442981.00522304],
                               [ 72210.71134101, 440302.69739058],
                               [ 74405.43581484, 438889.14603455],
                               [ 75632.99357138, 437401.19723874],
                               [ 79018.07708186, 435169.27404501],
                               [ 81324.39771538, 434536.89580679],
                               [ 82923.94267088, 434611.29324658],
                               [ 84449.09018659, 435132.07532512],
                               [ 86606.61594052, 434685.69068637],
                               [ 88689.74425466, 435355.26764449],
                               [ 90772.8725688 , 434983.28044554],
                               [ 91926.03288556, 435132.07532512]])
        val_ylim = None
        clim_bl = [-10,10]
        clim_sal = None
        crs = "EPSG:28992"
        file_nc_fou = os.path.join(dir_testinput,r'DFM_fou_RMM\RMM_dflowfm_0000_fou.nc')
    else:
        raise Exception('ERROR: no settings provided for this mapfile')
    
    file_nc_star = file_nc.replace('_0000_','_0*_') #TODO: make this the default
    data_xr = xr.open_dataset(file_nc)
    data_frommap_merged = dfmt.open_partitioned_dataset(file_nc_star)
    
    #get ugrid data, vars informatin and grid units (latter from bedlevel coordinates)
    vars_pd = dfmt.get_ncvarproperties(file_nc=file_nc)
    
    print('calculating and plotting cross section')
    runtime_tstart = dt.datetime.now() #start timer
    #intersect function, find crossed cell numbers (gridnos) and coordinates of intersection (2 per crossed cell)
    intersect_pd = dfmt.polygon_intersect(data_frommap_merged, line_array, optimize_dist=False, calcdist_fromlatlon=calcdist_fromlatlon)
    #derive vertices from cross section (distance from first point)
    crs_verts, crs_plotdata = dfmt.get_xzcoords_onintersection(data_frommap_merged, varname='mesh2d_sa1', intersect_pd=intersect_pd, timestep=timestep)
    fig, ax = plt.subplots()
    pc = dfmt.plot_netmapdata(crs_verts, values=crs_plotdata, ax=ax, cmap='jet')#, linewidth=0.5, edgecolor='k')
    fig.colorbar(pc, ax=ax)
    ax.set_ylim(val_ylim)
    plt.savefig(os.path.join(dir_output,f'{basename}_crossect'))
    runtime_timedelta = (dt.datetime.now()-runtime_tstart)
    print(f'calculating and plotting cross section finished in {runtime_timedelta}')