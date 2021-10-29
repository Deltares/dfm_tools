# -*- coding: utf-8 -*-
"""
Created on Wed Oct 27 14:21:28 2021

@author: veenstra
"""

import os
import matplotlib.pyplot as plt
plt.close('all')
import numpy as np
import datetime as dt

from dfm_tools.get_nc import get_netdata, get_ncmodeldata, get_xzcoords_onintersection, plot_netmapdata
from dfm_tools.io.polygon import LineBuilder#, Polygon

dir_testinput = r'c:\DATA\dfm_tools_testdata'
dir_output = '.'

file_nc_list = [os.path.join(dir_testinput,'DFM_sigma_curved_bend\\DFM_OUTPUT_cb_3d\\cb_3d_map.nc'), #sigmalayer
            os.path.join(dir_testinput,'DFM_3D_z_Grevelingen','computations','run01','DFM_OUTPUT_Grevelingen-FM','Grevelingen-FM_0000_map.nc'), #zlayer
            r'p:\1204257-dcsmzuno\2006-2012\3D-DCSM-FM\A18b_ntsu1\DFM_OUTPUT_DCSM-FM_0_5nm\DCSM-FM_0_5nm_0000_map.nc', #fullgrid
            r'p:\11203379-005-mwra-updated-bem\03_model\03_runs\A86_erik\DFM_OUTPUT_MB_02\MB_02_0000_map.nc', #fullgrid
            r'p:\11205258-006-kpp2020_rmm-g6\C_Work\08_RMM_FMmodel\computations\run_183\DFM_OUTPUT_RMM_dflowfm\RMM_dflowfm_0000_map.nc', #2D model
            ]

for file_nc in file_nc_list:    
    if 'cb_3d_map' in file_nc:
        timestep = 72
        layno = 5
        calcdist_fromlatlon = None
        multipart = None
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
        #optimize_dist = None
        #data_frommap_bl = get_ncmodeldata(file_nc=file_nc, varname='mesh2d_ucmag', timestep=10, layer=0, multipart=multipart)
    elif 'Grevelingen' in file_nc:
        timestep = 3
        layno = 35
        calcdist_fromlatlon = None
        multipart = None
        line_array = np.array([[ 56267.59146475, 415644.67447155],
                               [ 64053.73427496, 419407.58239502]])
        line_array = np.array([[ 53181.96942503, 424270.83361629],
                               [ 55160.15232593, 416913.77136685]])
        #line_array = np.array([[ 52787.21854294, 424392.10414528],
        #                       [ 55017.72655174, 416403.77313703],
        #                       [ 65288.43784807, 419360.49305567]])
        val_ylim = [-25,5]
        clim_bl = None
        #optimize_dist = 150
    elif '-mwra-' in file_nc:
        timestep = 30
        layno = 5
        calcdist_fromlatlon = True
        multipart = None
        #provide xy order, so lonlat
        line_array = np.array([[-71.10395926,  42.3404146 ],
                               [-69.6762489 ,  42.38341792]])
        #line_array = np.array([[-69.95107235,  41.75440072], #dummy for partion 0000
        #   [-70.32764734,  42.04059772]])
        val_ylim = [-400,10]
        clim_bl = None
        #optimize_dist = None
    elif 'DCSM-FM_0_5nm' in file_nc:
        timestep = 365
        layno = 5
        calcdist_fromlatlon = True
        multipart = None
        #provide xy order, so lonlat
        line_array = np.array([[ 0.97452229, 51.13407643],
                               [ 1.89808917, 50.75191083]])
        line_array = np.array([[10.17702481, 57.03663877], #dummy for partition 0000
                               [12.38583134, 57.61284917]])
        line_array = np.array([[ 8.92659074, 56.91538014],
                               [ 8.58447136, 58.66874192]])
        val_ylim = [-600,1]
        clim_bl = [-500,0]
        #optimize_dist = 0.1
    elif 'DFM_OUTPUT_RMM_dflowfm' in file_nc:
        timestep = 365
        layno = None
        calcdist_fromlatlon = None
        multipart = None
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
        #optimize_dist = 150
    else:
        raise Exception('ERROR: no settings provided for this mapfile')
    
    data_frommap_bl = get_ncmodeldata(file_nc=file_nc, varname='mesh2d_flowelem_bl', multipart=multipart)
   
    ugrid = get_netdata(file_nc=file_nc, multipart=multipart)
    #get bed layer
    
    #create plot with ugrid and cross section line
    fig, ax_input = plt.subplots()
    pc = plot_netmapdata(ugrid.verts, values=data_frommap_bl, ax=ax_input, linewidth=0.5, edgecolors='face', cmap='jet')#, color='crimson', facecolor="None")
    pc.set_clim(clim_bl)
    fig.colorbar(pc, ax=ax_input)
    ax_input.set_aspect('equal')
    if 0: #click interactive polygon
        #pol_frominput = Polygon.frominteractive(ax) #this is old, does not work, use code below
        line, = ax_input.plot([], [],'o-')  # empty line
        linebuilder = LineBuilder(line) #after this click your line and then run the line below
        line_array = linebuilder.line_array
    ax_input.plot(line_array[0,0],line_array[0,1],'bx',linewidth=3,markersize=10)
    ax_input.plot(line_array[:,0],line_array[:,1],'b',linewidth=3)
    
    
    runtime_tstart = dt.datetime.now() #start timer
    #intersect function, find crossed cell numbers (gridnos) and coordinates of intersection (2 per crossed cell)
    intersect_gridnos, intersect_coords = ugrid.polygon_intersect(line_array, optimize_dist=False)
    #derive vertices from cross section (distance from first point)
    crs_verts, crs_plotdata = get_xzcoords_onintersection(file_nc=file_nc, varname='mesh2d_sa1', line_array=line_array, intersect_gridnos=intersect_gridnos, intersect_coords=intersect_coords, timestep=timestep, calcdist_fromlatlon=calcdist_fromlatlon, multipart=multipart)
    #for iIC, int_coord in enumerate(intersect_coords):
    #    #print(int_coord)
    #    ax_input.plot(int_coord[0,:],int_coord[1,:],'ro-')
    
    #plot crossed cells (gridnos) in first plot
    #print(layno)#data_frommap_flat = data_frommap[0,intersect_gridnos,layno]
    #pc = plot_netmapdata(ugrid.verts[intersect_gridnos,:,:], values=data_frommap_flat, ax=ax_input, linewidth=0.5, cmap="jet")
    plt.savefig(os.path.join(dir_output,'%s_gridbed'%(os.path.basename(file_nc).replace('.',''))))

    fig, ax = plt.subplots()
    pc = plot_netmapdata(crs_verts, values=crs_plotdata, ax=ax, linewidth=0.5, cmap='jet')#, edgecolor='k')
    fig.colorbar(pc, ax=ax)
    ax.set_ylim(val_ylim)
    plt.savefig(os.path.join(dir_output,'%s_crossect'%(os.path.basename(file_nc).replace('.',''))))
    
    runtime_tstop = dt.datetime.now()
    runtime_timedelta = (runtime_tstop-runtime_tstart).total_seconds()
    print('calculating and plotting cross section finished in %.1f seconds'%(runtime_timedelta))

