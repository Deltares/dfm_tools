# -*- coding: utf-8 -*-
"""
Created on Fri Mar 20 01:19:24 2020

@author: veenstra
"""

import pytest
import inspect
import os

dir_testinput = os.path.join(r'c:/DATA','dfm_tools_testdata')
from dfm_tools.testutils import getmakeoutputdir



@pytest.mark.acceptance
def test_workinprogress():
    ## WARNING: THIS TEST IS NOT YET FINISHED, WILL BE IMPROVED AND LINKED TO INTERNAL FUNCTIONS ASAP
    import os
    import matplotlib.pyplot as plt
    import numpy as np
    from netCDF4 import Dataset
    plt.close('all')
    
    from dfm_tools.get_nc import get_netdata, get_ncmodeldata, plot_netmapdata
    from dfm_tools.get_nc_helpers import get_ncvardimlist, get_hisstationlist#, get_varname_fromnc
    
    dir_output = getmakeoutputdir(__file__,inspect.currentframe().f_code.co_name)
    #dir_output = './test_output'

    # test Grevelingen (integrated example, where all below should move towards)
    file_nc = os.path.join(dir_testinput,r'DFM_3D_z_Grevelingen\computations\run01\DFM_OUTPUT_Grevelingen-FM\Grevelingen-FM_0000_map.nc')
    data_nc = Dataset(file_nc)
    ugrid = get_netdata(file_nc=file_nc)
    fig, ax = plt.subplots()
    plot_netmapdata(ugrid.verts, values=None, ax=None, linewidth=0.5, color="crimson", facecolor="None")
    ax.set_aspect('equal')
    
    #hirlam
    file_nc = r'p:\1204257-dcsmzuno\2014\data\meteo\HIRLAM72_2018\h72_201803.nc'
    vars_pd, dims_pd = get_ncvardimlist(file_nc=file_nc)
    data_nc = Dataset(file_nc)
    
    airp = data_nc.variables['air_pressure_fixed_height'][0,:,:]
    mesh2d_node_x = data_nc.variables['x'][:]
    mesh2d_node_y = data_nc.variables['y'][:]
    
    fig, ax = plt.subplots()
    ax.scatter(mesh2d_node_x,mesh2d_node_y,0.1,c='b')
    plt.savefig(os.path.join(dir_output,'hirlam_scatter'))

    fig, ax = plt.subplots()
    ax.plot(mesh2d_node_x,mesh2d_node_y,'-b',linewidth=0.2)
    ax.plot(mesh2d_node_x.T,mesh2d_node_y.T,'-b',linewidth=0.2)
    plt.savefig(os.path.join(dir_output,'hirlam_mesh'))

    fig, ax = plt.subplots()
    ax.pcolor(mesh2d_node_x,mesh2d_node_y,airp)
    #plt.pcolor(mesh2d_node_x,mesh2d_node_y,airp,linewidth=0.5)
    plt.savefig(os.path.join(dir_output,'hirlam_airp_pcolor'))
    

    #ERA5
    file_nc = r'p:\11200665-c3s-codec\2_Hydro\ECWMF_meteo\meteo\ERA-5\2000\ERA5_metOcean_atm_19991201_19991231.nc'
    
    data_lon = get_ncmodeldata(file_nc=file_nc, varname='longitude', multipart=False)
    data_lat = get_ncmodeldata(file_nc=file_nc, varname='latitude')
    data_psl = get_ncmodeldata(file_nc=file_nc, varname='msl',timestep=10, multipart=False)
    lons,lats = np.meshgrid(data_lon,data_lat)
    
    fig, ax = plt.subplots()
    ax.plot(lons, lats,'-b',linewidth=0.2)
    ax.plot(lons.T, lats.T,'-b',linewidth=0.2)
    plt.savefig(os.path.join(dir_output,'ERA5_mesh'))

    fig, ax = plt.subplots()
    ax.pcolor(lons, lats, data_psl[0,:,:])
    #plt.pcolor(mesh2d_node_x,mesh2d_node_y,airp,linewidth=0.5)
    plt.savefig(os.path.join(dir_output,'ERA5_msl_pcolor'))


    #SFINCS
    file_nc = r'p:\11202255-sfincs\Testbed\Original_runs\01_Implementation\14_restartfile\sfincs_map.nc'
    #file_nc = r'p:\11202255-sfincs\Testbed\Original_runs\03_Application\22_Tsunami_Japan_Sendai\sfincs_map.nc'
    
    data_fromnc_x = get_ncmodeldata(file_nc=file_nc, varname='x')
    data_fromnc_y = get_ncmodeldata(file_nc=file_nc, varname='y')
    data_fromnc_zs = get_ncmodeldata(file_nc=file_nc, varname='zs', timestep='all')
    fig, ax = plt.subplots()
    ax.plot(data_fromnc_x, data_fromnc_y,'-b',linewidth=0.2)
    ax.plot(data_fromnc_x.T, data_fromnc_y.T,'-b',linewidth=0.2)
    plt.savefig(os.path.join(dir_output,'SFINCS_mesh'))
    fig, ax = plt.subplots()
    ax.pcolor(data_fromnc_x, data_fromnc_y, data_fromnc_zs[0,:,:])
    plt.savefig(os.path.join(dir_output,'SFINCS_zs_pcolor'))

    data_fromnc_edgex = get_ncmodeldata(file_nc=file_nc, varname='edge_x')
    data_fromnc_edgey = get_ncmodeldata(file_nc=file_nc, varname='edge_y')
    data_fromnc_vmax = get_ncmodeldata(file_nc=file_nc, varname='u', timestep=10)
    fig, ax = plt.subplots()
    ax.plot(data_fromnc_edgex, data_fromnc_edgey,'-b',linewidth=0.2)
    ax.plot(data_fromnc_edgex.T, data_fromnc_edgey.T,'-b',linewidth=0.2)
    plt.savefig(os.path.join(dir_output,'SFINCS_meshedge'))
    fig, ax = plt.subplots()
    ax.pcolor(data_fromnc_edgex, data_fromnc_edgey, data_fromnc_vmax[0,:,:])
    plt.title('%s (%s)'%(data_fromnc_vmax.var_varname, data_fromnc_vmax.var_object.units))
    plt.savefig(os.path.join(dir_output,'SFINCS_u_pcolor'))

    #SFINCS HIS
    #file_nc = r'p:\11202255-sfincs\Testbed\Original_runs\01_Implementation\14_restartfile\sfincs_his.nc'
    file_nc = r'p:\11202255-sfincs\Testbed\Original_runs\03_Application\22_Tsunami_Japan_Sendai\sfincs_his.nc'
    
    station_names = get_hisstationlist(file_nc=file_nc)
    data_fromnc_his = get_ncmodeldata(file_nc=file_nc, varname='point_zs', station='all', timestep='all')

    fig, ax = plt.subplots()
    for iS,stat_name in enumerate(station_names):
        ax.plot(data_fromnc_his.var_times, data_fromnc_his[:,iS], label=stat_name)
    ax.legend()
    plt.savefig(os.path.join(dir_output,'SFINCS_hiszs'))
    





def test_trygetondepth():
    import numpy as np
    from netCDF4 import Dataset
    
    from dfm_tools.get_nc import get_ncmodeldata#, get_netdata
    from dfm_tools.get_nc_helpers import get_varname_fromnc
    
    #code from test_get_nc test d
    file_nc = os.path.join(dir_testinput,r'DFM_sigma_curved_bend\DFM_OUTPUT_cb_3d\cb_3d_map.nc')
    
    timestep = 72
    #layno = 5
    #calcdist_fromlatlon = None
    multipart = None
    #line_array = np.array([[ 104.15421399, 2042.7077107 ],
    #                       [2913.47878063, 2102.48057382]])
    #val_ylim = None
    #clim_bl = None
    #optimize_dist = None
    #ugrid = get_netdata(file_nc=file_nc, multipart=multipart)
    #intersect_gridnos, intersect_coords = ugrid.polygon_intersect(line_array, optimize_dist=None)
    
    
    #code from get_xzcoords_onintersection
    data_nc = Dataset(file_nc)
    
    varn_mesh2d_s1 = get_varname_fromnc(data_nc,'mesh2d_s1')
    data_frommap_wl3 = get_ncmodeldata(file_nc, varname=varn_mesh2d_s1, timestep=timestep, multipart=multipart)
    data_frommap_wl3 = data_frommap_wl3[0,:]
    #data_frommap_wl3_sel = data_frommap_wl3[0,intersect_gridnos]
    varn_mesh2d_flowelem_bl = get_varname_fromnc(data_nc,'mesh2d_flowelem_bl')
    data_frommap_bl = get_ncmodeldata(file_nc, varname=varn_mesh2d_flowelem_bl, multipart=multipart)
    #data_frommap_bl_sel = data_frommap_bl[intersect_gridnos]
    
    dimn_layer = get_varname_fromnc(data_nc,'nmesh2d_layer')
    if dimn_layer is None: #no layers, 2D model
        nlay = 1
    else:
        nlay = data_nc.dimensions[dimn_layer].size
    
    varn_layer_z = get_varname_fromnc(data_nc,'mesh2d_layer_z')
    if varn_layer_z is None:
        laytyp = 'sigmalayer'
        #zvals_cen = np.linspace(data_frommap_bl_sel,data_frommap_wl3_sel,nlay)
        #zvals_interface = np.linspace(data_frommap_bl_sel,data_frommap_wl3_sel,nlay+1)
        zvals_interface = np.linspace(data_frommap_bl,data_frommap_wl3,nlay+1)
    else:
        laytyp = 'zlayer'
        #zvals_cen = get_ncmodeldata(file_nc=file_map, varname='mesh2d_layer_z', lay='all')#, multipart=False)
        #zvals_interface = get_ncmodeldata(file_nc=file_map, varname='mesh2d_interface_z')#, multipart=False)
        zvals_interface = data_nc.variables['mesh2d_interface_z'][:]
    
    print(laytyp)
    depth = -1
    z_test_higher = np.argmax((zvals_interface > depth),axis=0)
    z_test_lower = np.argmin((zvals_interface < depth),axis=0)
    z_test_all = z_test_higher==z_test_lower
    







