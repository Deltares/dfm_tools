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
    dir_output = getmakeoutputdir(__file__,inspect.currentframe().f_code.co_name)
    """
    dir_output = './test_output'
    """
    
    import os
    import matplotlib.pyplot as plt
    import numpy as np
    plt.close('all')
    
    from dfm_tools.get_nc import get_netdata, get_ncmodeldata, plot_netmapdata
    from dfm_tools.get_nc_helpers import get_ncvardimlist, get_hisstationlist#, get_varname_fromnc
    

    # test Grevelingen (integrated example, where all below should move towards)
    file_nc = os.path.join(dir_testinput,r'DFM_3D_z_Grevelingen\computations\run01\DFM_OUTPUT_Grevelingen-FM\Grevelingen-FM_0000_map.nc')
    ugrid = get_netdata(file_nc=file_nc)
    fig, ax = plt.subplots()
    plot_netmapdata(ugrid.verts, values=None, ax=None, linewidth=0.5, color="crimson", facecolor="None")
    ax.set_aspect('equal')
    
    #hirlam
    file_nc = r'p:\1204257-dcsmzuno\2014\data\meteo\HIRLAM72_2018\h72_201803.nc'
    vars_pd, dims_pd = get_ncvardimlist(file_nc=file_nc)
    
    mesh2d_node_x = get_ncmodeldata(file_nc=file_nc, varname='x')
    mesh2d_node_y = get_ncmodeldata(file_nc=file_nc, varname='y')
    data_v = get_ncmodeldata(file_nc=file_nc, varname='northward_wind',timestep=0)[0,:,:]
    data_u = get_ncmodeldata(file_nc=file_nc, varname='eastward_wind',timestep=0)[0,:,:]
    #airp = get_ncmodeldata(file_nc=file_nc, varname='air_pressure_fixed_height',timestep=0)[0,:,:]
    magn = np.sqrt(data_u**2 + data_v**2)
    
    fig, ax = plt.subplots()
    ax.plot(mesh2d_node_x,mesh2d_node_y,'-b',linewidth=0.2)
    ax.plot(mesh2d_node_x.T,mesh2d_node_y.T,'-b',linewidth=0.2)
    plt.savefig(os.path.join(dir_output,'hirlam_mesh'))

    
    fig, ax = plt.subplots()
    ax.pcolor(mesh2d_node_x,mesh2d_node_y,magn)
    #plt.pcolor(mesh2d_node_x,mesh2d_node_y,airp,linewidth=0.5)
    plt.savefig(os.path.join(dir_output,'hirlam_magn_pcolor'))
    

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
    vars_pd, dims_pd = get_ncvardimlist(file_nc=file_nc)
    
    data_fromnc_x = get_ncmodeldata(file_nc=file_nc, varname='x')
    data_fromnc_y = get_ncmodeldata(file_nc=file_nc, varname='y')
    data_fromnc_zs = get_ncmodeldata(file_nc=file_nc, varname='zs', timestep='all')

    fig, ax = plt.subplots()
    ax.plot(data_fromnc_x, data_fromnc_y,'-b',linewidth=0.2)
    ax.plot(data_fromnc_x.T, data_fromnc_y.T,'-b',linewidth=0.2)
    plt.savefig(os.path.join(dir_output,'SFINCS_mesh'))    

    fig, axs = plt.subplots(3,1, figsize=(14,9))
    for iT, timestep in enumerate([0,1,10]):
        ax=axs[iT]
        pc = ax.pcolor(data_fromnc_x, data_fromnc_y, data_fromnc_zs[timestep,:,:],cmap='jet')
        pc.set_clim([0,0.15])
        cbar = fig.colorbar(pc, ax=ax)
        cbar.set_label('%s (%s)'%(data_fromnc_zs.var_varname, data_fromnc_zs.var_object.units))
        ax.set_title('t=%d (%s)'%(timestep, data_fromnc_zs.var_times.loc[timestep]))
        ax.set_aspect('equal')
    fig.tight_layout()
    plt.savefig(os.path.join(dir_output,'SFINCS_zs_pcolor'))


    data_fromnc_edgex = get_ncmodeldata(file_nc=file_nc, varname='edge_x')
    data_fromnc_edgey = get_ncmodeldata(file_nc=file_nc, varname='edge_y')
    data_fromnc_u = get_ncmodeldata(file_nc=file_nc, varname='u', timestep='all')
    data_fromnc_v = get_ncmodeldata(file_nc=file_nc, varname='v', timestep='all')    
    vel_magn = np.sqrt(data_fromnc_u**2 + data_fromnc_v**2)

    fig, ax = plt.subplots()
    ax.plot(data_fromnc_edgex, data_fromnc_edgey,'-b',linewidth=0.2)
    ax.plot(data_fromnc_edgex.T, data_fromnc_edgey.T,'-b',linewidth=0.2)
    plt.savefig(os.path.join(dir_output,'SFINCS_meshedge'))    

    fig, axs = plt.subplots(3,1, figsize=(14,9))
    for iT, timestep in enumerate([0,1,10]):
        ax=axs[iT]
        pc = ax.pcolor(data_fromnc_edgex, data_fromnc_edgey,vel_magn[timestep,:,:],cmap='jet')
        pc.set_clim([0,0.6])
        cbar = fig.colorbar(pc, ax=ax)
        cbar.set_label('velocity magnitude (%s)'%(data_fromnc_u.var_object.units))
        ax.set_title('t=%d (%s)'%(timestep, data_fromnc_u.var_times.loc[timestep]))
        ax.set_aspect('equal')
        thinning = 5
        ax.quiver(data_fromnc_edgex[::thinning,::thinning], data_fromnc_edgey[::thinning,::thinning], data_fromnc_u[timestep,::thinning,::thinning], data_fromnc_v[timestep,::thinning,::thinning], 
                  color='w')#,scale=3,width=0.005)#, edgecolor='face', cmap='jet')
    fig.tight_layout()
    plt.savefig(os.path.join(dir_output,'SFINCS_velocity_pcolorquiver'))


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
    


def test_delft3D_netcdf():
    dir_output = getmakeoutputdir(__file__,inspect.currentframe().f_code.co_name)
    """
    to get delft3D to write netCDF output instead of .dat files, add these lines to your mdf:
        FlNcdf= #maphis#
        ncFormat=4
    dir_output = './test_output'
    """
    
    import numpy as np
    import datetime as dt
    import matplotlib.pyplot as plt
    plt.close('all')
    
    from dfm_tools.get_nc import get_ncmodeldata#, get_netdata, plot_netmapdata
    from dfm_tools.get_nc_helpers import get_ncvardimlist, get_hisstationlist#, get_varname_fromnc
    from dfm_tools.regulargrid import uva2xymagdeg

    file_nc = r'p:\1220688-lake-kivu\3_modelling\1_FLOW\7_heatfluxinhis\063_netcdf\trim-thiery_002_coarse.nc'
    vars_pd, dims_pd = get_ncvardimlist(file_nc=file_nc)
    
    data_nc_XZ = get_ncmodeldata(file_nc=file_nc, varname='XZ')
    data_nc_YZ = get_ncmodeldata(file_nc=file_nc, varname='YZ')
    data_nc_ALFAS = get_ncmodeldata(file_nc=file_nc, varname='ALFAS') #contains rotation of all cells wrt real world
    data_nc_U1 = get_ncmodeldata(file_nc=file_nc, varname='U1',timestep='all',layer='all')
    data_nc_V1 = get_ncmodeldata(file_nc=file_nc, varname='V1',timestep='all',layer='all')
    #data_nc_S1 = get_ncmodeldata(file_nc=file_nc, varname='S1',timestep='all')
    data_nc_QNET = get_ncmodeldata(file_nc=file_nc, varname='QNET',timestep='all')
    #data_nc_QEVA = get_ncmodeldata(file_nc=file_nc, varname='QEVA',timestep='all')
    
    mask_XY = (data_nc_XZ==0) & (data_nc_YZ==0)
    mask_U = data_nc_U1==-999.
    mask_V = data_nc_V1==-999.
    #mask_U = (get_ncmodeldata(file_nc=file_nc, varname='KCU')==0)
    #mask_V = (get_ncmodeldata(file_nc=file_nc, varname='KCV')==0)
    data_nc_XZ[mask_XY] = np.nan
    data_nc_YZ[mask_XY] = np.nan
    data_nc_U1[mask_U] = np.nan
    data_nc_V1[mask_V] = np.nan
    #masking should work but quiver does not read masks for X and Y, so use own
    #data_nc_XZ.mask = mask_XY
    #data_nc_YZ.mask = mask_XY
    #data_nc_U1.mask = mask_U
    #data_nc_V1.mask = mask_V
    
    fig, ax = plt.subplots()
    ax.plot(data_nc_XZ,data_nc_YZ,'-b',linewidth=0.2)
    ax.plot(data_nc_XZ.T,data_nc_YZ.T,'-b',linewidth=0.2)
    ax.set_aspect('equal')
    plt.savefig(os.path.join(dir_output,'kivu_mesh'))
    
    fig, axs = plt.subplots(1,3, figsize=(16,7))
    for iT, timestep in enumerate([1,10,15]):
        ax=axs[iT]
        vel_x, vel_y, vel_magn, direction_naut_deg = uva2xymagdeg(u=data_nc_U1[timestep,90,:,:],v=data_nc_V1[timestep,90,:,:],alpha=data_nc_ALFAS)
        #pc = ax.pcolor(data_nc_XZ,data_nc_YZ,direction_naut_deg,cmap='jet')
        #pc.set_clim([0,360])
        pc = ax.pcolor(data_nc_XZ,data_nc_YZ,vel_magn,cmap='jet')
        pc.set_clim([0,0.15])
        cbar = fig.colorbar(pc, ax=ax)
        cbar.set_label('velocity magnitude (%s)'%(data_nc_U1.var_object.units))
        ax.set_title('t=%d (%s)'%(timestep, data_nc_U1.var_times.iloc[timestep]))
        ax.set_aspect('equal')
        ax.quiver(data_nc_XZ[::2,::2], data_nc_YZ[::2,::2], vel_x[::2,::2], vel_y[::2,::2], 
                  scale=3,color='w',width=0.005)#, edgecolor='face', cmap='jet')
    fig.tight_layout()
    plt.savefig(os.path.join(dir_output,'kivu_velocity_pcolor'))

    fig, axs = plt.subplots(1,3, figsize=(16,7))
    for iT, timestep in enumerate([1,10,15]):
        ax=axs[iT]
        vel_x, vel_y, vel_magn, direction_naut_deg = uva2xymagdeg(u=data_nc_U1[timestep,90,:,:],v=data_nc_V1[timestep,90,:,:],alpha=data_nc_ALFAS)
        #pc = ax.pcolor(data_nc_XZ,data_nc_YZ,vel_magn,cmap='jet')
        ax.set_title('t=%d (%s)'%(timestep, data_nc_U1.var_times.iloc[timestep]))
        ax.set_aspect('equal')
        pc = ax.quiver(data_nc_XZ[::2,::2], data_nc_YZ[::2,::2], vel_x[::2,::2], vel_y[::2,::2], vel_magn[::2,::2],
                  scale=3,color='w',width=0.005, edgecolor='face', cmap='jet')
        pc.set_clim([0,0.15])
        cbar = fig.colorbar(pc, ax=ax)
        cbar.set_label('velocity magnitude (%s)'%(data_nc_U1.var_object.units))
    fig.tight_layout()
    plt.savefig(os.path.join(dir_output,'kivu_velocity'))
    
    #QNET
    fig, axs = plt.subplots(1,3, figsize=(16,7))
    for iT, timestep in enumerate([1,10,15]):
        ax=axs[iT]
        pc = ax.pcolor(data_nc_XZ,data_nc_YZ,data_nc_QNET[iT,:,:],cmap='jet')
        pc.set_clim([-60,60])
        cbar = fig.colorbar(pc, ax=ax)
        cbar.set_label('%s (%s)'%(data_nc_QNET.var_varname, data_nc_QNET.var_object.units))
        ax.set_title('t=%d (%s)'%(timestep, data_nc_U1.var_times.iloc[timestep]))
        ax.set_aspect('equal')
    fig.tight_layout()
    plt.savefig(os.path.join(dir_output,'kivu_Qnet'))


    #FROM HIS data
    file_nc = r'p:\1220688-lake-kivu\3_modelling\1_FLOW\7_heatfluxinhis\063_netcdf\trih-thiery_002_coarse.nc'
    vars_pd, dims_pd = get_ncvardimlist(file_nc=file_nc)
    
    data_nc_NAMST = get_hisstationlist(file_nc=file_nc, varname_stat='NAMST')
    data_nc_ZWL = get_ncmodeldata(file_nc=file_nc, varname='ZWL',timestep='all',station='all')
    #data_nc_ZCURU = get_ncmodeldata(file_nc=file_nc, varname='ZCURU',timestep='all',station='all',layer='all')
    #data_nc_ZCURV = get_ncmodeldata(file_nc=file_nc, varname='ZCURV',timestep='all',station='all',layer='all')

    fig, ax = plt.subplots(figsize=(16,7))
    for iS in range(10):
        ax.plot(data_nc_ZWL.var_times,data_nc_ZWL[:,iS],label=data_nc_NAMST.iloc[iS], linewidth=1)
    ax.legend()
    ax.set_ylabel('%s (%s)'%(data_nc_ZWL.var_varname, data_nc_ZWL.var_object.units))
    ax.set_xlim([data_nc_ZWL.var_times[0],data_nc_ZWL.var_times[0]+dt.timedelta(days=14)])
    plt.savefig(os.path.join(dir_output,'kivu_his_ZWL'))




    #from MAP DATA CURVEDBEND
    file_nc = os.path.join(dir_testinput,'D3D_3D_sigma_curved_bend_nc\\trim-cb2-sal-added-3d.nc')
    vars_pd, dims_pd = get_ncvardimlist(file_nc=file_nc)
    
    data_nc_XZ = get_ncmodeldata(file_nc=file_nc, varname='XZ')
    data_nc_YZ = get_ncmodeldata(file_nc=file_nc, varname='YZ')
    data_nc_ALFAS = get_ncmodeldata(file_nc=file_nc, varname='ALFAS') #contains rotation of all cells wrt real world
    data_nc_U1 = get_ncmodeldata(file_nc=file_nc, varname='U1',timestep='all',layer='all')
    data_nc_V1 = get_ncmodeldata(file_nc=file_nc, varname='V1',timestep='all',layer='all')
    #data_nc_S1 = get_ncmodeldata(file_nc=file_nc, varname='S1',timestep='all')
    
    mask_XY = (data_nc_XZ==0) & (data_nc_YZ==0)
    #mask_U = data_nc_U1==-999.
    #mask_V = data_nc_V1==-999.
    #mask_U = (get_ncmodeldata(file_nc=file_nc, varname='KCU')==0)
    #mask_V = (get_ncmodeldata(file_nc=file_nc, varname='KCV')==0)
    data_nc_XZ[mask_XY] = np.nan
    data_nc_YZ[mask_XY] = np.nan
    #data_nc_U1[mask_U] = np.nan
    #data_nc_V1[mask_V] = np.nan
    #masking should work but quiver does not read masks for X and Y, so use own
    #data_nc_XZ.mask = mask_XY
    #data_nc_YZ.mask = mask_XY
    #data_nc_U1.mask = mask_U
    #data_nc_V1.mask = mask_V
    
    fig, ax = plt.subplots()
    ax.plot(data_nc_XZ,data_nc_YZ,'-b',linewidth=0.2)
    ax.plot(data_nc_XZ.T,data_nc_YZ.T,'-b',linewidth=0.2)
    ax.set_aspect('equal')
    plt.savefig(os.path.join(dir_output,'curvedbend_mesh'))
    
    fig, axs = plt.subplots(1,3, figsize=(16,5))
    for iT, timestep in enumerate([0,2,4]):
        ax=axs[iT]
        vel_x, vel_y, vel_magn, direction_naut_deg = uva2xymagdeg(u=data_nc_U1[timestep,9,:,:],v=data_nc_V1[timestep,9,:,:],alpha=data_nc_ALFAS)
        pc = ax.pcolor(data_nc_XZ,data_nc_YZ,vel_magn,cmap='jet')
        pc.set_clim([0,1.2])
        cbar = fig.colorbar(pc, ax=ax)
        cbar.set_label('velocity magnitude (%s)'%(data_nc_U1.var_object.units))
        ax.set_title('t=%d (%s)'%(timestep, data_nc_U1.var_times.iloc[timestep]))
        ax.set_aspect('equal')
        ax.quiver(data_nc_XZ[::2,::2], data_nc_YZ[::2,::2], vel_x[::2,::2], vel_y[::2,::2],
                  scale=8,color='w',width=0.005)#, edgecolor='face', cmap='jet')
    fig.tight_layout()
    plt.savefig(os.path.join(dir_output,'curvedbend_velocity_pcolor'))


    #FROM HIS data curvedbend
    file_nc = os.path.join(dir_testinput,'D3D_3D_sigma_curved_bend_nc\\trih-cb2-sal-added-3d.nc')
    vars_pd, dims_pd = get_ncvardimlist(file_nc=file_nc)
    
    data_nc_NAMST = get_hisstationlist(file_nc=file_nc, varname_stat='NAMST')
    data_nc_ZWL = get_ncmodeldata(file_nc=file_nc, varname='ZWL',timestep='all',station='all')
    #data_nc_ZCURU = get_ncmodeldata(file_nc=file_nc, varname='ZCURU',timestep='all',station='all',layer='all')
    #data_nc_ZCURV = get_ncmodeldata(file_nc=file_nc, varname='ZCURV',timestep='all',station='all',layer='all')

    fig, ax = plt.subplots(figsize=(16,7))
    for iS in range(5):
        ax.plot(data_nc_ZWL.var_times,data_nc_ZWL[:,iS],label=data_nc_NAMST.iloc[iS], linewidth=1)
    ax.legend()
    ax.set_ylabel('%s (%s)'%(data_nc_ZWL.var_varname, data_nc_ZWL.var_object.units))
    ax.set_xlim([data_nc_ZWL.var_times[0],data_nc_ZWL.var_times[0]+dt.timedelta(days=2)])
    plt.savefig(os.path.join(dir_output,'curvedbend_his_ZWL'))



    
 

#WARNING: this is excluded from the testbench, since Delft3D models that were converted with getdata.pl sometimes give corrupt variables (see comments in code for details). NEFIS conversion is not fully up to date in getdata.pl, whereas WAQUA conversion is
def EXCLUDE_test_delft3D_netcdf_convertedwith_getdata():
    dir_output = getmakeoutputdir(__file__,inspect.currentframe().f_code.co_name)
    """
    get the netcdf files via putty with:
        module load simona
        cd /p/1220688-lake-kivu/3_modelling/1_FLOW/7_heatfluxinhis/063
        getdata.pl -f trim-thiery_002_coarse.dat -v S1,U1,V1,ALFAS,QEVA -o netcdf
        http://simona.deltares.nl/release/doc/usedoc/getdata/getdata.pdf
        #double precision in trimfile causes this conversion to fail, use netCDF output in Delft3D instead
        
    get the netcdf files via putty with:
        module load simona
        cd ./D3D_3D_sigma_curved_bend
        getdata.pl -f trim-cb2-sal-added-3d.dat -v S1,U1,V1,ALFAS -o netcdf
        getdata.pl -f trih-cb2-sal-added-3d.dat -v ZWL,ZCURU,ZCURV,ALFAS,NAMST -o netcdf
        http://simona.deltares.nl/release/doc/usedoc/getdata/getdata.pdf
        #faulty data in NAMST variable, station names are not available

    dir_output = './test_output'
    """
    
    import numpy as np
    import datetime as dt
    import matplotlib.pyplot as plt
    plt.close('all')
    
    from dfm_tools.get_nc import get_ncmodeldata#, get_netdata, plot_netmapdata
    from dfm_tools.get_nc_helpers import get_ncvardimlist, get_hisstationlist#, get_varname_fromnc
    from dfm_tools.regulargrid import uva2xymagdeg
    
    #file_nc = r'p:\1220688-lake-kivu\3_modelling\1_FLOW\7_heatfluxinhis\063\trim-thiery_002_coarse.nc' #werkt niet
    file_nc = os.path.join(dir_testinput,'D3D_3D_sigma_curved_bend\\trim-cb2-sal-added-3d.nc')
    
    vars_pd, dims_pd = get_ncvardimlist(file_nc=file_nc)
    
    #set units attribute on time variable (workaround for bug)
    """
    from netCDF4 import Dataset
    data_nc = Dataset(file_nc,'r+')
    data_nc.variables['TIME'].units = 'minutes since 2012-01-01 00:00 +00:00'
    """
    data_nc_XZ = get_ncmodeldata(file_nc=file_nc, varname='map-const:XZ')
    data_nc_YZ = get_ncmodeldata(file_nc=file_nc, varname='map-const:YZ')
    data_nc_ALFAS = get_ncmodeldata(file_nc=file_nc, varname='map-const:ALFAS')
    #data_nc_TIME = get_ncmodeldata(file_nc=file_nc, varname='TIME',timestep='all')
    data_nc_U1 = get_ncmodeldata(file_nc=file_nc, varname='map-series:U1',timestep='all')
    data_nc_V1 = get_ncmodeldata(file_nc=file_nc, varname='map-series:V1',timestep='all')
    #data_nc_S1 = get_ncmodeldata(file_nc=file_nc, varname='map-series:S1',timestep='all')
    #data_nc_QNET = get_ncmodeldata(file_nc=file_nc, varname='map-series:QNET',timestep='all')
    
    mask_XY = (data_nc_XZ==0) & (data_nc_YZ==0)
    mask_U = data_nc_U1==-999.
    mask_V = data_nc_V1==-999.
    #mask_U = (get_ncmodeldata(file_nc=file_nc, varname='KCU')==0)
    #mask_V = (get_ncmodeldata(file_nc=file_nc, varname='KCV')==0)
    data_nc_XZ[mask_XY] = np.nan
    data_nc_YZ[mask_XY] = np.nan
    data_nc_U1[mask_U] = np.nan
    data_nc_V1[mask_V] = np.nan
    #masking should work but quiver does not read masks for X and Y, so use own
    #data_nc_XZ.mask = mask_XY
    #data_nc_YZ.mask = mask_XY
    #data_nc_U1.mask = mask_U
    #data_nc_V1.mask = mask_V
    
    fig, ax = plt.subplots()
    ax.plot(data_nc_XZ,data_nc_YZ,'-b',linewidth=0.2)
    ax.plot(data_nc_XZ.T,data_nc_YZ.T,'-b',linewidth=0.2)
    ax.set_aspect('equal')
    plt.savefig(os.path.join(dir_output,'curvedbend_mesh'))
    
    fig, axs = plt.subplots(1,4, figsize=(20,5))
    for iT, timestep in enumerate([0,1,2,4]):
        ax=axs[iT]
        vel_x, vel_y, vel_magn, direction_naut_deg = uva2xymagdeg(u=data_nc_U1[timestep,9,:,:],v=data_nc_V1[timestep,9,:,:],alpha=data_nc_ALFAS)
        #pc = ax.pcolor(data_nc_XZ,data_nc_YZ,direction_naut_deg,cmap='jet')
        #pc.set_clim([0,360])
        pc = ax.pcolor(data_nc_XZ,data_nc_YZ,vel_magn,cmap='jet')
        pc.set_clim([0,1.2])
        cbar = fig.colorbar(pc, ax=ax)
        cbar.set_label('velocity magnitude (%s)'%(data_nc_U1.var_object.units))
        ax.set_title('t=%d (%s)'%(timestep, data_nc_U1.var_times.iloc[timestep]))
        ax.set_aspect('equal')
        ax.quiver(data_nc_XZ[::2,::2], data_nc_YZ[::2,::2], vel_x[::2,::2], vel_y[::2,::2], 
                  scale=8,color='w',width=0.005)#, edgecolor='face', cmap='jet')
    fig.tight_layout()
    plt.savefig(os.path.join(dir_output,'curvedbend_velocity_pcolor'))

    fig, axs = plt.subplots(1,4, figsize=(20,5))
    for iT, timestep in enumerate([0,1,2,4]):
        ax=axs[iT]
        vel_x, vel_y, vel_magn, direction_naut_deg = uva2xymagdeg(u=data_nc_U1[timestep,9,:,:],v=data_nc_V1[timestep,9,:,:],alpha=data_nc_ALFAS)
        #pc = ax.pcolor(data_nc_XZ,data_nc_YZ,vel_magn,cmap='jet')
        ax.set_title('t=%d (%s)'%(timestep, data_nc_U1.var_times.iloc[timestep]))
        ax.set_aspect('equal')
        pc = ax.quiver(data_nc_XZ[::2,::2], data_nc_YZ[::2,::2], vel_x[::2,::2], vel_y[::2,::2], vel_magn[::2,::2],
                  scale=8,color='w',width=0.005, edgecolor='face', cmap='jet')
        pc.set_clim([0,1.2])
        cbar = fig.colorbar(pc, ax=ax)
        cbar.set_label('velocity magnitude (%s)'%(data_nc_U1.var_object.units))
    fig.tight_layout()
    plt.savefig(os.path.join(dir_output,'curvedbend_velocity'))
    

    #FROM HIS data
    file_nc = r'c:\DATA\dfm_tools_testdata\D3D_3D_sigma_curved_bend\trih-cb2-sal-added-3d.nc'
    vars_pd, dims_pd = get_ncvardimlist(file_nc=file_nc)
    
    data_nc_ZWL = get_ncmodeldata(file_nc=file_nc, varname='his-series:ZWL',timestep='all')#,station='all')
    #layers and stations are not yet taken care of properly (stations are incorrectly parsed so var/dim is commented in validvals list to avoid crash, 'Layers' dimension is not yet added to translation table)
    #data_nc_NAMST = get_ncmodeldata(file_nc=file_nc, varname='his-const:NAMST') #this should not work
    #data_nc_NAMST = get_hisstationlist(file_nc=file_nc, varname_stat='his-const:NAMST')
    #data_nc_ZCURU = get_ncmodeldata(file_nc=file_nc, varname='his-series:ZCURU',timestep='all')#,layer='all',station='all')
    #data_nc_ZCURV = get_ncmodeldata(file_nc=file_nc, varname='his-series:ZCURV',timestep='all')#,layer='all',station='all')

    fig, ax = plt.subplots(figsize=(16,7))
    for iS in range(6):
        ax.plot(data_nc_ZWL.var_times,data_nc_ZWL[:,iS], linewidth=1,label='unknown, broken NAMST variable')#,label=data_nc_NAMST.iloc[iS])
    ax.legend()
    ax.set_ylabel('%s (%s)'%(data_nc_ZWL.var_varname, data_nc_ZWL.var_object.units))
    ax.set_xlim([data_nc_ZWL.var_times[0],data_nc_ZWL.var_times[0]+dt.timedelta(days=2)])
    plt.savefig(os.path.join(dir_output,'curvedbend_his_ZWL'))






def test_waqua_netcdf_convertedwith_getdata():
    dir_output = getmakeoutputdir(__file__,inspect.currentframe().f_code.co_name)
    
    """
    get the netcdf files via putty with:
        module load simona
        cd /p/1204257-dcsmzuno/2019/DCSMv6/A01
        getdata.pl -f SDS-A01 -v SEP,VELU,VELV -o netcdf -d SDS-A01_map
        getdata.pl -f SDS-A01 -v ZWL,ZCURU,ZCURV,NAMWL,NAMC -o netcdf -d SDS-A01_his
        http://simona.deltares.nl/release/doc/usedoc/getdata/getdata.pdf
        
    get the netcdf files via putty with:
        module load simona
        cd /p/archivedprojects/1230049-zoutlastbeperking/Gaten_langsdam/Simulaties/OSR-model_GatenLangsdam/berekeningen/run7
        getdata.pl -f SDS-nsctri -v SEP,VELU,VELV -o netcdf -d SDS-nsctri_map
        getdata.pl -f SDS-nsctri -v ZWL,ZCURU,ZCURV,NAMWL,NAMC -o netcdf -d SDS-nsctri_his
        http://simona.deltares.nl/release/doc/usedoc/getdata/getdata.pdf

    get the netcdf files via putty with:
        module load simona
        cd /p/11205258-006-kpp2020_rmm-g6/C_Work/07_WAQUAresultaten/j15
        getdata.pl -f SDS-riv_tba -v SEP,VELU,VELV -o netcdf -d SDS-riv_tba_map
        getdata.pl -f SDS-riv_tba -v ZWL,ZCURU,ZCURV,NAMWL,NAMC -o netcdf -d SDS-riv_tba_his
        http://simona.deltares.nl/release/doc/usedoc/getdata/getdata.pdf
        
    dir_output = './test_output'
    """
    import datetime as dt
    import numpy as np
    import matplotlib.pyplot as plt
    plt.close('all')
    
    from dfm_tools.get_nc import get_ncmodeldata#, get_netdata, plot_netmapdata
    from dfm_tools.get_nc_helpers import get_ncvardimlist, get_hisstationlist#, get_varname_fromnc
    
    #MAP ZUNO
    file_nc = r'p:\1204257-dcsmzuno\2019\DCSMv6\A01\SDS-A01_map.nc'
    vars_pd, dims_pd = get_ncvardimlist(file_nc=file_nc)
    
    data_nc_x = get_ncmodeldata(file_nc=file_nc, varname='grid_x')
    data_nc_y = get_ncmodeldata(file_nc=file_nc, varname='grid_y')
    data_nc_SEP = get_ncmodeldata(file_nc=file_nc, varname='SEP',timestep='all')
    #data_nc_VELU = get_ncmodeldata(file_nc=file_nc, varname='VELU',timestep='all')
    #data_nc_VELV = get_ncmodeldata(file_nc=file_nc, varname='VELV',timestep='all')
    data_nc_xcen = np.mean(data_nc_x, axis=2)
    data_nc_ycen = np.mean(data_nc_y, axis=2)
    
    fig, ax = plt.subplots()
    #vel_x, vel_y, vel_magn, direction_naut_deg = uva2xymagdeg(u=data_nc_U1[timestep,90,:,:],v=data_nc_V1[timestep,90,:,:],alpha=data_nc_ALFAS)
    #pc = ax.pcolor(data_nc_XZ,data_nc_YZ,direction_naut_deg,cmap='jet')
    #pc.set_clim([0,360])
    timestep=0
    pc = ax.pcolor(data_nc_xcen,data_nc_ycen,data_nc_SEP[timestep,:,:],cmap='jet')
    pc.set_clim([-0.1,0.1])
    cbar = fig.colorbar(pc, ax=ax)
    cbar.set_label('%s (%s)'%(data_nc_SEP.var_varname, data_nc_SEP.var_object.units))
    ax.set_title('t=%d (%s)'%(timestep, data_nc_SEP.var_times.iloc[timestep]))
    ax.set_aspect('equal')
    #ax.quiver(data_nc_XZ[::2,::2], data_nc_YZ[::2,::2], vel_x[::2,::2], vel_y[::2,::2], 
    #          scale=3,color='w',width=0.005)#, edgecolor='face', cmap='jet')
    fig.tight_layout()
    plt.savefig(os.path.join(dir_output,'waqua_DCSM_map_wl'))

    
    #HIS ZUNO
    file_nc = r'p:\1204257-dcsmzuno\2019\DCSMv6\A01\SDS-A01_his.nc'
    vars_pd, dims_pd = get_ncvardimlist(file_nc=file_nc)
    data_nc_NAMWL = get_hisstationlist(file_nc=file_nc, varname_stat='NAMWL')
    #data_nc_NAMC = get_hisstationlist(file_nc=file_nc, varname_stat='NAMC')
    data_nc_ZWL = get_ncmodeldata(file_nc=file_nc, varname='ZWL',timestep='all',station='all')
    #data_nc_ZCURU = get_ncmodeldata(file_nc=file_nc, varname='ZCURU',timestep='all',station='all')
    #data_nc_ZCURV = get_ncmodeldata(file_nc=file_nc, varname='ZCURV',timestep='all',station='all')

    fig, ax = plt.subplots(figsize=(16,7))
    for iS in range(10):
        ax.plot(data_nc_ZWL.var_times,data_nc_ZWL[:,iS],label=data_nc_NAMWL.iloc[iS], linewidth=1)
    ax.legend()
    ax.set_ylabel('%s (%s)'%(data_nc_ZWL.var_varname, data_nc_ZWL.var_object.units))
    ax.set_xlim([data_nc_ZWL.var_times[0],data_nc_ZWL.var_times[0]+dt.timedelta(days=14)])
    plt.savefig(os.path.join(dir_output,'waqua_DSCM_his_ZWL'))
    
    
    
    #MAP OSR
    file_nc = r'p:\11205258-006-kpp2020_rmm-g6\C_Work\ZZ_Jelmer\SDS-nsctri_map.nc'
    vars_pd, dims_pd = get_ncvardimlist(file_nc=file_nc)
    
    data_nc_x = get_ncmodeldata(file_nc=file_nc, varname='grid_x')
    data_nc_y = get_ncmodeldata(file_nc=file_nc, varname='grid_y')
    data_nc_xcen = np.mean(data_nc_x, axis=2)
    data_nc_ycen = np.mean(data_nc_y, axis=2)
    
    timestep = 10
    data_nc_SEP = get_ncmodeldata(file_nc=file_nc, varname='SEP',timestep=timestep)
    data_nc_VELU = get_ncmodeldata(file_nc=file_nc, varname='VELU',timestep=timestep, layer=9)
    data_nc_VELV = get_ncmodeldata(file_nc=file_nc, varname='VELV',timestep=timestep, layer=9)
    
    fig, ax = plt.subplots()
    pc = ax.pcolor(data_nc_xcen,data_nc_ycen,data_nc_SEP[0,:,:],cmap='jet')
    pc.set_clim([0,2])
    cbar = fig.colorbar(pc, ax=ax)
    cbar.set_label('%s (%s)'%(data_nc_SEP.var_varname, data_nc_SEP.var_object.units))
    ax.set_title('t=%d (%s)'%(timestep, data_nc_SEP.var_times.loc[timestep]))
    ax.set_aspect('equal')
    fig.tight_layout()
    plt.savefig(os.path.join(dir_output,'waqua_OSR_map_wl'))

    fig, ax = plt.subplots()
    vel_magn = np.sqrt(data_nc_VELU**2 + data_nc_VELV**2)
    pc = ax.pcolor(data_nc_xcen,data_nc_ycen,vel_magn[0,:,:,0],cmap='jet')
    pc.set_clim([0,1])
    cbar = fig.colorbar(pc, ax=ax)
    cbar.set_label('velocity magnitude (%s)'%(data_nc_VELU.var_object.units))
    ax.set_title('t=%d (%s)'%(timestep, data_nc_VELU.var_times.loc[timestep]))
    ax.set_aspect('equal')
    thinning = 10
    ax.quiver(data_nc_xcen[::thinning,::thinning], data_nc_ycen[::thinning,::thinning], data_nc_VELU[0,::thinning,::thinning,0], data_nc_VELV[0,::thinning,::thinning,0], 
              color='w',scale=10)#,width=0.005)#, edgecolor='face', cmap='jet')
    ax.set_xlim([58000, 66000])
    ax.set_ylim([442000, 448000])
    fig.tight_layout()
    plt.savefig(os.path.join(dir_output,'waqua_OSR_map_vel'))
    
    #HIS OSR
    file_nc = r'p:\11205258-006-kpp2020_rmm-g6\C_Work\ZZ_Jelmer\SDS-nsctri_his.nc'
    vars_pd, dims_pd = get_ncvardimlist(file_nc=file_nc)
    data_nc_NAMWL = get_hisstationlist(file_nc=file_nc, varname_stat='NAMWL')
    #data_nc_NAMC = get_hisstationlist(file_nc=file_nc, varname_stat='NAMC')
    data_nc_ZWL = get_ncmodeldata(file_nc=file_nc, varname='ZWL',timestep='all',station='all')
    #data_nc_ZCURU = get_ncmodeldata(file_nc=file_nc, varname='ZCURU',timestep='all',station='all',layer='all')
    #data_nc_ZCURV = get_ncmodeldata(file_nc=file_nc, varname='ZCURV',timestep='all',station='all',layer='all')

    fig, ax = plt.subplots(figsize=(16,7))
    for iS in range(10):
        ax.plot(data_nc_ZWL.var_times,data_nc_ZWL[:,iS],label=data_nc_NAMWL.iloc[iS], linewidth=1)
    ax.legend()
    ax.set_ylabel('%s (%s)'%(data_nc_ZWL.var_varname, data_nc_ZWL.var_object.units))
    ax.set_xlim([data_nc_ZWL.var_times[0],data_nc_ZWL.var_times[0]+dt.timedelta(days=14)])
    plt.savefig(os.path.join(dir_output,'waqua_OSR_his_ZWL'))
    
    
    
    #MAP RMM
    file_nc = r'p:\11205258-006-kpp2020_rmm-g6\C_Work\07_WAQUAresultaten\j15\SDS-riv_tba_map.nc'
    vars_pd, dims_pd = get_ncvardimlist(file_nc=file_nc)
    
    data_nc_x = get_ncmodeldata(file_nc=file_nc, varname='grid_x')
    data_nc_y = get_ncmodeldata(file_nc=file_nc, varname='grid_y')
    data_nc_xcen = np.mean(data_nc_x, axis=2)
    data_nc_ycen = np.mean(data_nc_y, axis=2)
    
    timestep = 1
    data_nc_SEP = get_ncmodeldata(file_nc=file_nc, varname='SEP',timestep=timestep)
    data_nc_VELU = get_ncmodeldata(file_nc=file_nc, varname='VELU',timestep=timestep, layer=0)
    data_nc_VELV = get_ncmodeldata(file_nc=file_nc, varname='VELV',timestep=timestep, layer=0)
    
    fig, ax = plt.subplots(figsize=(16,7))
    pc = ax.pcolor(data_nc_xcen,data_nc_ycen,data_nc_SEP[0,:,:],cmap='jet')
    pc.set_clim([0,3])
    cbar = fig.colorbar(pc, ax=ax)
    cbar.set_label('%s (%s)'%(data_nc_SEP.var_varname, data_nc_SEP.var_object.units))
    ax.set_title('t=%d (%s)'%(timestep, data_nc_SEP.var_times.loc[timestep]))
    ax.set_aspect('equal')
    fig.tight_layout()
    plt.savefig(os.path.join(dir_output,'waqua_RMM_map_wl'))

    fig, ax = plt.subplots()
    vel_magn = np.sqrt(data_nc_VELU**2 + data_nc_VELV**2)
    pc = ax.pcolor(data_nc_xcen,data_nc_ycen,vel_magn[0,:,:,0],cmap='jet')
    pc.set_clim([0,1])
    cbar = fig.colorbar(pc, ax=ax)
    cbar.set_label('velocity magnitude (%s)'%(data_nc_VELU.var_object.units))
    ax.set_title('t=%d (%s)'%(timestep, data_nc_VELU.var_times.loc[timestep]))
    ax.set_aspect('equal')
    thinning = 10
    ax.quiver(data_nc_xcen[::thinning,::thinning], data_nc_ycen[::thinning,::thinning], data_nc_VELU[0,::thinning,::thinning,0], data_nc_VELV[0,::thinning,::thinning,0], 
              color='w',scale=10)#,width=0.005)#, edgecolor='face', cmap='jet')
    ax.set_xlim([61000, 72000])
    ax.set_ylim([438000, 446000])
    fig.tight_layout()
    plt.savefig(os.path.join(dir_output,'waqua_RMM_map_vel'))
    
    
    #HIS RMM
    file_nc = r'p:\11205258-006-kpp2020_rmm-g6\C_Work\07_WAQUAresultaten\j15\SDS-riv_tba_his.nc'
    vars_pd, dims_pd = get_ncvardimlist(file_nc=file_nc)
    data_nc_NAMWL = get_hisstationlist(file_nc=file_nc, varname_stat='NAMWL')
    #data_nc_NAMC = get_hisstationlist(file_nc=file_nc, varname_stat='NAMC')
    data_nc_ZWL = get_ncmodeldata(file_nc=file_nc, varname='ZWL',timestep='all',station='all')
    #data_nc_ZCURU = get_ncmodeldata(file_nc=file_nc, varname='ZCURU',timestep='all',station='all',layer='all')
    #data_nc_ZCURV = get_ncmodeldata(file_nc=file_nc, varname='ZCURV',timestep='all',station='all',layer='all')

    fig, ax = plt.subplots(figsize=(16,7))
    for iS in range(10):
        ax.plot(data_nc_ZWL.var_times,data_nc_ZWL[:,iS],label=data_nc_NAMWL.iloc[iS], linewidth=1)
    ax.legend()
    ax.set_ylabel('%s (%s)'%(data_nc_ZWL.var_varname, data_nc_ZWL.var_object.units))
    ax.set_xlim([data_nc_ZWL.var_times[0],data_nc_ZWL.var_times[0]+dt.timedelta(days=14)])
    plt.savefig(os.path.join(dir_output,'waqua_RMM_his_ZWL'))

    