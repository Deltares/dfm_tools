# -*- coding: utf-8 -*-
"""
Created on Thu Feb 13 23:10:51 2020

@author: veenstra
"""

import pytest
import inspect
import os

dir_testinput = os.path.join(r'c:/DATA','dfm_tools_testdata')
from tests.TestUtils import getmakeoutputdir


@pytest.mark.unittest    
def test_getvarnamemapnc():
    """
    this test tests if a netcdf varname can be retrieved from the 'dictionary' and if the variable can be retrieved from de netcdf
    """
    
    from netCDF4 import Dataset
    
    from dfm_tools.get_nc_helpers import get_varname_mapnc
    
    file_map = os.path.join(dir_testinput,r'DFM_3D_z_Grevelingen\computations\run01\DFM_OUTPUT_Grevelingen-FM\Grevelingen-FM_0000_map.nc')
    #file_net = os.path.join(dir_testinput,r'DFM_3D_z_Grevelingen\computations\run01\Grevelingen_FM_grid_20190603_net.nc')
    
    file_nc = file_map
    data_nc = Dataset(file_nc)
    varname_requested = 'NetNode_y' #is actually in file, so this is not a good test
    
    varname = get_varname_mapnc(data_nc,varname_requested)
    data_nc_var = data_nc.variables[varname]
    dimname = data_nc_var.dimensions[0]
    
    # Vefiry expectations
    assert varname == 'mesh2d_node_y'
    assert dimname == 'nmesh2d_node'





@pytest.mark.unittest    
def test_getncmatchingvarlist():
    """
    this test tests retrieves a pandas list of variable long names in a netcdf that match the pattern, useful for waq variables.
    """
    from dfm_tools.get_nc_helpers import get_ncvardimlist
    
    file_nc = os.path.join(dir_testinput,r'DFM_3D_z_Grevelingen\computations\run01\DFM_OUTPUT_Grevelingen-FM\Grevelingen-FM_0000_map.nc')
    vars_pd, dims_pd = get_ncvardimlist(file_nc=file_nc)

    pattern = 'Flow .*component'
    vars_pd_matching = vars_pd[vars_pd.loc[:,'nc_varlongnames'].str.match(pattern)] #does not have to stop after pattern
    #vars_pd_matching = vars_pd[vars_pd.loc[:,'nc_varlongnames'].str.startswith('Flow') & vars_pd.loc[:,'nc_varlongnames'].str.endswith('component')]
    varkeys_list_matching = list(vars_pd_matching['nc_varkeys'])
    
    assert varkeys_list_matching == ['mesh2d_ucx', 'mesh2d_ucy', 'mesh2d_ucz', 'mesh2d_ucxa', 'mesh2d_ucya']




@pytest.mark.parametrize("file_nc, expected_size", [pytest.param(os.path.join(dir_testinput,r'DFM_3D_z_Grevelingen\computations\run01\DFM_OUTPUT_Grevelingen-FM\Grevelingen-FM_0000_map.nc'), 5599, id='from 1 map partion Grevelingen'),
                                                    #pytest.param(r'p:\11205258-006-kpp2020_rmm-g6\C_Work\01_Rooster\final_totaalmodel\rooster_rmm_v1p5_net.nc', 44804?, id='fromnet RMM'),
                                                    pytest.param(os.path.join(dir_testinput,r'DFM_3D_z_Grevelingen\computations\run01\Grevelingen_FM_grid_20190603_net.nc'), 44804, id='fromnet Grevelingen')])
@pytest.mark.unittest
def test_UGrid(file_nc, expected_size):
    from dfm_tools.ugrid import UGrid
    
    ugrid = UGrid.fromfile(file_nc)
    
    assert ugrid.verts.shape[0] == expected_size


@pytest.mark.parametrize("file_nc, expected_size", [pytest.param(os.path.join(dir_testinput,r'DFM_3D_z_Grevelingen\computations\run01\DFM_OUTPUT_Grevelingen-FM\Grevelingen-FM_0000_map.nc'), 44796, id='from partitioned map Grevelingen'),
                                                    #pytest.param(r'p:\11205258-006-kpp2020_rmm-g6\C_Work\01_Rooster\final_totaalmodel\rooster_rmm_v1p5_net.nc', 44804?, id='fromnet RMM'),
                                                    pytest.param(os.path.join(dir_testinput,r'DFM_3D_z_Grevelingen\computations\run01\Grevelingen_FM_grid_20190603_net.nc'), 44804, id='fromnet Grevelingen')])
@pytest.mark.unittest
def test_getnetdata(file_nc, expected_size):
    from dfm_tools.get_nc import get_netdata
    
    ugrid = get_netdata(file_nc)
    
    assert ugrid.verts.shape[0] == expected_size


@pytest.mark.unittest
def test_getncmodeldata_timeid():
    from dfm_tools.get_nc import get_ncmodeldata
    
    file_map1 = os.path.join(dir_testinput,r'DFM_sigma_curved_bend\DFM_OUTPUT_cb_3d\cb_3d_map.nc')
    data_frommap = get_ncmodeldata(file_nc=file_map1, varname='mesh2d_sa1', timestep=1, layer=5)#, multipart=False)
    
    assert (data_frommap.data[0,0,0] - 31. ) < 1E-9
    



@pytest.mark.unittest
def test_getncmodeldata_datetime():
    import numpy as np
    import datetime as dt
    
    from dfm_tools.get_nc import get_ncmodeldata

    file_map1 = os.path.join(dir_testinput,r'DFM_sigma_curved_bend\DFM_OUTPUT_cb_3d\cb_3d_map.nc')
    data_frommap = get_ncmodeldata(file_nc=file_map1, varname='mesh2d_sa1', timestep=np.arange(dt.datetime(2001,1,1),dt.datetime(2001,1,2),dt.timedelta(hours=1)), layer=5)#, multipart=False)
    
    assert (data_frommap.data[0,0,0] - 31. ) < 1E-9
    


@pytest.mark.systemtest
def test_getplotfoudata():
    
    dir_output = getmakeoutputdir(__file__,inspect.currentframe().f_code.co_name)
    #dir_output = '.'

    import matplotlib.pyplot as plt
    plt.close('all')
    
    from dfm_tools.get_nc import get_netdata, get_ncmodeldata, plot_netmapdata
    
    file_nc = os.path.join(dir_testinput,r'DFM_fou_RMM\RMM_dflowfm_0000_fou.nc')

    ugrid = get_netdata(file_nc=file_nc)
    data_fromfou = get_ncmodeldata(file_nc=file_nc, varname='mesh2d_fourier003_mean')#, multipart=False)
    
    fig, ax = plt.subplots()
    pc = plot_netmapdata(ugrid.verts, values=data_fromfou, ax=None, linewidth=0.5, color="crimson", facecolor="None")
    pc.set_clim([0,10])
    fig.colorbar(pc)
    ax.set_aspect('equal')
    plt.savefig(os.path.join(dir_output,os.path.basename(file_nc).replace('.nc','')))

    assert ugrid.verts.shape[0] == data_fromfou.shape[0]




    
@pytest.mark.parametrize("file_nc", [pytest.param(os.path.join(dir_testinput,r'DFM_3D_z_Grevelingen\computations\run01\Grevelingen_FM_grid_20190603_net.nc'), id='Grevelingen'),
                                     pytest.param(os.path.join(dir_testinput,'vanNithin','myortho3_RGFGRID_net.nc'), id='Nithin'),
                                     pytest.param(r'p:\11205258-006-kpp2020_rmm-g6\C_Work\01_Rooster\final_totaalmodel\rooster_rmm_v1p5_net.nc', id='RMM')])
@pytest.mark.acceptance
def test_getnetdata_plotnet(file_nc):
    """
    this test retrieves grid data and plots it
    
    file_nc = os.path.join(dir_testinput,'DFM_3D_z_Grevelingen','computations','run01','Grevelingen_FM_grid_20190603_net.nc')
    file_nc = 'p:\\11205258-006-kpp2020_rmm-g6\\C_Work\\01_Rooster\\final_totaalmodel\\rooster_rmm_v1p5_net.nc'
    """
    
    dir_output = getmakeoutputdir(__file__,inspect.currentframe().f_code.co_name)
    #dir_output = '.'

    import matplotlib.pyplot as plt
    plt.close('all')

    from dfm_tools.get_nc import get_netdata, plot_netmapdata

    print('plot only grid from net.nc')
    ugrid = get_netdata(file_nc=file_nc)
    fig, ax = plt.subplots()
    plot_netmapdata(ugrid.verts, values=None, ax=None, linewidth=0.5, color="crimson", facecolor="None")
    ax.set_aspect('equal')
    plt.savefig(os.path.join(dir_output,os.path.basename(file_nc).replace('.nc','')))
    

    
    
    
@pytest.mark.acceptance
def test_getsobekmodeldata():
    """
    this test retrieves sobek observation data and plots it
    """
    dir_output = getmakeoutputdir(__file__,inspect.currentframe().f_code.co_name)
    #dir_output = '.'

    import matplotlib.pyplot as plt
    plt.close('all')

    from dfm_tools.get_nc import get_ncmodeldata
    #from dfm_tools.get_nc_helpers import get_hisstationlist
    
    file_nc = os.path.join(dir_testinput,'KenmerkendeWaarden','observations.nc')
    
    #station_names = get_hisstationlist(file_nc=file_nc, varname_stat="observation_id")
    data_fromsobek = get_ncmodeldata(file_nc=file_nc, varname='water_level', station=['Maasmond','HKVHLD','MO_1035.00'], timestep='all')
    
    fig, ax = plt.subplots()
    for iL in range(data_fromsobek.shape[1]):
        ax.plot(data_fromsobek.var_times,data_fromsobek[:,iL],'-', label=data_fromsobek.var_stations.iloc[iL])
    ax.legend()
    plt.savefig(os.path.join(dir_output,'%s_waterlevel'%(os.path.basename(file_nc).replace('.nc',''))))
    
    
    
    
    
    
@pytest.mark.parametrize("file_nc", [pytest.param(os.path.join(dir_testinput,'vanNithin','tttz_0000_his.nc'), id='Nithin'),
                                     pytest.param(os.path.join(dir_testinput,'DFM_3D_z_Grevelingen\\computations\\run01\\DFM_OUTPUT_Grevelingen-FM\\Grevelingen-FM_0000_his.nc'), id='Grevelingen')])
@pytest.mark.acceptance
def test_gethismodeldata(file_nc):
    """
    this test retrieves his data and plots it
    file_nc = os.path.join(dir_testinput,'DFM_3D_z_Grevelingen\\computations\\run01\\DFM_OUTPUT_Grevelingen-FM\\Grevelingen-FM_0000_his.nc')
    """
    dir_output = getmakeoutputdir(__file__,inspect.currentframe().f_code.co_name)
    #dir_output = '.'

    import numpy as np
    import matplotlib.pyplot as plt
    plt.close('all')
    
    from dfm_tools.get_nc import get_ncmodeldata
    
    def cen2cor(time_cen):
        #convert time centers2corners (more accurate representation in zt-plot, but can also be skipped)
        import pandas as pd
        time_int = (time_cen.iloc[2]-time_cen.iloc[1])
        time_cor = data_fromhis_temp.var_times-time_int/2
        time_cor = time_cor.append(pd.Series([time_cor.iloc[-1]+time_int]))
        return time_cor

    if 'Grevelingen-FM_0000' in file_nc:
        #file_nc = os.path.join(dir_testinput,r'DFM_3D_z_Grevelingen\computations\run01\DFM_OUTPUT_Grevelingen-FM\Grevelingen-FM_0000_his.nc')
        station = 'all'
        station_zt = 'GTSO-02'
    elif 'tttz' in file_nc: #NITHIN
        #file_nc = os.path.join(dir_testinput,'vanNithin','tttz_0000_his.nc')
        station = ['Peiraias', 'Ovrios_2','Ovrios','Ovrios']
        station_zt = 'Ortholithi'
    
    print('plot bedlevel from his')
    data_fromhis = get_ncmodeldata(file_nc=file_nc, varname='bedlevel', station=station)#, multipart=False)
    fig, ax = plt.subplots()
    ax.plot(data_fromhis.var_stations,data_fromhis,'-')
    ax.tick_params('x',rotation=90)
    plt.savefig(os.path.join(dir_output,'%s_bedlevel'%(os.path.basename(file_nc).replace('.nc',''))))

    print('plot waterlevel from his')
    data_fromhis = get_ncmodeldata(file_nc=file_nc, varname='waterlevel', timestep='all', station=station)#, multipart=False)
    fig, ax = plt.subplots()
    ax.plot(data_fromhis.var_times,data_fromhis,'-')
    plt.savefig(os.path.join(dir_output,'%s_waterlevel'%(os.path.basename(file_nc).replace('.nc',''))))
    
    print('plot salinity from his')
    data_fromhis = get_ncmodeldata(file_nc=file_nc, varname='salinity', timestep='all', layer=5, station=station)#, multipart=False)
    data_fromhis_flat = data_fromhis[:,:,0]
    fig, ax = plt.subplots()
    ax.plot(data_fromhis.var_times,data_fromhis_flat,'-')
    plt.savefig(os.path.join(dir_output,'%s_salinity'%(os.path.basename(file_nc).replace('.nc',''))))

    print('plot salinity over depth')
    #depth retrieval is probably wrong
    data_fromhis_depth = get_ncmodeldata(file_nc=file_nc, varname='zcoordinate_c', timestep=4, layer='all', station=station)
    data_fromhis = get_ncmodeldata(file_nc=file_nc, varname='salinity', timestep=4, layer='all', station=station)
    fig, ax = plt.subplots()
    ax.plot(data_fromhis[0,:,:].T, data_fromhis_depth[0,:,:].T,'-')
    ax.legend(data_fromhis.var_stations)
    plt.savefig(os.path.join(dir_output,'%s_salinityoverdepth'%(os.path.basename(file_nc).replace('.nc',''))))

    print('zt temperature plot')
    #WARNING: layers in dfowfm hisfile are currently incorrect, check your figures carefully
    data_fromhis_zcen = get_ncmodeldata(file_nc=file_nc, varname='zcoordinate_c', timestep=range(40,100), layer= 'all', station=station_zt)
    data_fromhis_zcor = get_ncmodeldata(file_nc=file_nc, varname='zcoordinate_w', timestep=range(40,100), station=station_zt)
    data_fromhis_temp = get_ncmodeldata(file_nc=file_nc, varname='temperature', timestep=range(40,100), layer= 'all', station=station_zt)
    #time_cor = cen2cor(data_fromhis_temp.var_times)
    time_cen = data_fromhis_temp.var_times
    time_int = (time_cen.iloc[2]-time_cen.iloc[1])
    time_cor = data_fromhis_temp.var_times-time_int/2
    # generate 2 2d grids for the x & y bounds (you can also give one 2D array as input in case of eg time varying z coordinates)
    time_mesh_cor = np.tile(time_cor,(data_fromhis_zcor.shape[-1],1)).T
    time_mesh_cen = np.tile(data_fromhis_temp.var_times,(data_fromhis_zcen.shape[-1],1)).T
    fig, ax = plt.subplots(figsize=(12,5))
    c = ax.pcolormesh(time_mesh_cor, data_fromhis_zcor[:,0,:], data_fromhis_temp[:,0,:],cmap='jet')
    fig.colorbar(c)
    #contour
    CS = ax.contour(time_mesh_cen,data_fromhis_zcen[:,0,:],data_fromhis_temp[:,0,:],6, colors='k', linewidths=0.8, linestyles='solid')
    ax.clabel(CS, fontsize=10)
    plt.savefig(os.path.join(dir_output,'%s_zt_temp'%(os.path.basename(file_nc).replace('.nc',''))))

    

    
    

@pytest.mark.parametrize("file_nc", [pytest.param(os.path.join(dir_testinput,r'DFM_sigma_curved_bend\DFM_OUTPUT_cb_3d\cb_3d_map.nc'), id='curvibend'),
                                     pytest.param(os.path.join(dir_testinput,r'DFM_3D_z_Grevelingen\computations\run01\DFM_OUTPUT_Grevelingen-FM\Grevelingen-FM_0000_map.nc'), id='Grevelingen'),
                                     pytest.param(r'p:\11205258-006-kpp2020_rmm-g6\C_Work\08_RMM_FMmodel\computations\run_156\DFM_OUTPUT_RMM_dflowfm\RMM_dflowfm_0000_map.nc', id='RMM')])
@pytest.mark.acceptance
def test_getnetdata_getmapmodeldata_plotnetmapdata(file_nc):
    """
    this test retrieves grid data, retrieves map data, and plots it
    file_nc = os.path.join(dir_testinput,'DFM_3D_z_Grevelingen','computations','run01','DFM_OUTPUT_Grevelingen-FM','Grevelingen-FM_0000_map.nc')
    """
    dir_output = getmakeoutputdir(__file__,inspect.currentframe().f_code.co_name)
    #dir_output = '.'

    import matplotlib.pyplot as plt
    plt.close('all')
    
    from dfm_tools.get_nc import get_netdata, get_ncmodeldata, plot_netmapdata
    
    if 'cb_3d_map' in file_nc:
        timestep = 3
        layer = 5
        clim_bl = None
        clim_wl = [-0.5,1]
        clim_sal = None
        clim_tem = None
    elif 'Grevelingen-FM_0000_map' in file_nc:
        timestep = 3
        layer = 33
        clim_bl = None
        clim_wl = [-0.5,1]
        clim_sal = [28,30.2]
        clim_tem = [4,10]
    elif 'RMM_dflowfm_0000_map' in file_nc:
        timestep = 50
        layer = None
        clim_bl = [-10,10]
        clim_wl = [-2,2]
        clim_sal = None
        clim_tem = None
    else:
        raise Exception('ERROR: no settings provided for this mapfile')
        
    
    #PLOT GRID
    print('plot only grid from mapdata')
    ugrid_all = get_netdata(file_nc=file_nc)#,multipart=False)
    fig, ax = plt.subplots()
    pc = plot_netmapdata(ugrid_all.verts, values=None, ax=None, linewidth=0.5, color="crimson", facecolor="None")
    ax.set_aspect('equal')
    plt.savefig(os.path.join(dir_output,'%s_grid'%(os.path.basename(file_nc).replace('.nc',''))))


    #PLOT bedlevel
    if not 'cb_3d_map' in file_nc:
        print('plot grid and values from mapdata (constantvalue, 1 dim)')
        ugrid = get_netdata(file_nc=file_nc)#,multipart=False)
        #iT = 3 #for iT in range(10):
        data_frommap = get_ncmodeldata(file_nc=file_nc, varname='mesh2d_flowelem_bl')#, multipart=False)
        data_frommap_flat = data_frommap.flatten()
        fig, ax = plt.subplots()
        pc = plot_netmapdata(ugrid.verts, values=data_frommap_flat, ax=None, linewidth=0.5, cmap="jet")
        pc.set_clim(clim_bl)
        fig.colorbar(pc, ax=ax)
        ax.set_aspect('equal')
        plt.savefig(os.path.join(dir_output,'%s_mesh2d_flowelem_bl'%(os.path.basename(file_nc).replace('.nc',''))))
        
    
    #PLOT water level on map
    print('plot grid and values from mapdata (waterlevel, 2dim)')
    data_frommap = get_ncmodeldata(file_nc=file_nc, varname='mesh2d_s1', timestep=timestep)#, multipart=False)
    data_frommap_flat = data_frommap.flatten()
    fig, ax = plt.subplots()
    pc = plot_netmapdata(ugrid_all.verts, values=data_frommap_flat, ax=None, linewidth=0.5, cmap="jet")
    pc.set_clim(clim_wl)
    fig.colorbar(pc, ax=ax)
    ax.set_aspect('equal')
    plt.savefig(os.path.join(dir_output,'%s_mesh2d_s1'%(os.path.basename(file_nc).replace('.nc',''))))

    #PLOT var layer on map
    if not 'RMM_dflowfm_0000_map' in file_nc:
        print('plot grid and values from mapdata (salinity on layer, 3dim)')
        data_frommap = get_ncmodeldata(file_nc=file_nc, varname='mesh2d_sa1', timestep=timestep, layer=layer)#, multipart=False)
        data_frommap_flat = data_frommap.flatten()
        fig, ax = plt.subplots()
        pc = plot_netmapdata(ugrid_all.verts, values=data_frommap_flat, ax=None, linewidth=0.5, cmap="jet")
        pc.set_clim(clim_sal)
        fig.colorbar(pc, ax=ax)
        ax.set_aspect('equal')
        plt.savefig(os.path.join(dir_output,'%s_mesh2d_sa1'%(os.path.basename(file_nc).replace('.nc',''))))
    
        print('plot grid and values from mapdata (temperature on layer, 3dim)')
        data_frommap = get_ncmodeldata(file_nc=file_nc, varname='mesh2d_tem1', timestep=timestep, layer=layer)#, multipart=False)
        data_frommap_flat = data_frommap.flatten()
        fig, ax = plt.subplots()
        pc = plot_netmapdata(ugrid_all.verts, values=data_frommap_flat, ax=None, linewidth=0.5, cmap="jet")
        pc.set_clim(clim_tem)
        fig.colorbar(pc, ax=ax)
        ax.set_aspect('equal')
        plt.savefig(os.path.join(dir_output,'%s_mesh2d_tem1'%(os.path.basename(file_nc).replace('.nc',''))))





@pytest.mark.parametrize("file_nc", [pytest.param('p:\\11201806-sophie\\Oosterschelde\\WAQ\\r03\\postprocessing\\oost_tracer_map_backup.nc', id='oost_tracer_map'),
                                     pytest.param('p:\\11201806-sophie\\Oosterschelde\\WAQ\\r02\\postprocessing\\oost_tracer_2_map.nc', id='oost_tracer_2_map')])
@pytest.mark.acceptance
def test_getplotmapWAQOS(file_nc):
    dir_output = getmakeoutputdir(__file__,inspect.currentframe().f_code.co_name)
    """
    file_nc = 'p:\\11201806-sophie\\Oosterschelde\\WAQ\\r03\\postprocessing\\oost_tracer_map_backup.nc'
    dir_output = '.'
    """

    import matplotlib.pyplot as plt
    plt.close('all')
    
    from dfm_tools.get_nc import get_netdata, get_ncmodeldata, plot_netmapdata

    ugrid = get_netdata(file_nc=file_nc)

    print('plot grid and values from mapdata (constantvalue, 1 dim)')
    if 'oost_tracer_map' in file_nc:
        var_names = ['mesh2d_FColi_1','mesh2d_HIWAI_1','mesh2d_mspaf_1','mesh2d_Pharma_1'] #nieuwe file, te veel dimensies
        var_clims = [None,[0,100000000000],None,[0,10000]]
    elif 'oost_tracer_2_map' in file_nc:
        var_names = ['mesh2d_FColi','mesh2d_HIWAI','mesh2d_Pharma'] #oude file
        var_clims = [None,[0,100000000000],[0,10000]]
    else:
        raise Exception('ERROR: no settings provided for this mapfile')

    for var_name, var_clim in zip(var_names, var_clims):
        fig, ax = plt.subplots()
        data_frommap = get_ncmodeldata(file_nc=file_nc, varname=var_name)#, multipart=False)
        pc = plot_netmapdata(ugrid.verts, values=data_frommap, ax=None, linewidth=0.5, cmap="jet")
        if var_clim != None:
            pc.set_clim(var_clim)
        fig.colorbar(pc, ax=ax)
        ax.set_aspect('equal')
        ax.set_xlabel(var_name)
        plt.savefig(os.path.join(dir_output,'%s_%s'%(os.path.basename(file_nc).replace('.nc',''),var_name)))
        




@pytest.mark.parametrize("file_nc", [pytest.param(os.path.join(dir_testinput,r'DFM_sigma_curved_bend\DFM_OUTPUT_cb_3d\cb_3d_map.nc'), id='cb_3d_map'),
                                     pytest.param(os.path.join(dir_testinput,r'DFM_3D_z_Grevelingen\computations\run01\DFM_OUTPUT_Grevelingen-FM\Grevelingen-FM_0000_map.nc'), id='Grevelingen-FM_0000_map'),
                                     #pytest.param(r'p:\11203379-mwra-new-bem-model\waq_model\simulations\A31_1year_20191219\DFM_OUTPUT_MB_02_waq\MB_02_waq_0000_map.nc', id='MB_02_waq_0000_map'),
                                     pytest.param(r'p:\1204257-dcsmzuno\2013-2017\3D-DCSM-FM\A17b\DFM_OUTPUT_DCSM-FM_0_5nm\DCSM-FM_0_5nm_0000_map.nc', id='DCSM-FM_0_5nm_0000_map'),
                                     pytest.param(r'p:\11205258-006-kpp2020_rmm-g6\C_Work\08_RMM_FMmodel\computations\run_156\DFM_OUTPUT_RMM_dflowfm\RMM_dflowfm_0000_map.nc', id='RMM_dflowfm_0000_map')])
@pytest.mark.acceptance
def test_getxzcoordsonintersection_plotcrossect(file_nc):

    dir_output = getmakeoutputdir(__file__,inspect.currentframe().f_code.co_name)
    """
    #manual test variables (run this script first to get the variable dir_testoutput)
    dir_output = '.'
    file_nc = os.path.join(dir_testinput,'DFM_sigma_curved_bend\\DFM_OUTPUT_cb_3d\\cb_3d_map.nc')
    file_nc = os.path.join(dir_testinput,'DFM_3D_z_Grevelingen','computations','run01','DFM_OUTPUT_Grevelingen-FM','Grevelingen-FM_0000_map.nc')
    file_nc = 'p:\\1204257-dcsmzuno\\2013-2017\\3D-DCSM-FM\\A17b\\DFM_OUTPUT_DCSM-FM_0_5nm\\DCSM-FM_0_5nm_0000_map.nc'
    file_nc = 'p:\\11205258-006-kpp2020_rmm-g6\\C_Work\\08_RMM_FMmodel\\computations\\run_156\\DFM_OUTPUT_RMM_dflowfm\\RMM_dflowfm_0000_map.nc'
    """
    
    import matplotlib.pyplot as plt
    plt.close('all')
    import numpy as np
    import datetime as dt
    
    from dfm_tools.get_nc import get_netdata, get_ncmodeldata, get_xzcoords_onintersection, plot_netmapdata
    from dfm_tools.polygon import LineBuilder
    
    
    if 'cb_3d_map' in file_nc:
        timestep = 72
        layno = 5
        calcdist_fromlatlon = None
        multipart = None
        line_array = np.array([[ 185.08667065, 2461.11775254],
                               [2934.63837418, 1134.16019127]])
        line_array = np.array([[ 104.15421399, 2042.7077107 ],
                               [2913.47878063, 2102.48057382]])
        val_ylim = None
        clim_bl = None
        #optimize_dist = None
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
    elif 'MB_02_waq_0000_map' in file_nc:
        timestep = 30
        layno = 5
        calcdist_fromlatlon = True
        multipart = None
        #provide xy order, so lonlat
        line_array = np.array([[-71.10395926,  42.3404146 ],
                               [-69.6762489 ,  42.38341792]])
        #line_array = np.array([[-70.87382752,  42.39103758], #dummy for partition 0000
        #                       [-70.42078633,  42.24876018]])
        val_ylim = None
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
        val_ylim = None
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
        #line_array = np.array([[ 88851.05823362, 413359.68286755], #dummy for partition 0000
        #                       [ 96948.34387646, 412331.45611925]])
        line_array = np.array([[129830.71514789, 425739.69372125], #waal
                               [131025.04347471, 425478.43439976],
                               [132126.06490098, 425758.35510136],
                               [133227.08632726, 426299.53512444],
                               [133824.25049067, 426504.81030561],
                               [134981.25605726, 426355.51926476],
                               [136810.07130769, 425329.14335891],
                               [137668.49479259, 425049.22265731],
                               [139534.63280323, 425403.78887934],
                               [140281.08800748, 425403.78887934],
                               [142464.46947993, 424620.01091487],
                               [143434.86124547, 424694.65643529],
                               [146271.39102164, 425534.41854008],
                               [148566.74077473, 426094.25994327]])
        val_ylim = None
        clim_bl = [-10,10]
        #optimize_dist = 150
    else:
        raise Exception('ERROR: no settings provided for this mapfile')
    
    
    ugrid = get_netdata(file_nc=file_nc, multipart=multipart)
    #get bed layer
    data_frommap_bl = get_ncmodeldata(file_nc=file_nc, varname='mesh2d_flowelem_bl', multipart=multipart)
    
    #create plot with ugrid and cross section line
    fig, ax_input = plt.subplots()
    pc = plot_netmapdata(ugrid.verts, values=data_frommap_bl, ax=ax_input, linewidth=0.5, edgecolors='face', cmap='jet')#, color='crimson', facecolor="None")
    pc.set_clim(clim_bl)
    fig.colorbar(pc, ax=ax_input)
    ax_input.set_aspect('equal')
    if 0: #click interactive polygon
        #pol_frominput = Polygon.frominteractive(ax)
        line, = ax_input.plot([], [],'o-')  # empty line
        linebuilder = LineBuilder(line)
        line_array = linebuilder.line_array
    ax_input.plot(line_array[:,0],line_array[:,1],'b',linewidth=3)
    
    
    runtime_tstart = dt.datetime.now() #start timer
    #intersect function, find crossed cell numbers (gridnos) and coordinates of intersection (2 per crossed cell)
    intersect_gridnos, intersect_coords = ugrid.polygon_intersect(line_array, optimize_dist=None)
    #derive vertices from cross section (distance from first point)
    crs_verts = get_xzcoords_onintersection(file_nc=file_nc, line_array=line_array, intersect_gridnos=intersect_gridnos, intersect_coords=intersect_coords, timestep=timestep, calcdist_fromlatlon=calcdist_fromlatlon, multipart=multipart)
    
    #get data to plot
    if 'DFM_OUTPUT_RMM_dflowfm' in file_nc:
        data_frommap = get_ncmodeldata(file_nc=file_nc, varname='mesh2d_sa1', timestep=timestep, multipart=multipart)
    else:
        data_frommap = get_ncmodeldata(file_nc=file_nc, varname='mesh2d_sa1', timestep=timestep, layer='all', multipart=multipart)
    
    #plot crossed cells (gridnos) in first plot
    print(layno)#data_frommap_flat = data_frommap[0,intersect_gridnos,layno]
    #pc = plot_netmapdata(ugrid.verts[intersect_gridnos,:,:], values=data_frommap_flat, ax=ax_input, linewidth=0.5, cmap="jet")
    plt.savefig(os.path.join(dir_output,'%s_gridbed'%(os.path.basename(file_nc).replace('.nc',''))))

    #plot cross section
    if len(data_frommap.shape) == 3:
        data_frommap_sel = data_frommap[0,intersect_gridnos,:]
        data_frommap_sel_flat = data_frommap_sel.T.flatten()
    elif len(data_frommap.shape) == 2: #for 2D models, no layers 
        data_frommap_sel = data_frommap[0,intersect_gridnos]
        data_frommap_sel_flat = data_frommap_sel
    fig, ax = plt.subplots()
    pc = plot_netmapdata(crs_verts, values=data_frommap_sel_flat, ax=ax, linewidth=0.5, cmap='jet')
    fig.colorbar(pc, ax=ax)
    ax.set_ylim(val_ylim)
    plt.savefig(os.path.join(dir_output,'%s_crossect'%(os.path.basename(file_nc).replace('.nc',''))))
    
    runtime_tstop = dt.datetime.now()
    runtime_timedelta = (runtime_tstop-runtime_tstart).total_seconds()
    print('calculating and plotting cross section finished in %.1f seconds'%(runtime_timedelta))





