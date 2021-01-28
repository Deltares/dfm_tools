# -*- coding: utf-8 -*-
"""
Created on Thu Feb 13 23:10:51 2020

@author: veenstra
"""

import pytest
import inspect
import os

from dfm_tools.testutils import getmakeoutputdir, gettestinputdir
dir_testinput = gettestinputdir()



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
def test_getncmodeldata_indexcountmetadata():
    from dfm_tools.get_nc import get_ncmodeldata
    
    #check if retrieving 1 index of data from 1 dimensional variable works (does not work if indices are np.arrays, so conversion to list in get_nc.py)
    file_his = r'p:\11205258-006-kpp2020_rmm-g6\C_Work\08_RMM_FMmodel\computations\run_161\DFM_OUTPUT_RMM_dflowfm\RMM_dflowfm_0000_his.nc'
    data_statcoord = get_ncmodeldata(file_nc=file_his, varname='station_x_coordinate',station='NM_1005.26_R_HBR-Cl_VCS-Nieuwe-Maas-85m')
    assert len(data_statcoord.var_stations) == 1
    
    file_his = os.path.join(dir_testinput,r'DFM_sigma_curved_bend\DFM_OUTPUT_cb_3d\cb_3d_his.nc')
    data_fromhis = get_ncmodeldata(file_nc=file_his, varname='x_velocity', timestep=1, layer=5, station='Innersouth boundary')#, multipart=False)
    assert len(data_fromhis.var_times) == 1
    assert len(data_fromhis.var_layers) == 1
    assert len(data_fromhis.var_stations) == 1
    
    #data_fromhis_all = get_ncmodeldata(file_nc=file_his, varname='x_velocity', timestep='all', layer='all', station='all')
    data_fromhis = get_ncmodeldata(file_nc=file_his, varname='x_velocity', timestep=[1,0,-5,1,4,3,2,0,1,-1], layer=[5,-2,3,0,5], station=[4,-1,0,0])
    assert len(data_fromhis.var_times) == 7
    assert data_fromhis.var_times.index.tolist() == [0,1,2,3,4,2156,2160]
    assert len(data_fromhis.var_layers) == 4
    assert data_fromhis.var_layers == [0,3,5,8]
    assert len(data_fromhis.var_stations) == 3
    assert data_fromhis.var_stations.index.tolist() == [0,4,5]
    
    



@pytest.mark.unittest
def test_getncmodeldata_datetime():
    import numpy as np
    import datetime as dt
    import pandas as pd
    
    from dfm_tools.get_nc import get_ncmodeldata

    #retrieve all
    file_nc = os.path.join(dir_testinput,r'DFM_sigma_curved_bend\DFM_OUTPUT_cb_3d\cb_3d_map.nc')
    data_frommap = get_ncmodeldata(file_nc=file_nc, varname='mesh2d_sa1', timestep='all', layer=5)#, multipart=False)
    assert data_frommap.shape[0] == len(data_frommap.var_times)
    
    #retrieve with numpy datetime array
    data_frommap = get_ncmodeldata(file_nc=file_nc, varname='mesh2d_sa1', 
                                   timestep=np.arange(dt.datetime(2001,1,1),dt.datetime(2001,1,2),dt.timedelta(hours=1)), layer=5)#, multipart=False)
    assert (data_frommap.data[0,0,0] - 31. ) < 1E-9
    assert data_frommap.shape[0] == len(data_frommap.var_times)
    
    #retrieve with pandas date_range
    data_frommap = get_ncmodeldata(file_nc=file_nc, varname='mesh2d_sa1', 
                                   timestep=pd.date_range(dt.datetime(2001,1,1),dt.datetime(2001,1,2),freq='30min'), layer=5)#, multipart=False)
    assert (data_frommap.data[0,0,0] - 31. ) < 1E-9
    assert data_frommap.shape[0] == len(data_frommap.var_times)
    


@pytest.mark.acceptance
def test_getplotfourstdata():
    dir_output = getmakeoutputdir(__file__,inspect.currentframe().f_code.co_name)
    """
    dir_output = './test_output'
    """
    import numpy as np
    import matplotlib.pyplot as plt
    plt.close('all')
    
    from dfm_tools.get_nc import get_netdata, get_ncmodeldata, plot_netmapdata
    from dfm_tools.get_nc_helpers import get_ncvardimlist
    from dfm_tools.regulargrid import scatter_to_regulargrid
    
    #RMM foufile met quivers
    file_nc = os.path.join(dir_testinput,r'DFM_fou_RMM\RMM_dflowfm_0000_fou.nc')
    vars_pd, dims_pd = get_ncvardimlist(file_nc=file_nc)
    #stations_pd = get_hisstationlist(file_nc,varname='waterlevel')
    
    ugrid_all = get_netdata(file_nc=file_nc)
    ux_mean = get_ncmodeldata(file_nc=file_nc, varname='mesh2d_fourier001_mean')
    uy_mean = get_ncmodeldata(file_nc=file_nc, varname='mesh2d_fourier002_mean')
    magn_mean = np.sqrt(ux_mean**2+uy_mean**2)
    #uc_mean = get_ncmodeldata(file_nc=file_nc, varname='mesh2d_fourier003_mean')
    #uc_max = get_ncmodeldata(file_nc=file_nc, varname='mesh2d_fourier004_max')
    facex = get_ncmodeldata(file_nc=file_nc, varname='mesh2d_face_x')
    facey = get_ncmodeldata(file_nc=file_nc, varname='mesh2d_face_y')
    
    X,Y,U = scatter_to_regulargrid(xcoords=facex, ycoords=facey, ncellx=60, ncelly=35, values=ux_mean)
    X,Y,V = scatter_to_regulargrid(xcoords=facex, ycoords=facey, ncellx=60, ncelly=35, values=uy_mean)
    
    #thinning = 3
    fig1,ax1 = plt.subplots(figsize=(9,5))
    pc1 = plot_netmapdata(ugrid_all.verts, magn_mean, edgecolor='face')
    #ax1.quiver(facex[::thinning], facey[::thinning], ux_mean[::thinning], uy_mean[::thinning], color='w',scale=20)#,width=0.005)#, edgecolor='face', cmap='jet')
    ax1.quiver(X,Y,U,V, color='w',scale=5)#,width=0.005)#, edgecolor='face', cmap='jet')
    pc1.set_clim([0,0.10])
    #ax1.set_title('sqrt(x^2+y^2)\nx=%s\ny=%s'%(ux_mean.var_ncvarobject.long_name,uy_mean.var_ncvarobject.long_name))
    ax1.set_aspect('equal')
    ax1.set_xlabel('RD x [m]')
    ax1.set_ylabel('RD y [m]')
    cbar = fig1.colorbar(pc1)
    cbar.set_label('residuele stroming [m/s]')
    fig1.tight_layout()
    fig1.savefig(os.path.join(dir_output,os.path.basename(file_nc).replace('.','')))
    
    
    #RMM rst file
    file_nc = os.path.join(dir_testinput,r'DFM_fou_RMM\RMM_dflowfm_0006_20131127_000000_rst.nc')
    vars_pd, dims_pd = get_ncvardimlist(file_nc=file_nc)
    #ugrid = get_netdata(file_nc=file_nc) #does not work, so scatter has to be used
    ugrid_FlowElem_xzw = get_ncmodeldata(file_nc=file_nc, varname='FlowElem_xzw', multipart=False)
    ugrid_FlowElem_yzw = get_ncmodeldata(file_nc=file_nc, varname='FlowElem_yzw', multipart=False)
    data_s1 = get_ncmodeldata(file_nc=file_nc, varname='s1',timestep=0, multipart=False)
    
    fig, ax = plt.subplots(figsize=(10,4))
    pc = plt.scatter(ugrid_FlowElem_xzw,ugrid_FlowElem_yzw,1, data_s1[0,:],cmap='jet')
    pc.set_clim([0,2])
    fig.colorbar(pc)
    ax.set_aspect('equal')
    fig.tight_layout()
    plt.savefig(os.path.join(dir_output,os.path.basename(file_nc).replace('.','')))
    
    #keep this test to check flow link handling by get_nc.py (does not have domains)
    #ugrid_FlowLink_xu = get_ncmodeldata(file_nc=file_nc, varname='FlowLink_xu', multipart=False)
    #ugrid_FlowLink_yu = get_ncmodeldata(file_nc=file_nc, varname='FlowLink_yu', multipart=False)
    data_q1 = get_ncmodeldata(file_nc=file_nc, varname='q1',timestep=0, multipart=False)
    
    """
    fig, ax = plt.subplots(figsize=(10,4))
    pc = plt.scatter(ugrid_FlowLink_xu,ugrid_FlowLink_yu,1, data_q1[0,:],cmap='jet')
    pc.set_clim(None)
    fig.colorbar(pc)
    ax.set_aspect('equal')
    fig.tight_layout()
    plt.savefig(os.path.join(dir_output,os.path.basename(file_nc).replace('.','')))
    """

    assert len(data_q1[0,:]) == 138756






    
@pytest.mark.parametrize("file_nc", [pytest.param(os.path.join(dir_testinput,r'DFM_3D_z_Grevelingen\computations\run01\Grevelingen_FM_grid_20190603_net.nc'), id='Grevelingen'),
                                     pytest.param(os.path.join(dir_testinput,'vanNithin','myortho3_RGFGRID_net.nc'), id='Nithin'),
                                     pytest.param(r'p:\11205258-006-kpp2020_rmm-g6\C_Work\01_Rooster\final_totaalmodel\rooster_rmm_v1p5_net.nc', id='RMM')])
@pytest.mark.acceptance
def test_getnetdata_plotnet(file_nc):
    dir_output = getmakeoutputdir(__file__,inspect.currentframe().f_code.co_name)
    """
    this test retrieves grid data and plots it
    
    file_nc = os.path.join(dir_testinput,'DFM_3D_z_Grevelingen','computations','run01','Grevelingen_FM_grid_20190603_net.nc')
    file_nc = 'p:\\11205258-006-kpp2020_rmm-g6\\C_Work\\01_Rooster\\final_totaalmodel\\rooster_rmm_v1p5_net.nc'
    dir_output = './test_output'
    """
    

    import matplotlib.pyplot as plt
    plt.close('all')

    from dfm_tools.get_nc import get_netdata, plot_netmapdata

    print('plot only grid from net.nc')
    ugrid = get_netdata(file_nc=file_nc)
    fig, ax = plt.subplots()
    plot_netmapdata(ugrid.verts, values=None, ax=None, linewidth=0.5, color="crimson", facecolor="None")
    ax.set_aspect('equal')
    plt.savefig(os.path.join(dir_output,os.path.basename(file_nc).replace('.','')))
    
    




@pytest.mark.parametrize("file_nc", [pytest.param(r'p:\11203869-morwaqeco3d\05-Tidal_inlet\02_FM_201910\FM_MF10_Max_30s\wave\wavm-inlet.nc', id='Tidal_inlet_wavm'),
                                     #pytest.param(r'p:\11200665-c3s-codec\2_Hydro\ECWMF_meteo\meteo\ERA-5\2000\ERA5_metOcean_atm_19991201_19991231.nc', id='ERA5_meteo'),
                                     pytest.param(r'p:\1204257-dcsmzuno\2014\data\meteo\HIRLAM72_2018\h72_201803.nc', id='hirlam_meteo'),
                                     pytest.param(r'p:\11202255-sfincs\Testbed\Original_runs\01_Implementation\08_restartfile\sfincs_map.nc', id='sfincs_map')])
@pytest.mark.acceptance
def test_getnetdata_plotnet_regular(file_nc):

    dir_output = getmakeoutputdir(__file__,inspect.currentframe().f_code.co_name)
    """
    dir_output = './test_output'
    file_nc = 'p:\\11203869-morwaqeco3d\\05-Tidal_inlet\\02_FM_201910\\FM_MF10_Max_30s\\wave\\wavm-inlet.nc'
    file_nc = 'p:\\11200665-c3s-codec\\2_Hydro\\ECWMF_meteo\\meteo\\ERA-5\\2000\\ERA5_metOcean_atm_19991201_19991231.nc'
    file_nc = 'p:\\1204257-dcsmzuno\\2014\\data\\meteo\\HIRLAM72_2018\\h72_201803.nc'
    file_nc = r'p:\11202255-sfincs\Testbed\Original_runs\01_Implementation\08_restartfile\sfincs_map.nc'
    """
    
    import numpy as np
    import matplotlib.pyplot as plt
    plt.close('all')
    
    from dfm_tools.get_nc import get_ncmodeldata, plot_netmapdata#, get_xzcoords_onintersection
    from dfm_tools.regulargrid import meshgridxy2verts, center2corner
    from dfm_tools.get_nc_helpers import get_ncvardimlist
    
    
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
    


    
    
    
@pytest.mark.acceptance
def test_getsobekmodeldata():
    dir_output = getmakeoutputdir(__file__,inspect.currentframe().f_code.co_name)
    """
    this test retrieves sobek observation data and plots it
    
    dir_output = './test_output'
    """

    import matplotlib.pyplot as plt
    plt.close('all')

    from dfm_tools.get_nc import get_ncmodeldata
    #from dfm_tools.get_nc_helpers import get_hisstationlist
    
    file_nc = os.path.join(dir_testinput,'KenmerkendeWaarden','observations.nc')
    
    #station_names = get_hisstationlist(file_nc=file_nc, varname='observation_id')
    #station_names = get_hisstationlist(file_nc=file_nc, varname='water_level')
    data_fromsobek = get_ncmodeldata(file_nc=file_nc, varname='water_level', station=['Maasmond','HKVHLD','MO_1035.00'], timestep='all')
    
    fig, ax = plt.subplots()
    for iL in range(data_fromsobek.shape[1]):
        ax.plot(data_fromsobek.var_times,data_fromsobek[:,iL],'-', label=data_fromsobek.var_stations.iloc[iL,0])
    ax.legend()
    plt.savefig(os.path.join(dir_output,'%s_waterlevel'%(os.path.basename(file_nc).replace('.',''))))
    
    
    
    
    
    
@pytest.mark.parametrize("file_nc", [pytest.param(os.path.join(dir_testinput,'vanNithin','tttz_0000_his.nc'), id='Nithin'),
                                     pytest.param(os.path.join(dir_testinput,'DFM_3D_z_Grevelingen\\computations\\run01\\DFM_OUTPUT_Grevelingen-FM\\Grevelingen-FM_0000_his.nc'), id='Grevelingen')])
@pytest.mark.acceptance
def test_gethismodeldata(file_nc):
    dir_output = getmakeoutputdir(__file__,inspect.currentframe().f_code.co_name)
    """
    this test retrieves his data and plots it
    file_nc = os.path.join(dir_testinput,'vanNithin','tttz_0000_his.nc')
    file_nc = os.path.join(dir_testinput,'DFM_3D_z_Grevelingen\\computations\\run01\\DFM_OUTPUT_Grevelingen-FM\\Grevelingen-FM_0000_his.nc')
    file_nc = r'p:\11202512-h2020_impaqt\Mediterranean_model\MedSea_impaqt_model\computations\r003_test\DFM_OUTPUT_MedSea_impaqt_FM\MedSea_impaqt_FM_0000_his.nc'
    dir_output = './test_output'
    """

    import numpy as np
    import matplotlib.pyplot as plt
    plt.close('all')
    
    from dfm_tools.get_nc import get_ncmodeldata, plot_ztdata
    
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
        station = ['Peiraias', 'Ovrios_2','Ovrios','Ovrios','Ortholithi']
        station_zt = 'Ortholithi'
    elif 'impaqt' in file_nc:
        station_zt = 'IOC_thes' #['IOC_thes','farm_impaqt']

    print('plot bedlevel from his')
    data_fromhis = get_ncmodeldata(file_nc=file_nc, varname='bedlevel', station=station)#, multipart=False)
    fig, ax = plt.subplots()
    ax.plot(data_fromhis.var_stations.iloc[:,0],data_fromhis,'-')
    ax.tick_params('x',rotation=90)
    fig.savefig(os.path.join(dir_output,'%s_bedlevel'%(os.path.basename(file_nc).replace('.',''))))

    print('plot waterlevel from his')
    data_fromhis = get_ncmodeldata(file_nc=file_nc, varname='waterlevel', timestep='all', station=station)#, multipart=False)
    fig, ax = plt.subplots()
    ax.plot(data_fromhis.var_times,data_fromhis,'-')
    fig.savefig(os.path.join(dir_output,'%s_waterlevel'%(os.path.basename(file_nc).replace('.',''))))
    
    print('plot salinity from his')
    data_fromhis = get_ncmodeldata(file_nc=file_nc, varname='salinity', timestep='all', layer=5, station=station)#, multipart=False)
    data_fromhis_flat = data_fromhis[:,:,0]
    fig, ax = plt.subplots()
    ax.plot(data_fromhis.var_times,data_fromhis_flat,'-')
    fig.savefig(os.path.join(dir_output,'%s_salinity'%(os.path.basename(file_nc).replace('.',''))))

    print('plot salinity over depth')
    #depth retrieval is probably wrong
    data_fromhis_depth = get_ncmodeldata(file_nc=file_nc, varname='zcoordinate_c', timestep=4, layer='all', station=station)
    data_fromhis = get_ncmodeldata(file_nc=file_nc, varname='salinity', timestep=4, layer='all', station=station)
    fig, ax = plt.subplots()
    ax.plot(data_fromhis[0,:,:].T, data_fromhis_depth[0,:,:].T,'-')
    ax.legend(data_fromhis.var_stations.iloc[:,0])
    fig.savefig(os.path.join(dir_output,'%s_salinityoverdepth'%(os.path.basename(file_nc).replace('.',''))))

    print('zt temperature plot and wl')
    data_fromhis_temp = get_ncmodeldata(file_nc=file_nc, varname='temperature', timestep=range(40,100), layer= 'all', station=station_zt)
    data_fromhis_wl = get_ncmodeldata(file_nc=file_nc, varname='waterlevel', timestep=range(40,100), station=station_zt)#, multipart=False)
    fig, (axwl,ax1) = plt.subplots(2,1,figsize=(12,7),gridspec_kw={'height_ratios':[1,2]},sharex=True,sharey=True)
    statid_subset = data_fromhis_wl.var_stations['station_id'].tolist().index(station_zt)
    axwl.plot(data_fromhis_wl.var_times.iloc[[0,-1]],[0,0],'k-',linewidth=0.5)
    ax1.plot(data_fromhis_wl.var_times.iloc[[0,-1]],[0,0],'k-',linewidth=0.5)
    axwl.plot(data_fromhis_wl.var_times,data_fromhis_wl[:,statid_subset],'-',label='wl %s'%(station_zt))
    c = plot_ztdata(dfmtools_hisvar=data_fromhis_temp, ax=ax1, statid_subset=statid_subset, cmap='jet')
    fig.colorbar(c,ax=axwl)
    fig.colorbar(c,ax=ax1)
    #contour
    CS = plot_ztdata(dfmtools_hisvar=data_fromhis_temp, ax=ax1, statid_subset=statid_subset, only_contour=True, levels=6, colors='k', linewidths=0.8, linestyles='solid')
    ax1.clabel(CS, fontsize=10)
    fig.tight_layout()
    fig.savefig(os.path.join(dir_output,'%s_zt_temp'%(os.path.basename(file_nc).replace('.',''))))
    axwl.set_ylim(-2,0.5)
    fig.savefig(os.path.join(dir_output,'%s_zt_temp_zoomwl'%(os.path.basename(file_nc).replace('.',''))))

    
    

@pytest.mark.parametrize("file_nc", [pytest.param(os.path.join(dir_testinput,r'DFM_sigma_curved_bend\DFM_OUTPUT_cb_3d\cb_3d_map.nc'), id='curvibend'),
                                     pytest.param(os.path.join(dir_testinput,r'DFM_3D_z_Grevelingen\computations\run01\DFM_OUTPUT_Grevelingen-FM\Grevelingen-FM_0000_map.nc'), id='Grevelingen'),
                                     pytest.param(r'p:\11205258-006-kpp2020_rmm-g6\C_Work\08_RMM_FMmodel\computations\run_180\DFM_OUTPUT_RMM_dflowfm\RMM_dflowfm_0000_map.nc', id='RMM')])
@pytest.mark.acceptance
def test_getnetdata_getmapmodeldata_plotnetmapdata(file_nc):
    dir_output = getmakeoutputdir(__file__,inspect.currentframe().f_code.co_name)
    """
    this test retrieves grid data, retrieves map data, and plots it
    file_nc = os.path.join(dir_testinput,'DFM_3D_z_Grevelingen','computations','run01','DFM_OUTPUT_Grevelingen-FM','Grevelingen-FM_0000_map.nc')
    file_nc = r'p:\11205258-006-kpp2020_rmm-g6\C_Work\08_RMM_FMmodel\computations\run_180\DFM_OUTPUT_RMM_dflowfm\RMM_dflowfm_0000_map.nc'
    dir_output = './test_output'
    """

    import matplotlib.pyplot as plt
    plt.close('all')
    
    from dfm_tools.get_nc import get_netdata, get_ncmodeldata, plot_netmapdata
    from dfm_tools.get_nc_helpers import get_ncvardimlist
    
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
        
    
    vars_pd, dims_pd = get_ncvardimlist(file_nc=file_nc)

    #PLOT GRID
    print('plot only grid from mapdata')
    ugrid_all = get_netdata(file_nc=file_nc)#,multipart=False)
    fig, ax = plt.subplots()
    pc = plot_netmapdata(ugrid_all.verts, values=None, ax=None, linewidth=0.5, color="crimson", facecolor="None")
    ax.set_xlabel('x-direction')
    ax.set_ylabel('y-direction')
    ax.set_aspect('equal')
    fig.tight_layout()
    fig.savefig(os.path.join(dir_output,'%s_grid'%(os.path.basename(file_nc).replace('.',''))))


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
        cbar = fig.colorbar(pc, ax=ax)
        cbar.set_label('%s [%s]'%(data_frommap.var_ncvarobject.long_name, data_frommap.var_ncvarobject.units))
        varcoords_x = data_frommap.var_ncobject.variables.get(data_frommap.var_ncvarobject.coordinates.split()[0])
        varcoords_y = data_frommap.var_ncobject.variables.get(data_frommap.var_ncvarobject.coordinates.split()[1])
        ax.set_xlabel('%s [%s]'%(varcoords_x.long_name,varcoords_x.units))
        ax.set_ylabel('%s [%s]'%(varcoords_y.long_name,varcoords_y.units))
        ax.set_aspect('equal')
        fig.tight_layout()
        fig.savefig(os.path.join(dir_output,'%s_mesh2d_flowelem_bl'%(os.path.basename(file_nc).replace('.',''))))
        
    
    #PLOT water level on map
    print('plot grid and values from mapdata (waterlevel, 2dim)')
    data_frommap = get_ncmodeldata(file_nc=file_nc, varname='mesh2d_s1', timestep=timestep)#, multipart=False)
    data_frommap_flat = data_frommap.flatten()
    fig, ax = plt.subplots()
    pc = plot_netmapdata(ugrid_all.verts, values=data_frommap_flat, ax=None, linewidth=0.5, cmap="jet")
    pc.set_clim(clim_wl)
    cbar = fig.colorbar(pc, ax=ax)
    cbar.set_label('%s [%s]'%(data_frommap.var_ncvarobject.long_name, data_frommap.var_ncvarobject.units))
    varcoords_x = data_frommap.var_ncobject.variables.get(data_frommap.var_ncvarobject.coordinates.split()[0])
    varcoords_y = data_frommap.var_ncobject.variables.get(data_frommap.var_ncvarobject.coordinates.split()[1])
    ax.set_xlabel('%s [%s]'%(varcoords_x.long_name,varcoords_x.units))
    ax.set_ylabel('%s [%s]'%(varcoords_y.long_name,varcoords_y.units))
    ax.set_aspect('equal')
    fig.tight_layout()
    fig.savefig(os.path.join(dir_output,'%s_mesh2d_s1'%(os.path.basename(file_nc).replace('.',''))))

    #PLOT var layer on map
    if 'RMM_dflowfm_0000_map' in file_nc:
        print('plot grid and values from mapdata (wind x velocity on cell edges)')
        data_frommap = get_ncmodeldata(file_nc=file_nc, varname='mesh2d_windxu', timestep=timestep, layer=layer)#, multipart=False)
        data_frommap_flat = data_frommap.flatten()
        fig, ax = plt.subplots()
        pc = plot_netmapdata(ugrid_all.edge_verts, values=data_frommap_flat, ax=None, linewidth=0.5, cmap="jet", edgecolor='face')
        #pc.set_clim(0,5)
        cbar = fig.colorbar(pc, ax=ax)
        cbar.set_label('%s [%s]'%(data_frommap.var_ncvarobject.long_name, data_frommap.var_ncvarobject.units))
        varcoords_x = data_frommap.var_ncobject.variables.get(data_frommap.var_ncvarobject.coordinates.split()[0])
        varcoords_y = data_frommap.var_ncobject.variables.get(data_frommap.var_ncvarobject.coordinates.split()[1])
        ax.set_xlabel('%s [%s]'%(varcoords_x.long_name,varcoords_x.units))
        ax.set_ylabel('%s [%s]'%(varcoords_y.long_name,varcoords_y.units))
        ax.set_aspect('equal')
        fig.tight_layout()
        fig.savefig(os.path.join(dir_output,'%s_mesh2d_windxu_edges'%(os.path.basename(file_nc).replace('.',''))))
    else:
        print('plot grid and values from mapdata (salinity on layer, 3dim)')
        data_frommap = get_ncmodeldata(file_nc=file_nc, varname='mesh2d_sa1', timestep=timestep, layer=layer)#, multipart=False)
        data_frommap_flat = data_frommap.flatten()
        fig, ax = plt.subplots()
        pc = plot_netmapdata(ugrid_all.verts, values=data_frommap_flat, ax=None, linewidth=0.5, cmap="jet")
        pc.set_clim(clim_sal)
        cbar = fig.colorbar(pc, ax=ax)
        cbar.set_label('%s [%s]'%(data_frommap.var_ncvarobject.long_name, data_frommap.var_ncvarobject.units))
        varcoords_x = data_frommap.var_ncobject.variables.get(data_frommap.var_ncvarobject.coordinates.split()[0])
        varcoords_y = data_frommap.var_ncobject.variables.get(data_frommap.var_ncvarobject.coordinates.split()[1])
        ax.set_xlabel('%s [%s]'%(varcoords_x.long_name,varcoords_x.units))
        ax.set_ylabel('%s [%s]'%(varcoords_y.long_name,varcoords_y.units))
        ax.set_aspect('equal')
        fig.tight_layout()
        fig.savefig(os.path.join(dir_output,'%s_mesh2d_sa1'%(os.path.basename(file_nc).replace('.',''))))
    
        print('plot grid and values from mapdata (temperature on layer, 3dim)')
        data_frommap = get_ncmodeldata(file_nc=file_nc, varname='mesh2d_tem1', timestep=timestep, layer=layer)#, multipart=False)
        data_frommap_flat = data_frommap.flatten()
        fig, ax = plt.subplots()
        pc = plot_netmapdata(ugrid_all.verts, values=data_frommap_flat, ax=None, linewidth=0.5, cmap="jet")
        pc.set_clim(clim_tem)
        cbar = fig.colorbar(pc, ax=ax)
        cbar.set_label('%s [%s]'%(data_frommap.var_ncvarobject.long_name, data_frommap.var_ncvarobject.units))
        varcoords_x = data_frommap.var_ncobject.variables.get(data_frommap.var_ncvarobject.coordinates.split()[0])
        varcoords_y = data_frommap.var_ncobject.variables.get(data_frommap.var_ncvarobject.coordinates.split()[1])
        ax.set_xlabel('%s [%s]'%(varcoords_x.long_name,varcoords_x.units))
        ax.set_ylabel('%s [%s]'%(varcoords_y.long_name,varcoords_y.units))
        ax.set_aspect('equal')
        fig.tight_layout()
        fig.savefig(os.path.join(dir_output,'%s_mesh2d_tem1'%(os.path.basename(file_nc).replace('.',''))))
  
    
  
    
  
  
    
@pytest.mark.acceptance
def test_contextily_addbasemap():
    dir_output = getmakeoutputdir(__file__,inspect.currentframe().f_code.co_name)
    
    """
    https://contextily.readthedocs.io/en/latest/reference.html
    https://contextily.readthedocs.io/en/latest/intro_guide.html
    ctx.add_basemap() defaults:
        source: None defaults to ctx.providers.Stamen.Terrain
        crs: coordinate reference system (CRS). If None (default), no warping is performed and the original Spherical Mercator (EPSG:3857) is used.
    
    dir_output = './test_output'
    """
    
    import matplotlib.pyplot as plt
    plt.close('all')

    from dfm_tools.testutils import try_importmodule
    try_importmodule(modulename='contextily') #check if contextily was installed since it is an optional module, also happens in plot_cartopybasemap()
    import contextily as ctx

    from dfm_tools.get_nc import get_netdata, get_ncmodeldata, plot_netmapdata

    file_nc_map = os.path.join(dir_testinput,'DFM_3D_z_Grevelingen\\computations\\run01\\DFM_OUTPUT_Grevelingen-FM\\Grevelingen-FM_0000_map.nc')
    ugrid = get_netdata(file_nc=file_nc_map)
    data_frommap_bl = get_ncmodeldata(file_nc=file_nc_map, varname='mesh2d_flowelem_bl')
    
    source_list = [ctx.providers.Stamen.Terrain, #default source
                   ctx.providers.Esri.WorldImagery,
                   ctx.providers.CartoDB.Voyager,
                   #ctx.providers.NASAGIBS.ViirsEarthAtNight2012,
                   ctx.providers.Stamen.Watercolor]
    
    for source_ctx in source_list:
        source_name = source_ctx['name'].replace('.','_')
        fig, ax = plt.subplots(1,1,figsize=(10,6))
        pc = plot_netmapdata(ugrid.verts, values=data_frommap_bl, ax=ax, linewidth=0.5, cmap='jet')
        fig.colorbar(pc, ax=ax)
        fig.tight_layout()
        ctx.add_basemap(ax, source=source_ctx, crs="EPSG:28992", attribution_size=5)
        fig.savefig(os.path.join(dir_output,'contextily_grevelingen_RD_%s'%(source_name)))
    





@pytest.mark.systemtest
def test_cartopy_epsg():
    
    from dfm_tools.testutils import try_importmodule
    try_importmodule(modulename='cartopy') #check if cartopy was installed since it is an optional module, also happens in plot_cartopybasemap()
    
    from dfm_tools.get_nc import plot_background
    
    #this one crashes if the dummy in plot_background() is not created
    plot_background(ax=None, projection=28992, google_style='satellite', resolution=5, features='land', nticks=6, latlon_format=False, gridlines=False)





    
    
@pytest.mark.acceptance
def test_cartopy_satellite_coastlines():
    dir_output = getmakeoutputdir(__file__,inspect.currentframe().f_code.co_name)
    
    """
    dir_output = './test_output'
    """
    import matplotlib.pyplot as plt
    plt.close('all')
    import numpy as np

    from dfm_tools.testutils import try_importmodule
    try_importmodule(modulename='cartopy') #check if cartopy was installed since it is an optional module, also happens in plot_cartopybasemap()
    import cartopy.crs as ccrs

    from dfm_tools.get_nc import get_netdata, get_ncmodeldata, plot_netmapdata, plot_background
    from dfm_tools.get_nc_helpers import get_ncvardimlist
    
    #HIRLAM
    file_nc = r'p:\1204257-dcsmzuno\2014\data\meteo\HIRLAM72_2018\h72_201803.nc'
    vars_pd, dims_pd = get_ncvardimlist(file_nc=file_nc)
    
    timestep = 0
    mesh2d_node_x = get_ncmodeldata(file_nc=file_nc, varname='x')
    mesh2d_node_y = get_ncmodeldata(file_nc=file_nc, varname='y')
    mesh2d_node_x_sel = mesh2d_node_x[::2,::2]
    mesh2d_node_y_sel = mesh2d_node_y[::2,::2]
    data_v = get_ncmodeldata(file_nc=file_nc, varname='northward_wind',timestep=timestep)
    data_u = get_ncmodeldata(file_nc=file_nc, varname='eastward_wind',timestep=timestep)
    #airp = get_ncmodeldata(file_nc=file_nc, varname='air_pressure_fixed_height',timestep=0)[0,:,:]
    magn = np.sqrt(data_u**2 + data_v**2)[0,::2,::2]
    
    fig, ax = plt.subplots()
    ax.pcolor(mesh2d_node_x_sel,mesh2d_node_y_sel,magn)
    plt.savefig(os.path.join(dir_output,'cartopy_hirlam_raw'))
    
    fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()})
    ax.pcolor(mesh2d_node_x_sel,mesh2d_node_y_sel,magn)#, transform=ccrs.PlateCarree())
    plt.savefig(os.path.join(dir_output,'cartopy_hirlam_aspect'))
    
    fig, ax = plt.subplots(figsize=(9,5),subplot_kw={'projection': ccrs.PlateCarree()}) #provide axis projection on initialisation, cannot be edited later on
    pc = ax.pcolor(mesh2d_node_x_sel,mesh2d_node_y_sel,magn)#, transform=ccrs.PlateCarree())
    cbar = fig.colorbar(pc, ax=ax)
    cbar.set_label('velocity magnitude (%s)'%(data_v.var_ncvarobject.units))
    plot_background(ax=ax, resolution=1, google_style='street', features=['countries_highres'], linewidth=0.5, edgecolor='gray', facecolor='none', latlon_format=True)
    plot_background(ax=ax, google_style=None, features=['coastlines_highres'], linewidth=0.5, latlon_format=True)
    plt.savefig(os.path.join(dir_output,'cartopy_hirlam_moreoptions'))
    
    fig, ax = plt.subplots(figsize=(6,7),subplot_kw={'projection': ccrs.EuroPP()}) #provide axis projection on initialisation, cannot be edited later on
    pc = ax.pcolor(mesh2d_node_x_sel[:100,:100],mesh2d_node_y_sel[:100,:100],magn[:100,:100], transform=ccrs.PlateCarree()) #take subset of dataset to speed up coordinate transformation
    plot_background(ax=ax, google_style=None, features=['coastlines_highres'], latlon_format=True, gridlines=True)
    plt.savefig(os.path.join(dir_output,'cartopy_hirlam_curvedgridlines'))
    
        
    #GREVELINGEN
    file_nc_map = os.path.join(dir_testinput,'DFM_3D_z_Grevelingen\\computations\\run01\\DFM_OUTPUT_Grevelingen-FM\\Grevelingen-FM_0000_map.nc')
    ugrid = get_netdata(file_nc=file_nc_map)
    data_frommap_bl = get_ncmodeldata(file_nc=file_nc_map, varname='mesh2d_flowelem_bl')
    
    fig, ax = plt.subplots(1,1, subplot_kw={'projection': ccrs.epsg(28992)}) #provide axis projection on initialisation, cannot be edited later on
    pc = plot_netmapdata(ugrid.verts, values=data_frommap_bl, ax=ax, linewidth=0.5, cmap='jet')
    plot_background(ax=ax, resolution=12, features=['coastlines_highres'], linewidth=0.5)
    plt.savefig(os.path.join(dir_output,'cartopy_grevelingen_RD'))

    







@pytest.mark.parametrize("file_nc", [pytest.param('p:\\11201806-sophie\\Oosterschelde\\WAQ\\r02\\postprocessing\\oost_tracer_2_map.nc', id='oost_tracer_2_map')])
@pytest.mark.acceptance
def test_getplotmapWAQOS(file_nc):
    dir_output = getmakeoutputdir(__file__,inspect.currentframe().f_code.co_name)
    """
    file_nc = 'p:\\11201806-sophie\\Oosterschelde\\WAQ\\r03\\postprocessing\\oost_tracer_map.nc' #constantly changes name and dimensions, removed from testbank
    file_nc = 'p:\\11201806-sophie\\Oosterschelde\\WAQ\\r02\\postprocessing\\oost_tracer_2_map.nc'
    dir_output = './test_output'
    """

    import matplotlib.pyplot as plt
    plt.close('all')
    
    from dfm_tools.get_nc import get_netdata, get_ncmodeldata, plot_netmapdata

    ugrid = get_netdata(file_nc=file_nc)

    print('plot grid and values from mapdata (constantvalue, 1 dim)')
    if 'oost_tracer_map' in file_nc:
        var_names = ['FColi1','HIWAI1','mspaf1','Pharma1'] #nieuwe file, te veel dimensies
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
        plt.savefig(os.path.join(dir_output,'%s_%s'%(os.path.basename(file_nc).replace('.',''),var_name)))
        




@pytest.mark.parametrize("file_nc", [pytest.param(os.path.join(dir_testinput,r'DFM_sigma_curved_bend\DFM_OUTPUT_cb_3d\cb_3d_map.nc'), id='cb_3d_map'),
                                     pytest.param(os.path.join(dir_testinput,r'DFM_3D_z_Grevelingen\computations\run01\DFM_OUTPUT_Grevelingen-FM\Grevelingen-FM_0000_map.nc'), id='Grevelingen-FM_0000_map'),
                                     #pytest.param(r'p:\11203379-mwra-new-bem-model\waq_model\simulations\A31_1year_20191219\DFM_OUTPUT_MB_02_waq\MB_02_waq_0000_map.nc', id='MB_02_waq_0000_map'),
                                     pytest.param(r'p:\1204257-dcsmzuno\2013-2017\3D-DCSM-FM\A19\DFM_OUTPUT_DCSM-FM_0_5nm\DCSM-FM_0_5nm_0000_map.nc', id='DCSM-FM_0_5nm_0000_map'),
                                     pytest.param(r'p:\11205258-006-kpp2020_rmm-g6\C_Work\08_RMM_FMmodel\computations\run_180\DFM_OUTPUT_RMM_dflowfm\RMM_dflowfm_0000_map.nc', id='RMM_dflowfm_0000_map')])
@pytest.mark.acceptance
def test_getxzcoordsonintersection_plotcrossect(file_nc):

    dir_output = getmakeoutputdir(__file__,inspect.currentframe().f_code.co_name)
    """
    #manual test variables (run this script first to get the variable dir_testoutput)
    dir_output = './test_output'
    file_nc = os.path.join(dir_testinput,'DFM_sigma_curved_bend\\DFM_OUTPUT_cb_3d\\cb_3d_map.nc')
    file_nc = os.path.join(dir_testinput,'DFM_3D_z_Grevelingen','computations','run01','DFM_OUTPUT_Grevelingen-FM','Grevelingen-FM_0000_map.nc')
    file_nc = 'p:\\1204257-dcsmzuno\\2013-2017\\3D-DCSM-FM\\A19\\DFM_OUTPUT_DCSM-FM_0_5nm\\DCSM-FM_0_5nm_0000_map.nc'
    file_nc = 'p:\\11205258-006-kpp2020_rmm-g6\\C_Work\\08_RMM_FMmodel\\computations\\run_180\\DFM_OUTPUT_RMM_dflowfm\\RMM_dflowfm_0000_map.nc'
    """
    
    import matplotlib.pyplot as plt
    plt.close('all')
    import numpy as np
    import datetime as dt
    
    from dfm_tools.get_nc import get_netdata, get_ncmodeldata, get_xzcoords_onintersection, plot_netmapdata
    from dfm_tools.io.polygon import LineBuilder
    
    
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
    plt.savefig(os.path.join(dir_output,'%s_gridbed'%(os.path.basename(file_nc).replace('.',''))))

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
    plt.savefig(os.path.join(dir_output,'%s_crossect'%(os.path.basename(file_nc).replace('.',''))))
    
    runtime_tstop = dt.datetime.now()
    runtime_timedelta = (runtime_tstop-runtime_tstart).total_seconds()
    print('calculating and plotting cross section finished in %.1f seconds'%(runtime_timedelta))


 
    






@pytest.mark.acceptance
def test_readme_example_usage():
    #import statements
    import os
    import matplotlib.pyplot as plt
    plt.close('all')
    from dfm_tools.get_nc import get_netdata, get_ncmodeldata, plot_netmapdata
    from dfm_tools.get_nc_helpers import get_ncvardimlist, get_timesfromnc, get_hisstationlist
    
    #uncomment the line below, copy data locally and change this path to increase performance
    #dir_testinput = os.path.join(r'n:\Deltabox\Bulletin\veenstra\info dfm_tools\test_input')
    file_nc_map = os.path.join(dir_testinput,'DFM_sigma_curved_bend','DFM_OUTPUT_cb_3d','cb_3d_map.nc')
    file_nc_his = os.path.join(dir_testinput,'DFM_sigma_curved_bend','DFM_OUTPUT_cb_3d','cb_3d_his.nc')
    
    #get lists with vars/dims, times, station/crs/structures
    vars_pd, dims_pd = get_ncvardimlist(file_nc=file_nc_map)
    times_pd = get_timesfromnc(file_nc=file_nc_map)
    statlist_pd = get_hisstationlist(file_nc=file_nc_his, varname='station_name')
    
    #retrieve his data
    data_fromhis_wl = get_ncmodeldata(file_nc=file_nc_his, varname='waterlevel', station='all', timestep= 'all')
    fig, ax = plt.subplots(1,1,figsize=(10,5))
    for iP, station in enumerate(data_fromhis_wl.var_stations['station_name']):
        ax.plot(data_fromhis_wl.var_times,data_fromhis_wl[:,iP],'-', label=station)
    ax.legend()
    ax.set_ylabel('%s (%s)'%(data_fromhis_wl.var_varname, data_fromhis_wl.var_ncvarobject.units))
    
    #plot net/grid
    ugrid_all = get_netdata(file_nc=file_nc_map)#,multipart=False)
    fig, ax = plt.subplots()
    pc = plot_netmapdata(ugrid_all.verts, values=None, ax=None, linewidth=0.5, color="crimson", facecolor="None")
    ax.set_aspect('equal')

    #plot water level on map
    data_frommap_wl = get_ncmodeldata(file_nc=file_nc_map, varname='mesh2d_s1', timestep=3)#, multipart=False)
    fig, ax = plt.subplots()
    pc = plot_netmapdata(ugrid_all.verts, values=data_frommap_wl[0,:], ax=None, linewidth=0.5, cmap="jet")
    pc.set_clim([-0.5,1])
    fig.colorbar(pc, ax=ax)
    ax.set_title('%s (%s)'%(data_frommap_wl.var_varname, data_frommap_wl.var_ncvarobject.units))
    ax.set_aspect('equal')

    #plot salinity on map
    data_frommap_sal = get_ncmodeldata(file_nc=file_nc_map, varname='mesh2d_sa1', timestep=2, layer=5)#, multipart=False)
    fig, ax = plt.subplots()
    pc = plot_netmapdata(ugrid_all.verts, values=data_frommap_sal[0,:,0], ax=None, linewidth=0.5, cmap="jet")
    fig.colorbar(pc, ax=ax)
    ax.set_title('%s (%s)'%(data_frommap_sal.var_varname, data_frommap_sal.var_ncvarobject.units))
    ax.set_aspect('equal')
    
    #print contents of retrieved data withing data_frommap_sal variable
    print_var = data_frommap_sal
    print('++++++\nthe data in the variable %s is:\n%s\n'%(print_var.var_varname, print_var))
    print('++++++\nthe time indices and times in the variable %s are:\n%s\n'%(print_var.var_varname, print_var.var_times))
    print('++++++\nthe station indices and station names in the variable %s are:\n%s\n'%(print_var.var_varname, print_var.var_stations))
    print('++++++\nthe layer indices in the variable %s are:\n%s\n'%(print_var.var_varname, print_var.var_layers))
    print('++++++\nthe shape of the variable %s is:\n%s\n'%(print_var.var_varname, print_var.shape))
    print('++++++\nthe dimensions of the variable %s are (copied from netCDF variable):\n%s\n'%(print_var.var_varname, print_var.var_dimensions))
    print('++++++\nthe netCDF variable where the data in variable %s comes from is:\n%s\n'%(print_var.var_varname, print_var.var_ncvarobject))
    print('++++++\nsome example contents of this netCDF variable:')
    print('\tthe dimension names of the netCDF variable %s are:\n\t\t%s'%(print_var.var_varname, print_var.var_ncvarobject.dimensions))
    print('\tthe shape of the netCDF variable %s is:\n\t\t%s'%(print_var.var_varname, print_var.var_ncvarobject.shape))
    print('\tthe units of the netCDF variable %s are:\n\t\t%s'%(print_var.var_varname, print_var.var_ncvarobject.units))
    print('\tthe long_name of the netCDF variable %s is:\n\t\t%s'%(print_var.var_varname, print_var.var_ncvarobject.long_name))
    print('\tthe standard_name of the netCDF variable %s is:\n\t\t%s'%(print_var.var_varname, print_var.var_ncvarobject.standard_name))







