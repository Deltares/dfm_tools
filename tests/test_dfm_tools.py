#!/usr/bin/env python

"""Tests for `dfm_tools` package."""

import pytest


def Test_mdu(self):
    #from netCDF4 import Dataset
    from dfm_tools.mdu import read_deltares_ini, write_deltares_ini
    
    filename_mdu = r'n:\My Documents\werkmap\vanJulien_shortmodelfiles\DFM_3D_z_Grevelingen\computations\run01\Grevelingen-FM.mdu'
    filename_mdu_out = r'n:\My Documents\werkmap\vanJulien_shortmodelfiles\DFM_3D_z_Grevelingen\computations\run01\Grevelingen-FM_out.mdu'
    
    data_mdu = read_deltares_ini(filename_mdu)
    write_deltares_ini(data_mdu, filename_mdu_out)
    
    
    assert 1==1


def Test_getvarnamemapnc(self):
    from netCDF4 import Dataset
    
    from dfm_tools.get_varname_mapnc import get_varname_mapnc
    
    # 1. define test data
    file_map = r'n:\My Documents\werkmap\vanJulien_shortmodelfiles\DFM_3D_z_Grevelingen\computations\run01\DFM_OUTPUT_Grevelingen-FM\Grevelingen-FM_0000_map.nc'
    #file_net = r'n:\My Documents\werkmap\vanJulien_shortmodelfiles\DFM_3D_z_Grevelingen\computations\run01\Grevelingen_FM_grid_20190603_net.nc'
    
    file_nc = file_map
    data_nc = Dataset(file_nc)
    varname_requested = 'NetNode_y' #is actually in file, so this is not a good test
    
    
    # 2. define initial expectations
    
    
    # 3. run test
    varname = get_varname_mapnc(data_nc,varname_requested)
    data_nc_var = data_nc.variables[varname]
    print(data_nc_var)
    
    # 4. Vefiry final expectations
    assert 1==1


def Test_grid_UGrid(self):
    #from netCDF4 import Dataset
    from dfm_tools.grid import UGrid
    
    #file_map = r'n:\My Documents\werkmap\vanJulien_shortmodelfiles\DFM_3D_z_Grevelingen\computations\run01\DFM_OUTPUT_Grevelingen-FM\Grevelingen-FM_0000_map.nc'
    file_net = r'n:\My Documents\werkmap\vanJulien_shortmodelfiles\DFM_3D_z_Grevelingen\computations\run01\Grevelingen_FM_grid_20190603_net.nc'
    #file_net = r'p:\11205258-006-kpp2020_rmm-g6\C_Work\01_Rooster\final_totaalmodel\rooster_rmm_v1p5_net.nc'
    
    data_ncUG = UGrid.fromfile(file_net)
    test = data_ncUG.cellcoords()
    print(test)
    assert 1==1


    
    
    
def Test_grid_getnetdata_plotnet(self):
    """
    this test retrieves grid data and plots it
    """
    import matplotlib.pyplot as plt
    plt.close('all')

    from dfm_tools.grid import get_netdata, plot_netmapdata

    file_net = r'n:\My Documents\werkmap\vanJulien_shortmodelfiles\DFM_3D_z_Grevelingen\computations\run01\Grevelingen_FM_grid_20190603_net.nc'
    file_net = r'c:\DATA\werkmap\vanJulien_shortmodelfiles\DFM_3D_z_Grevelingen\computations\run01\Grevelingen_FM_grid_20190603_net.nc'
    file_net_rmm = r'p:\11205258-006-kpp2020_rmm-g6\C_Work\01_Rooster\final_totaalmodel\rooster_rmm_v1p5_net.nc'

    print('plot only grid from net.nc')
    ugrid = get_netdata(file_nc=file_net)
    fig, ax = plt.subplots()
    plot_netmapdata(ugrid.verts, values=None, ax=None, linewidth=0.5, color="crimson", facecolor="None")
    ax.set_aspect('equal')

    print('plot only grid from net.nc (RMM)')
    ugrid = get_netdata(file_net_rmm)
    fig, ax = plt.subplots()
    plot_netmapdata(ugrid.verts, values=None, ax=None, linewidth=0.5, color="crimson", facecolor="None")
    ax.set_aspect('equal')
 
    

def Test_grid_gethismodeldata_Nithin(self):
    """
    this test retrieves his data#, and plots it
    """
    import matplotlib.pyplot as plt
    plt.close('all')
    
    from dfm_tools.grid import get_netdata, plot_netmapdata, get_hismapmodeldata
    
    file_net = r'c:\DATA\werkmap\vanNithin_shortmodelfiles\myortho3_net.nc'
    #file_net = r'n:\My Documents\werkmap\vanNithin_shortmodelfiles\myortho3_withcellinfo_net.nc'
    file_net = r'c:\DATA\werkmap\vanNithin_shortmodelfiles\myortho3_RGFGRID_net.nc'
    file_his = r'c:\DATA\werkmap\vanNithin_shortmodelfiles\tttz_0000_his.nc'
    
    station = ['Peiraias', 'Ovrios_2','Ovrios','Ovrios']
    #station = ['Peiraias']

    print('plot only grid from net.nc')
    ugrid = get_netdata(file_nc=file_net)
    fig, ax = plt.subplots()
    plot_netmapdata(ugrid.verts, values=None, ax=None, linewidth=0.5, color="crimson", facecolor="None")
    ax.set_aspect('equal')
    
    #NITHIN
    data_fromhis = get_hismapmodeldata(file_nc=file_his, varname='bedlevel', station=station)#, multipart=False)
    fig, ax = plt.subplots()
    ax.plot(data_fromhis.var_stations,data_fromhis,'-')
    ax.tick_params('x',rotation=30)

    data_fromhis = get_hismapmodeldata(file_nc=file_his, varname='waterlevel', timestep='all', station=station)#, multipart=False)
    fig, ax = plt.subplots()
    ax.plot(data_fromhis.var_times,data_fromhis,'-')
    
    data_fromhis = get_hismapmodeldata(file_nc=file_his, varname='salinity', timestep='all', lay=5, station=station)#, multipart=False)
    data_fromhis_flat = data_fromhis[:,:,0]
    fig, ax = plt.subplots()
    ax.plot(data_fromhis.var_times,data_fromhis_flat,'-')
    
    

def Test_grid_gethismodeldata(self):
    """
    this test retrieves his data#, and plots it
    """
    import matplotlib.pyplot as plt
    plt.close('all')
    
    from dfm_tools.grid import get_hismapmodeldata
    
    file_his = r'c:\DATA\werkmap\vanJulien_shortmodelfiles\DFM_3D_z_Grevelingen\computations\run01\DFM_OUTPUT_Grevelingen-FM\Grevelingen-FM_0000_his.nc'
    
    #GREVELINGEN
    print('plot bedlevel from his')
    data_fromhis = get_hismapmodeldata(file_nc=file_his, varname='bedlevel', station='all')#, multipart=False)
    fig, ax = plt.subplots()
    ax.plot(data_fromhis.var_stations,data_fromhis,'-')
    ax.tick_params('x',rotation=30)

    print('plot waterlevel from his')
    data_fromhis = get_hismapmodeldata(file_nc=file_his, varname='waterlevel', timestep='all', station='all')#, multipart=False)
    fig, ax = plt.subplots()
    ax.plot(data_fromhis.var_times,data_fromhis,'-')
    
    print('plot salinity from his')
    data_fromhis = get_hismapmodeldata(file_nc=file_his, varname='salinity', timestep='all', lay=5, station='all')#, multipart=False)
    data_fromhis_flat = data_fromhis[:,:,0]
    fig, ax = plt.subplots()
    ax.plot(data_fromhis.var_times,data_fromhis_flat,'-')

    print('plot salinity,bedlevel')
    #depth retrieval is probably wrong
    data_fromhis_depth = get_hismapmodeldata(file_nc=file_his, varname='zcoordinate_c', timestep=4, lay='all', station='Bommenede')#, multipart=False)
    data_fromhis_depth_flat = data_fromhis_depth[0,0,:]
    data_fromhis = get_hismapmodeldata(file_nc=file_his, varname='salinity', timestep=4, lay='all', station='Bommenede')#, multipart=False)
    data_fromhis_flat = data_fromhis[0,0,:]
    fig, ax = plt.subplots()
    ax.plot(data_fromhis_flat, data_fromhis_depth_flat,'-')
    
    
    
def Test_grid_getnetdata_getmapmodeldata_plotnetmapdata(self):
    """
    this test retrieves grid data, retrieves map data, and plots it
    """
    import matplotlib.pyplot as plt
    plt.close('all')
    import datetime as dt
    import numpy as np
    
    from dfm_tools.grid import get_netdata, get_hismapmodeldata, plot_netmapdata
    
    file_map1 = r'c:\DATA\werkmap\vanJulien_shortmodelfiles\DFM_sigma_curved_bend\DFM_OUTPUT_cb_3d\cb_3d_map.nc'
    file_map8 = r'n:\My Documents\werkmap\vanJulien_shortmodelfiles\DFM_3D_z_Grevelingen\computations\run01\DFM_OUTPUT_Grevelingen-FM\Grevelingen-FM_0000_map.nc'
    file_map8 = r'c:\DATA\werkmap\vanJulien_shortmodelfiles\DFM_3D_z_Grevelingen\computations\run01\DFM_OUTPUT_Grevelingen-FM\Grevelingen-FM_0000_map.nc'
    file_map_rmm = r'p:\11205258-006-kpp2020_rmm-g6\C_Work\08_RMM_FMmodel\computations\run_156\DFM_OUTPUT_RMM_dflowfm\RMM_dflowfm_0000_map.nc'
    
    #CURVIBEND (datetime)
    print('plot grid and values from mapdata (constantvalue, 1 dim)')
    ugrid = get_netdata(file_nc=file_map1)#,multipart=False)
    #iT = 3 #for iT in range(10):
    data_frommap = get_hismapmodeldata(file_nc=file_map1, varname='mesh2d_sa1', timestep=np.arange(dt.datetime(2001,1,1),dt.datetime(2001,1,2),dt.timedelta(hours=1)), lay=5)#, multipart=False)
    data_frommap_flat = data_frommap[0,:,0]
    #data_frommap_depth = get_hismapmodeldata(file_nc=file_map1, varname='mesh2d_layer_sigma', lay='all')#, multipart=False)
    fig, ax = plt.subplots()
    pc = plot_netmapdata(ugrid.verts, values=data_frommap_flat, ax=None, linewidth=0.5, cmap="jet")
    #pc.set_clim([28,30.2])
    fig.colorbar(pc, ax=ax)
    ax.set_aspect('equal')
    
    #CURVIBEND
    print('plot grid and values from mapdata (constantvalue, 1 dim)')
    ugrid = get_netdata(file_nc=file_map1)#,multipart=False)
    #iT = 3 #for iT in range(10):
    data_frommap = get_hismapmodeldata(file_nc=file_map1, varname='mesh2d_sa1', timestep=3, lay=5)#, multipart=False)
    data_frommap_flat = data_frommap.flatten()
    fig, ax = plt.subplots()
    pc = plot_netmapdata(ugrid.verts, values=data_frommap_flat, ax=None, linewidth=0.5, cmap="jet")
    #pc.set_clim([28,30.2])
    fig.colorbar(pc, ax=ax)
    ax.set_aspect('equal')
        
    #GREVELINGEN
    print('plot only grid from mapdata')
    ugrid_all = get_netdata(file_nc=file_map8)#,multipart=False)
    fig, ax = plt.subplots()
    pc = plot_netmapdata(ugrid_all.verts, values=None, ax=None, linewidth=0.5, color="crimson", facecolor="None")
    ax.set_aspect('equal')
    
    print('plot grid and values from mapdata (constantvalue, 1 dim)')
    data_frommap = get_hismapmodeldata(file_nc=file_map8, varname='mesh2d_flowelem_bl')#, multipart=False)
    data_frommap_flat = data_frommap.flatten()
    fig, ax = plt.subplots()
    pc = plot_netmapdata(ugrid_all.verts, values=data_frommap_flat, ax=None, linewidth=0.5, cmap="jet")
    #pc.set_clim([28,30.2])
    fig.colorbar(pc, ax=ax)
    ax.set_aspect('equal')

    print('plot grid and values from mapdata (waterlevel, 2dim)')
    data_frommap = get_hismapmodeldata(file_nc=file_map8, varname='mesh2d_s1', timestep=3)#, multipart=False)
    data_frommap_flat = data_frommap.flatten()
    fig, ax = plt.subplots()
    pc = plot_netmapdata(ugrid_all.verts, values=data_frommap_flat, ax=None, linewidth=0.5, cmap="jet")
    pc.set_clim([-0.5,1])
    fig.colorbar(pc, ax=ax)
    ax.set_aspect('equal')

    print('plot grid and values from mapdata (salinity on layer, 3dim)')
    data_frommap = get_hismapmodeldata(file_nc=file_map8, varname='mesh2d_sa1', timestep=3, lay=33)#, multipart=False)
    data_frommap_flat = data_frommap.flatten()
    fig, ax = plt.subplots()
    pc = plot_netmapdata(ugrid_all.verts, values=data_frommap_flat, ax=None, linewidth=0.5, cmap="jet")
    pc.set_clim([28,30.2])
    fig.colorbar(pc, ax=ax)
    ax.set_aspect('equal')

    print('plot grid and values from mapdata (temperature on layer, 3dim)')
    data_frommap = get_hismapmodeldata(file_nc=file_map8, varname='mesh2d_tem1', timestep=3, lay=33)#, multipart=False)
    data_frommap_flat = data_frommap.flatten()
    fig, ax = plt.subplots()
    pc = plot_netmapdata(ugrid_all.verts, values=data_frommap_flat, ax=None, linewidth=0.5, cmap="jet")
    pc.set_clim([4,10])
    fig.colorbar(pc, ax=ax)
    ax.set_aspect('equal')

    #RMM
    print('plot only grid from mapdata (RMM)')
    ugrid_all = get_netdata(file_nc=file_map_rmm)#,multipart=False)
    fig, ax = plt.subplots()
    pc = plot_netmapdata(ugrid_all.verts, values=None, ax=None, linewidth=0.5, color="crimson", facecolor="None")
    ax.set_aspect('equal')

    print('plot grid and values from mapdata (RMM)')
    ugrid_all = get_netdata(file_nc=file_map_rmm)#,multipart=False)
    #data_frommap = get_hismapmodeldata(file_nc=file_map_rmm, varname='mesh2d_s1', timestep=50)#, multipart=False)
    data_frommap = get_hismapmodeldata(file_nc=file_map_rmm, varname='mesh2d_ucx', timestep=50)#, multipart=False)
    data_frommap_flat = data_frommap.flatten()
    fig, ax = plt.subplots()
    pc = plot_netmapdata(ugrid_all.verts, values=data_frommap_flat, ax=None, linewidth=0.5, cmap="jet")
    pc.set_clim([-1,1])
    fig.colorbar(pc, ax=ax)
    ax.set_aspect('equal')


def Test_maplora(self):
    import matplotlib.pyplot as plt
    plt.close('all')
    
    from dfm_tools.grid import get_netdata, plot_netmapdata, get_hismapmodeldata

    file_maplora = r'p:\11201806-sophie\Oosterschelde\WAQ\r03\postprocessing\oost_tracer_map.nc'
    file_maplora = r'p:\11201806-sophie\Oosterschelde\WAQ\r02\postprocessing\oost_tracer_2_map.nc'
    
    ugrid_lora = get_netdata(file_nc=file_maplora)

    print('plot grid and values from mapdata (constantvalue, 1 dim)')
    var_names = ['mesh2d_FColi','mesh2d_HIWAI','mesh2d_mspaf','mesh2d_Pharma']
    var_clims = [None,[0,100000000000],None,[0,10000]]
    var_names = ['mesh2d_FColi','mesh2d_HIWAI','mesh2d_Pharma']
    var_clims = [None,[0,100000000000],[0,10000]]
    #var_names = ['mesh2d_Pharma']
    #var_clims = [[0,10000]]
    for var_name, var_clim in zip(var_names, var_clims):
        fig, ax = plt.subplots()
        if 'oost_tracer_2_map' in file_maplora:
            data_frommap = get_hismapmodeldata(file_nc=file_maplora, varname=var_name)#, multipart=False)
        else:
            data_frommap = get_hismapmodeldata(file_nc=file_maplora, varname=var_name, timestep='all', lay=5)#, multipart=False)
            data_frommap = data_frommap.flatten()
        pc = plot_netmapdata(ugrid_lora.verts, values=data_frommap, ax=None, linewidth=0.5, cmap="jet")
        if var_clim != None:
            pc.set_clim(var_clim)
        fig.colorbar(pc, ax=ax)
        ax.set_aspect('equal')
        ax.set_xlabel(var_name)




def Test_grid_get_modeldata_onintersection(self):
    import matplotlib.pyplot as plt
    plt.close('all')
    import numpy as np
    
    from dfm_tools.grid import get_netdata, get_hismapmodeldata, get_modeldata_onintersection, plot_netmapdata
    from dfm_tools.polygon import LineBuilder#, Polygon
    
    file_map = r'c:\DATA\werkmap\vanJulien_shortmodelfiles\DFM_sigma_curved_bend\DFM_OUTPUT_cb_3d\cb_3d_map.nc'
    file_map = r'c:\DATA\werkmap\vanJulien_shortmodelfiles\DFM_3D_z_Grevelingen\computations\run01\DFM_OUTPUT_Grevelingen-FM\Grevelingen-FM_0000_map.nc'
    #file_map = r'p:\11203379-mwra-new-bem-model\waq_model\simulations\A31_1year_20191219\DFM_OUTPUT_MB_02_waq\MB_02_waq_0000_map.nc'
    
    if 'cb_3d_map' in file_map:
        timestep = 72
        layno = 5
        convert2merc = None
        line_array = np.array([[ 185.08667065, 2461.11775254],
                               [2934.63837418, 1134.16019127]])
        line_array = np.array([[ 104.15421399, 2042.7077107 ],
                               [2913.47878063, 2102.48057382]])
        clim = None
    elif 'Grevelingen' in file_map:
        timestep = 3
        layno = 35
        convert2merc = None
        line_array = np.array([[ 56267.59146475, 415644.67447155],
                               [ 64053.73427496, 419407.58239502]])
        line_array = np.array([[ 53181.96942503, 424270.83361629],
                               [ 55160.15232593, 416913.77136685]])
        line_array = np.array([[ 52787.21854294, 424392.10414528],
                               [ 55017.72655174, 416403.77313703],
                               [ 65288.43784807, 419360.49305567]])
        clim = [-25,5]
    elif 'DFM_OUTPUT_MB_02_waq' in file_map:
        timestep = 30
        layno = 5
        convert2merc = True
        line_array = np.array([[42.38, 42.38],
                               [-70.98, -70.25]])
        clim = None
    else:
        raise Exception('ERROR: no settings provided for this mapfile')
    
    
    ugrid = get_netdata(file_nc=file_map)#,multipart=False)
    
    #create plot with ugrid and cross section line
    fig, ax_input = plt.subplots()
    pc = plot_netmapdata(ugrid.verts, values=None, ax=ax_input, linewidth=0.5, color='crimson', facecolor="None")
    #pc.set_clim([28,30.2])
    #fig.colorbar(pc, ax=ax)
    ax_input.set_aspect('equal')
    if 0: #click interactive polygon
        #pol_frominput = Polygon.frominteractive(ax)
        line, = ax_input.plot([], [],'o-')  # empty line
        linebuilder = LineBuilder(line)
        line_array = linebuilder.line_array
    ax_input.plot(line_array[:,0],line_array[:,1])
    
    
    #intersect function, find crossed cell numbers (gridnos) and coordinates of intersection (2 per crossed cell)
    intersect_gridnos, intersect_coords = ugrid.polygon_intersect(line_array)
    
    #derive vertices from cross section (distance from first point)
    crs_verts = get_modeldata_onintersection(file_nc=file_map, line_array=line_array, intersect_gridnos=intersect_gridnos, intersect_coords=intersect_coords, timestep=timestep, convert2merc=convert2merc)
    
    #get bed layer
    data_frommap_bl = get_hismapmodeldata(file_nc=file_map, varname='mesh2d_flowelem_bl')#, multipart=False)
    pc = plot_netmapdata(ugrid.verts, values=data_frommap_bl, ax=ax_input, linewidth=0.5, cmap="jet")

    #get data to plot and select data for cross section (gridnos)
    data_frommap = get_hismapmodeldata(file_nc=file_map, varname='mesh2d_sa1', timestep=timestep, lay='all')#, multipart=False)
    data_frommap_sel = data_frommap[0,intersect_gridnos,:]
    data_frommap_sel_flat = data_frommap_sel.T.flatten()
    
    #plot crossed cells (gridnos) in first plot
    data_frommap_flat = data_frommap[0,intersect_gridnos,layno]
    pc = plot_netmapdata(ugrid.verts[intersect_gridnos,:,:], values=data_frommap_flat, ax=ax_input, linewidth=0.5, cmap="jet")
    
    #plot cross section
    fig, ax = plt.subplots()
    pc = plot_netmapdata(crs_verts, values=data_frommap_sel_flat, ax=ax, linewidth=0.5, cmap='jet')
    fig.colorbar(pc, ax=ax)
    ax.set_ylim(clim)
    



    
@pytest.fixture
def response():
    """Sample pytest fixture.

    See more at: http://doc.pytest.org/en/latest/fixture.html
    """
    # import requests
    # return requests.get('https://github.com/audreyr/cookiecutter-pypackage')


def test_content(response):
    """Sample pytest test function with the pytest fixture as an argument."""
    # from bs4 import BeautifulSoup
    # assert 'GitHub' in BeautifulSoup(response.content).title.string
