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
   


def Test_grid_gethismodeldata(self):
    """
    this test retrieves his data#, and plots it
    """
    import matplotlib.pyplot as plt
    plt.close('all')
    
    from dfm_tools.grid import get_hismapmodeldata, plot_netmapdata
    
    file_his = r'c:\DATA\werkmap\vanJulien_shortmodelfiles\DFM_3D_z_Grevelingen\computations\run01\DFM_OUTPUT_Grevelingen-FM\Grevelingen-FM_0000_his.nc'
    
    #GREVELINGEN
    print('plot grid and values from mapdata (constantvalue, 1 dim)')
    data_fromhis = get_hismapmodeldata(file_nc=file_his, var_values='salinity', timestep=np.arange(0,-1), lay=5)#, multipart=False)
    fig, ax = plt.subplots()
    
    #pc = plot_netmapdata(ugrid.verts, values=data_frommap_flat, ax=None, linewidth=0.5, cmap="jet")
    #pc.set_clim([28,30.2])
    fig.colorbar(pc, ax=ax)
    ax.set_aspect('equal')
    
    
def Test_grid_getnetdata_getmapmodeldata_plotnetmapdata(self):
    """
    this test retrieves grid data, retrieves map data, and plots it
    """
    import matplotlib.pyplot as plt
    plt.close('all')
    
    from dfm_tools.grid import get_netdata, get_hismapmodeldata, plot_netmapdata
    
    file_map1 = r'c:\DATA\werkmap\vanJulien_shortmodelfiles\DFM_sigma_curved_bend\DFM_OUTPUT_cb_3d\cb_3d_map.nc'
    file_map8 = r'n:\My Documents\werkmap\vanJulien_shortmodelfiles\DFM_3D_z_Grevelingen\computations\run01\DFM_OUTPUT_Grevelingen-FM\Grevelingen-FM_0000_map.nc'
    file_map8 = r'c:\DATA\werkmap\vanJulien_shortmodelfiles\DFM_3D_z_Grevelingen\computations\run01\DFM_OUTPUT_Grevelingen-FM\Grevelingen-FM_0000_map.nc'
    file_map_rmm = r'p:\11205258-006-kpp2020_rmm-g6\C_Work\08_RMM_FMmodel\computations\run_156\DFM_OUTPUT_RMM_dflowfm\RMM_dflowfm_0000_map.nc'
    
    #CURVIBEND
    print('plot grid and values from mapdata (constantvalue, 1 dim)')
    ugrid = get_netdata(file_nc=file_map1)#,multipart=False)
    iT = 3 #for iT in range(10):
    data_frommap = get_hismapmodeldata(file_nc=file_map1, var_values='mesh2d_sa1', timestep='all', lay=5)#, multipart=False)
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
    data_frommap = get_hismapmodeldata(file_nc=file_map8, var_values='mesh2d_flowelem_bl', timestep=3, lay=33)#, multipart=False)
    data_frommap_flat = data_frommap.flatten()
    fig, ax = plt.subplots()
    pc = plot_netmapdata(ugrid_all.verts, values=data_frommap_flat, ax=None, linewidth=0.5, cmap="jet")
    #pc.set_clim([28,30.2])
    fig.colorbar(pc, ax=ax)
    ax.set_aspect('equal')

    print('plot grid and values from mapdata (waterlevel, 2dim)')
    data_frommap = get_hismapmodeldata(file_nc=file_map8, var_values='mesh2d_s1', timestep=3, lay=33)#, multipart=False)
    data_frommap_flat = data_frommap.flatten()
    fig, ax = plt.subplots()
    pc = plot_netmapdata(ugrid_all.verts, values=data_frommap_flat, ax=None, linewidth=0.5, cmap="jet")
    pc.set_clim([-0.5,1])
    fig.colorbar(pc, ax=ax)
    ax.set_aspect('equal')

    print('plot grid and values from mapdata (salinity on layer, 3dim)')
    data_frommap = get_hismapmodeldata(file_nc=file_map8, var_values='mesh2d_sa1', timestep=3, lay=33)#, multipart=False)
    data_frommap_flat = data_frommap.flatten()
    fig, ax = plt.subplots()
    pc = plot_netmapdata(ugrid_all.verts, values=data_frommap_flat, ax=None, linewidth=0.5, cmap="jet")
    pc.set_clim([28,30.2])
    fig.colorbar(pc, ax=ax)
    ax.set_aspect('equal')

    print('plot grid and values from mapdata (temperature on layer, 3dim)')
    data_frommap = get_hismapmodeldata(file_nc=file_map8, var_values='mesh2d_tem1', timestep=3, lay=33)#, multipart=False)
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
    #data_frommap = get_hismapmodeldata(file_nc=file_map_rmm, var_values='mesh2d_s1', timestep=50, lay=0)#, multipart=False)
    data_frommap = get_hismapmodeldata(file_nc=file_map_rmm, var_values='mesh2d_ucx', timestep=50, lay=0)#, multipart=False)
    data_frommap_flat = data_frommap.flatten()
    fig, ax = plt.subplots()
    pc = plot_netmapdata(ugrid_all.verts, values=data_frommap_flat, ax=None, linewidth=0.5, cmap="jet")
    pc.set_clim([-1,1])
    fig.colorbar(pc, ax=ax)
    ax.set_aspect('equal')



def Test_maplora(self):
    from netCDF4 import Dataset
    import matplotlib.pyplot as plt
    plt.close('all')
    
    from dfm_tools.grid import get_netdata, get_hismapmodeldata, plot_netmapdata

    file_maplora = r'p:\11201806-sophie\Oosterschelde\WAQ\r02\postprocessing\oost_tracer_2_map.nc'
    
    ugrid_lora = get_netdata(file_maplora)

    data_nc = Dataset(file_maplora)
    list(data_nc.variables.keys())


    
    var_names = ['mesh2d_FColi','mesh2d_HIWAI','mesh2d_Pharma']
    var_clims = [None,[0,100000000000],[0,10000]]
    for var_name, var_clim in zip(var_names, var_clims):
        data_fromnc = data_nc.variables[var_name]
        
        fig, ax = plt.subplots()
        pc = plot_netmapdata(ugrid_lora.verts, values=data_fromnc, ax=None, linewidth=0.5, cmap="jet")
        if var_clim != None:
            pc.set_clim(var_clim)
        fig.colorbar(pc, ax=ax)
        ax.set_aspect('equal')
        ax.set_xlabel(var_name)

    #TODO: check if this now also works
    print('plot grid and values from mapdata (constantvalue, 1 dim)')
    data_frommap = get_hismapmodeldata(file_nc=file_maplora, var_values='mesh2d_FColi', timestep=3, lay=33)#, multipart=False)
    data_frommap_flat = data_frommap.flatten()
    fig, ax = plt.subplots()
    pc = plot_netmapdata(ugrid_lora.verts, values=data_frommap_flat, ax=None, linewidth=0.5, cmap="jet")
    #pc.set_clim([28,30.2])
    fig.colorbar(pc, ax=ax)
    ax.set_aspect('equal')
    

    
    
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
