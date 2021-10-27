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
    """
    file_nc = os.path.join(dir_testinput,'DFM_3D_z_Grevelingen','computations','run01','DFM_OUTPUT_Grevelingen-FM','Grevelingen-FM_0000_map.nc')
    expected_size = 44796
    """
    from dfm_tools.get_nc import get_netdata
    
    ugrid = get_netdata(file_nc)
    
    assert ugrid.verts.shape[0] == expected_size


@pytest.mark.parametrize("file_nc, expected_size", [pytest.param(os.path.join(dir_testinput,'DFM_3D_z_Grevelingen','computations','run01','DFM_OUTPUT_Grevelingen-FM','Grevelingen-FM_0000_map.nc'), (1,44796), id='from partitioned map Grevelingen')])
@pytest.mark.unittest
def test_getmapdata(file_nc, expected_size):
    """
    Checks whether ghost cells are properly taken care of.
    
    file_nc = os.path.join(dir_testinput,'DFM_3D_z_Grevelingen','computations','run01','DFM_OUTPUT_Grevelingen-FM','Grevelingen-FM_0000_map.nc')
    expected_size = (1, 44796)
    """
    from dfm_tools.get_nc import get_ncmodeldata
    
    data_nc = get_ncmodeldata(file_nc=file_nc, varname='mesh2d_s1',timestep=2)
    
    assert data_nc.shape == expected_size


@pytest.mark.parametrize("file_nc, varname, expected_size", [pytest.param(os.path.join(dir_testinput,'DFM_3D_z_Grevelingen','computations','run01','DFM_OUTPUT_Grevelingen-FM','Grevelingen-FM_0000_map.nc'), 'mesh2d_sa1', (1, 44796, 1), id='from partitioned map Grevelingen'),
                                                             pytest.param(os.path.join(r'p:\11203850-coastserv\06-Model\waq_model\simulations\run0_20200319\DFM_OUTPUT_kzn_waq', 'kzn_waq_0000_map.nc'), 'Chlfa', (1, 17385, 1), id='from partitioned waq map coastserv')])
@pytest.mark.unittest
def test_getmapbottomdata(file_nc, varname, expected_size):
    """
    Checks whether ghost cells are properly taken care and if the mask comes trough (needed for selecting the bottom layer).
    It will fail on get_ncmodeldata, the assertions are just extras
    
    file_nc = os.path.join(dir_testinput,'DFM_3D_z_Grevelingen','computations','run01','DFM_OUTPUT_Grevelingen-FM','Grevelingen-FM_0000_map.nc')
    varname = 'mesh2d_sa1'
    expected_size = (1, 44796, 1)
    
    
    file_nc = os.path.join('p:\\11203850-coastserv\\06-Model\\waq_model\\simulations\\run0_20200319\\DFM_OUTPUT_kzn_waq', 'kzn_waq_0000_map.nc')
    varname = 'Chlfa'
    expected_size = (1, 17385, 1)
    """
    from dfm_tools.get_nc import get_ncmodeldata
    
    data_fromnc_bot = get_ncmodeldata(file_nc=file_nc, varname=varname, timestep=3, layer='bottom')
    
    assert data_fromnc_bot.shape == expected_size
    assert data_fromnc_bot.mask.shape == expected_size


@pytest.mark.unittest
def test_gethisbottomdata():
    """
    checks whether time dimension is indexed once, resulted in shape of (10,10,1)
    """
    from dfm_tools.get_nc import get_ncmodeldata
    file_nc_his = os.path.join(dir_testinput,'MWRA','MB_02_0000_his.nc')
    
    moddata_dfm = get_ncmodeldata(file_nc_his,varname='salinity',station='MWRA_F22',timestep=range(10),layer='bottom')
    
    assert moddata_dfm.shape == (10, 1, 1)


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
    file_his = os.path.join(dir_testinput,r'DFM_fou_RMM\RMM_dflowfm_0000_his.nc')
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


@pytest.mark.systemtest
def test_cartopy_epsg():
    
    from dfm_tools.testutils import try_importmodule
    try_importmodule(modulename='cartopy') #check if cartopy was installed since it is an optional module, also happens in plot_cartopybasemap()
    
    from dfm_tools.get_nc import plot_background
    
    #this one crashes if the dummy in plot_background() is not created
    plot_background(ax=None, projection=28992, google_style='satellite', resolution=5, features='land', nticks=6, latlon_format=False, gridlines=False)


@pytest.mark.systemtest
def test_UGrid_polygon_intersect():
    import numpy as np
    from dfm_tools.get_nc import get_netdata
    
    file_nc = os.path.join(dir_testinput,'DFM_sigma_curved_bend\\DFM_OUTPUT_cb_3d\\cb_3d_map.nc')
    multipart = None
    line_array = np.array([[2084.67741935, 3353.02419355], #with linebend in cell en with line crossing same cell twice
       [2255.79637097, 3307.15725806],
       [2222.27822581, 3206.60282258],
       [2128.78024194, 3266.58266129]])

    ugrid = get_netdata(file_nc=file_nc, multipart=multipart)

    #intersect function, find crossed cell numbers (gridnos) and coordinates of intersection (2 per crossed cell)
    intersect_gridnos, intersect_coords = ugrid.polygon_intersect(line_array, optimize_dist=False)
    
    expected_intersectgridnos = np.array([ 91, 146, 146, 147, 147, 201, 201, 202], dtype=np.int64)
    expected_intersectcoords = np.array([[[2084.67741935, 2144.15041424],
                                            [3353.02419355, 3337.08297842]],
                                    
                                           [[2144.15041424, 2202.53662217],
                                            [3337.08297842, 3321.43306702]],
                                    
                                           [[2173.05750857, 2128.78024194],
                                            [3238.17837704, 3266.58266129]],
                                    
                                           [[2202.53662217, 2255.79637097],
                                            [3321.43306702, 3307.15725806]],
                                    
                                           [[2255.79637097, 2246.9810802 ],
                                            [3307.15725806, 3280.71138574]],
                                    
                                           [[2239.02015401, 2222.27822581],
                                            [3256.82860719, 3206.60282258]],
                                    
                                           [[2222.27822581, 2173.05750857],
                                            [3206.60282258, 3238.17837704]],
                                    
                                           [[2246.9810802 , 2239.02015401],
                                            [3280.71138574, 3256.82860719]]])
    
    assert (intersect_gridnos-expected_intersectgridnos<1e-9).all()
    assert (intersect_coords-expected_intersectcoords<1e-8).all()



