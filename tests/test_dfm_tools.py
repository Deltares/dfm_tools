#!/usr/bin/env python

"""Tests for dfm_tools package environment"""

import pytest
#import inspect
import os
import glob
from packaging import version

import dfm_tools as dfmt
import numpy as np
import datetime as dt
import pandas as pd
from netCDF4 import Dataset
import xarray as xr

dir_testinput = os.path.join(r'c:\DATA','dfm_tools_testdata')
dir_tests = os.path.dirname(__file__) #F9 doesnt work, only F5 (F5 also only method to reload external definition scripts)

# ACCEPTANCE TESTS VIA EXAMPLE SCRIPTS, these are the ones who are only meant to generate output files
list_configfiles = glob.glob(os.path.join(dir_tests,'examples','*.py')) + glob.glob(os.path.join(dir_tests,'examples_workinprogress','*.py'))
dir_output_general = os.path.join(dir_tests,'examples_output')
if not os.path.exists(dir_output_general):
    os.mkdir(dir_output_general)


@pytest.mark.acceptance
@pytest.mark.parametrize("file_config", [pytest.param(file_config, id=os.path.basename(file_config).replace('.py','')) for file_config in list_configfiles])
def test_run_examples(file_config):
    """
    file_config = os.path.join(dir_tests,'configfiles','predictie_2019_b02ex2_19Ycomp4Ydia_CUXHVN_test.py')
    """
    # 1. Set up test data
    dir_output = os.path.join(dir_output_general,os.path.basename(file_config).replace('.py',''))
    if not os.path.exists(dir_output):
        os.mkdir(dir_output)
    os.chdir(dir_output)
    #test = os.system("python {0} > {1}/FILE_DIAGNOSTICS.txt 2>&1".format(file_config, dir_output))#+ " & pause")
    test = os.system("python {0}".format(file_config))#+ " & pause")
    
    if test:
        raise Exception('execution did not finish properly')


##### UNITTESTS AND SYSTEMTESTS

modulename_list = ['os','sys','glob','shutil','scipy','numpy','datetime','pandas','matplotlib','netCDF4','click','shapely','shapely.geometry','cartopy','pyepsg'] #TODO: add xarray etc
@pytest.mark.parametrize("modulename", [pytest.param('%s'%(stat), id='%s'%(stat)) for stat in modulename_list])
@pytest.mark.unittest
def test_import_libraries(modulename):
    """
    tests whether shapely can be imported successfully, this is a problem in some environments
    in that case 'import shapely' works, but import 'shapely.geometry' fails
    """
    def try_importmodule(modulename=None):
        command = '\t- open command window (or anaconda prompt)\n\t- conda activate dfm_tools_env\n\t- conda install -c conda-forge %s'%(modulename)
    
        try:
            exec('import %s'%(modulename))
        except:
            raise Exception('ERROR: module %s not found, do the following:\n%s'%(modulename, command))
        if modulename == 'shapely':
            try:
                import shapely.geometry
            except:
                raise Exception('ERROR: cannot execute "import shapely.geometry", do the following:\n%s'%(command))
            if version.parse(shapely.__version__) < version.parse('1.7.0'):
                raise Exception(f'ERROR: incorrect shapely version ({shapely.__version__}), should be 1.7.0 or higher, do the following:\n{command}')
    
    try_importmodule(modulename=modulename)


@pytest.mark.parametrize("file_nc, varname, expected_size", [pytest.param(os.path.join(dir_testinput,r'DFM_3D_z_Grevelingen\computations\run01\DFM_OUTPUT_Grevelingen-FM\Grevelingen-FM_0*_map.nc'), 'mesh2d_sa1', (44796, 36), id='from partitioned map Grevelingen'),
                                                             pytest.param(os.path.join(r'p:\archivedprojects\11203850-coastserv\06-Model\waq_model\simulations\run0_20200319\DFM_OUTPUT_kzn_waq', 'kzn_waq_0*_map.nc'), 'Chlfa', (17385, 39), id='from partitioned waq map coastserv'),
                                                             #pytest.param(r'p:\11205258-006-kpp2020_rmm-g6\C_Work\01_Rooster\final_totaalmodel\rooster_rmm_v1p5_net.nc', 'mesh2d_face_x', (44804,), id='fromnet RMM'), #no time dimension
                                                             ])
@pytest.mark.unittest
def test_getmapdata(file_nc, varname, expected_size):
    """
    Checks whether ghost cells are properly taken care of (by asserting shape). And also whether varname can be found from attributes in case of Chlfa.
    
    file_nc = os.path.join(dir_testinput,'DFM_3D_z_Grevelingen','computations','run01','DFM_OUTPUT_Grevelingen-FM','Grevelingen-FM_0*_map.nc')
    expected_size = (44796,)
    """
    file_nc_nostar = file_nc.replace('0*','0000') #TODO: introduce support for * in dfm_tools definitions
    varname_found = dfmt.get_varnamefromattrs(file_nc_nostar,varname)
    
    data_xr_map = dfmt.open_partitioned_dataset(file_nc)
    data_varsel = data_xr_map[varname_found].isel(time=2)
    
    assert data_varsel.shape == expected_size


@pytest.mark.systemtest
def test_cartopy_epsg():
    
    #this one crashes if the dummy in plot_background() is not created
    dfmt.plot_background(ax=None, projection=28992, google_style='satellite', resolution=5, features='land', nticks=6, latlon_format=False, gridlines=False)


@pytest.mark.systemtest
def test_polygon_intersect(): #TODO: update to new xarray method

    file_nc = os.path.join(dir_testinput,'DFM_sigma_curved_bend\\DFM_OUTPUT_cb_3d\\cb_3d_map.nc')
    data_merged = dfmt.open_partitioned_dataset(file_nc)

    line_array = np.array([[2084.67741935, 3353.02419355], #with linebend in cell en with line crossing same cell twice
       [2255.79637097, 3307.15725806],
       [2222.27822581, 3206.60282258],
       [2128.78024194, 3266.58266129]])
    
    #intersect function, find crossed cell numbers (gridnos) and coordinates of intersection (2 per crossed cell)
    intersect_pd = dfmt.polygon_intersect(data_merged, line_array)
    intersect_pd = intersect_pd.sort_index()
    expected_intersectgridnos = np.array([ 91, 146, 146, 147, 147, 201, 201, 202], dtype=np.int64)
    expected_intersectcoords = np.array([[2084.67741935, 3353.02419355, 2144.15041424, 3337.08297842],
                                        [2144.15041424, 3337.08297842, 2202.53662217, 3321.43306702],
                                        [2173.05750857, 3238.17837704, 2128.78024194, 3266.58266129],
                                        [2202.53662217, 3321.43306702, 2255.79637097, 3307.15725806],
                                        [2255.79637097, 3307.15725806, 2246.9810802 , 3280.71138574],
                                        [2239.02015401, 3256.82860719, 2222.27822581, 3206.60282258],
                                        [2222.27822581, 3206.60282258, 2173.05750857, 3238.17837704],
                                        [2246.9810802 , 3280.71138574, 2239.02015401, 3256.82860719]])
    
    assert (intersect_pd.index-expected_intersectgridnos<1e-9).all()
    assert (intersect_pd.iloc[:,:4].values-expected_intersectcoords<1e-8).all()


@pytest.mark.unittest
def test_getvarnamemapnc():
    """
    this test tests if a netcdf varname can be retrieved from the 'dictionary' and if the variable can be retrieved from de netcdf
    """
    
    file_map = os.path.join(dir_testinput,r'DFM_3D_z_Grevelingen\computations\run01\DFM_OUTPUT_Grevelingen-FM\Grevelingen-FM_0000_map.nc')
    #file_net = os.path.join(dir_testinput,r'DFM_3D_z_Grevelingen\computations\run01\Grevelingen_FM_grid_20190603_net.nc')
    
    file_nc = file_map
    data_nc = Dataset(file_nc)
    varname_requested = 'NetNode_y' #is actually in file, so this is not a good test
    
    varname = dfmt.get_varname_fromnc(data_nc,varname_requested, vardim='var')
    data_nc_var = data_nc.variables[varname]
    dimname = data_nc_var.dimensions[0]
    
    # Verify expectations
    assert varname == 'mesh2d_node_y'
    assert dimname == 'nmesh2d_node'


#DEPRECATED TESTS
#DEPRECATED TESTS
#DEPRECATED TESTS
#DEPRECATED TESTS
#DEPRECATED TESTS
#DEPRECATED TESTS
#DEPRECATED TESTS
#DEPRECATED TESTS


@pytest.mark.parametrize("file_nc, expected_size", [pytest.param(os.path.join(dir_testinput,r'DFM_3D_z_Grevelingen\computations\run01\DFM_OUTPUT_Grevelingen-FM\Grevelingen-FM_0000_map.nc'), 5599, id='from 1 map partion Grevelingen'),
                                                    #pytest.param(r'p:\11205258-006-kpp2020_rmm-g6\C_Work\01_Rooster\final_totaalmodel\rooster_rmm_v1p5_net.nc', 44804?, id='fromnet RMM'),
                                                    pytest.param(os.path.join(dir_testinput,r'DFM_3D_z_Grevelingen\computations\run01\Grevelingen_FM_grid_20190603_net.nc'), 44804, id='fromnet Grevelingen')])
@pytest.mark.unittest
def SKIP_test_UGrid(file_nc, expected_size):
    
    ugrid = dfmt.UGrid.fromfile(file_nc)
    
    assert ugrid.verts.shape[0] == expected_size



@pytest.mark.parametrize("file_nc, expected_size", [pytest.param(os.path.join(dir_testinput,r'DFM_3D_z_Grevelingen\computations\run01\DFM_OUTPUT_Grevelingen-FM\Grevelingen-FM_0000_map.nc'), 44796, id='from partitioned map Grevelingen'),
                                                    #pytest.param(r'p:\11205258-006-kpp2020_rmm-g6\C_Work\01_Rooster\final_totaalmodel\rooster_rmm_v1p5_net.nc', 44804?, id='fromnet RMM'),
                                                    pytest.param(os.path.join(dir_testinput,r'DFM_3D_z_Grevelingen\computations\run01\Grevelingen_FM_grid_20190603_net.nc'), 44804, id='fromnet Grevelingen')])
@pytest.mark.unittest
def SKIP_test_getnetdata(file_nc, expected_size): #not relevant anymore since function will be phased out.
    """
    file_nc = os.path.join(dir_testinput,'DFM_3D_z_Grevelingen','computations','run01','DFM_OUTPUT_Grevelingen-FM','Grevelingen-FM_0000_map.nc')
    expected_size = 44796
    """
    
    ugrid = dfmt.get_netdata(file_nc)
    
    assert ugrid.verts.shape[0] == expected_size


@pytest.mark.parametrize("file_nc, varname, expected_size", [pytest.param(os.path.join(dir_testinput,'DFM_3D_z_Grevelingen','computations','run01','DFM_OUTPUT_Grevelingen-FM','Grevelingen-FM_0000_map.nc'), 'mesh2d_sa1', (1, 44796, 1), id='from partitioned map Grevelingen'),
                                                             pytest.param(os.path.join(r'p:\archivedprojects\11203850-coastserv\06-Model\waq_model\simulations\run0_20200319\DFM_OUTPUT_kzn_waq', 'kzn_waq_0000_map.nc'), 'Chlfa', (1, 17385, 1), id='from partitioned waq map coastserv')])
@pytest.mark.unittest
def SKIP_test_getmapbottomdata(file_nc, varname, expected_size): #skipping test since it is deprecated, retrieving bottom values is possible with ffill/bfill method like https://github.com/Deltares/dfm_tools/blob/main/tests/examples_workinprogress/workinprogress_exporttoshapefile.py#L50
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
    
    data_fromnc_bot = dfmt.get_ncmodeldata(file_nc=file_nc, varname=varname, timestep=3, layer='bottom')
    
    assert data_fromnc_bot.shape == expected_size
    assert data_fromnc_bot.mask.shape == expected_size


@pytest.mark.unittest
def SKIP_test_gethisbottomdata(): #skipping test since it is deprecated, retrieving bottom values is possible with ffill/bfill method like https://github.com/Deltares/dfm_tools/blob/main/tests/examples_workinprogress/workinprogress_exporttoshapefile.py#L50
    """
    checks whether time dimension is indexed once, resulted in shape of (10,10,1)
    """
    file_nc_his = os.path.join(dir_testinput,'MWRA','MB_02_0000_his.nc')
    
    moddata_dfm = dfmt.get_ncmodeldata(file_nc_his,varname='salinity',station='MWRA_F22',timestep=range(10),layer='bottom')
    
    assert moddata_dfm.shape == (10, 1, 1)


@pytest.mark.unittest
def SKIP_test_getncmodeldata_timeid(): #skipping test since it is deprecated
    
    file_map1 = os.path.join(dir_testinput,r'DFM_sigma_curved_bend\DFM_OUTPUT_cb_3d\cb_3d_map.nc')
    data_frommap = dfmt.get_ncmodeldata(file_nc=file_map1, varname='mesh2d_sa1', timestep=1, layer=5)#, multipart=False)
    
    assert (data_frommap.data[0,0,0] - 31. ) < 1E-9


@pytest.mark.unittest
def SKIP_test_getncmodeldata_indexcountmetadata(): #TODO: this is not valid nor necessary anymore when using xarray instead of dfm_tools for reading ncfiles (currently will crash because of station argument, so SKIPPING)
    
    #check if retrieving 1 index of data from 1 dimensional variable works (does not work if indices are np.arrays, so conversion to list in get_nc.py)
    file_his = os.path.join(dir_testinput,r'DFM_fou_RMM\RMM_dflowfm_0000_his.nc')
    data_statcoord = dfmt.get_ncmodeldata(file_nc=file_his, varname='station_x_coordinate',station='NM_1005.26_R_HBR-Cl_VCS-Nieuwe-Maas-85m')
    assert len(data_statcoord.var_stations) == 1
    
    file_his = os.path.join(dir_testinput,r'DFM_sigma_curved_bend\DFM_OUTPUT_cb_3d\cb_3d_his.nc')
    data_fromhis = dfmt.get_ncmodeldata(file_nc=file_his, varname='x_velocity', timestep=1, layer=5, station='Innersouth boundary')#, multipart=False)
    assert len(data_fromhis.var_times) == 1
    assert len(data_fromhis.var_layers) == 1
    assert len(data_fromhis.var_stations) == 1
    
    #data_fromhis_all = get_ncmodeldata(file_nc=file_his, varname='x_velocity', timestep='all', layer='all', station='all')
    data_fromhis = dfmt.get_ncmodeldata(file_nc=file_his, varname='x_velocity', timestep=[1,0,-5,1,4,3,2,0,1,-1], layer=[5,-2,3,0,5], station=[4,-1,0,0])
    assert len(data_fromhis.var_times) == 7
    assert data_fromhis.var_times.index.tolist() == [0,1,2,3,4,2156,2160]
    assert len(data_fromhis.var_layers) == 4
    assert data_fromhis.var_layers == [0,3,5,8]
    assert len(data_fromhis.var_stations) == 3
    assert data_fromhis.var_stations.index.tolist() == [0,4,5]


@pytest.mark.unittest
def SKIP_test_get_stationid_fromstationlist(): #skipping test since it is deprecated, xr.sel can be used if hisfile is opened with dfmt.preprocess_hisnc
    file_nc = os.path.join(dir_testinput,'DFM_3D_z_Grevelingen\\computations\\run01\\DFM_OUTPUT_Grevelingen-FM\\Grevelingen-FM_0000_his.nc')
    
    station = ['GTSO-02','GTSO-01','GTSO-03']
    
    data_xr = xr.open_dataset(file_nc)
    
    idx_stations = dfmt.get_stationid_fromstationlist(data_xr, stationlist=station)
    
    assert idx_stations==list([1,0,2])


@pytest.mark.unittest
def SKIP_test_getncmodeldata_datetime(): #skipping test since it is deprecated, datetime is easy with xr.sel()
    
    #retrieve all
    file_nc = os.path.join(dir_testinput,r'DFM_sigma_curved_bend\DFM_OUTPUT_cb_3d\cb_3d_map.nc')
    data_frommap = dfmt.get_ncmodeldata(file_nc=file_nc, varname='mesh2d_sa1', timestep='all', layer=5)#, multipart=False)
    assert data_frommap.shape[0] == len(data_frommap.var_times)
    
    #retrieve with numpy datetime array
    data_frommap = dfmt.get_ncmodeldata(file_nc=file_nc, varname='mesh2d_sa1', 
                                   timestep=np.arange(dt.datetime(2001,1,1),dt.datetime(2001,1,2),dt.timedelta(hours=1)), layer=5)#, multipart=False)
    assert (data_frommap.data[0,0,0] - 31. ) < 1E-9
    assert data_frommap.shape[0] == len(data_frommap.var_times)
    
    #retrieve with pandas date_range
    data_frommap = dfmt.get_ncmodeldata(file_nc=file_nc, varname='mesh2d_sa1', 
                                   timestep=pd.date_range(dt.datetime(2001,1,1),dt.datetime(2001,1,2),freq='30min'), layer=5)#, multipart=False)
    assert (data_frommap.data[0,0,0] - 31. ) < 1E-9
    assert data_frommap.shape[0] == len(data_frommap.var_times)


@pytest.mark.unittest
def SKIP_test_gettimesfromnc(): #skipping test since it is deprecated (time array extraction is easy+fast with xarray)
        
    #this file seems to have a reconstructable array, but contains gaps, so the time reader falls back to reading all
    file_nc = os.path.join(dir_testinput,'id1-LICHTELGRE.nc')
    data_fromhis_time = dfmt.get_ncmodeldata(file_nc,varname='time',timestep='all')


@pytest.mark.unittest
def SKIP_test_getncmatchingvarlist(): #skipping test since it is deprecated (use dfmt.get_varnamefromattrs(file_nc_nostar,varname) instead)
    """
    this test tests retrieves a pandas list of variable long names in a netcdf that match the pattern, useful for waq variables.
    """
    
    file_nc = os.path.join(dir_testinput,r'DFM_3D_z_Grevelingen\computations\run01\DFM_OUTPUT_Grevelingen-FM\Grevelingen-FM_0000_map.nc')
    vars_pd = dfmt.get_ncvarproperties(file_nc=file_nc)

    pattern = 'Flow .*component'
    vars_pd_matching = vars_pd[vars_pd.loc[:,'long_name'].str.match(pattern)] #does not have to stop after pattern
    #vars_pd_matching = vars_pd[vars_pd.loc[:,'long_name'].str.startswith('Flow') & vars_pd.loc[:,'long_name'].str.endswith('component')]
    varkeys_list_matching = list(vars_pd_matching.index)
    
    assert varkeys_list_matching == ['mesh2d_ucx', 'mesh2d_ucy', 'mesh2d_ucz', 'mesh2d_ucxa', 'mesh2d_ucya']


