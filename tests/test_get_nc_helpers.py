# -*- coding: utf-8 -*-
"""
Created on Fri Mar 20 08:57:42 2020

@author: veenstra
"""

import pytest
#import inspect
import os

dir_testinput = os.path.join(r'c:/DATA','dfm_tools_testdata')
#from tests.TestUtils import getmakeoutputdir



@pytest.mark.unittest    
def test_gettimesfromnc():

    from dfm_tools.get_nc_helpers import get_timesfromnc

    #file_nc = os.path.join(dir_testinput,r'DFM_3D_z_Grevelingen\computations\run01\DFM_OUTPUT_Grevelingen-FM\Grevelingen-FM_0000_map.nc')
    file_nc = r'p:\11205258-006-kpp2020_rmm-g6\C_Work\08_RMM_FMmodel\computations\run_156\DFM_OUTPUT_RMM_dflowfm\RMM_dflowfm_0000_map.nc'
    
    times_default = get_timesfromnc(file_nc=file_nc)
    times_noreconstruct = get_timesfromnc(file_nc=file_nc, force_noreconstruct=True)
    
    assert len(times_default) == len(times_noreconstruct)
    assert times_default.iloc[0] == times_noreconstruct.iloc[0]
    assert times_default.iloc[1] == times_noreconstruct.iloc[1]
    assert times_default.iloc[2] == times_noreconstruct.iloc[2]
    assert times_default.iloc[-3] == times_noreconstruct.iloc[-3]
    assert times_default.iloc[-2] == times_noreconstruct.iloc[-2]
    assert times_default.iloc[-1] == times_noreconstruct.iloc[-1]
    
    
    
    

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




