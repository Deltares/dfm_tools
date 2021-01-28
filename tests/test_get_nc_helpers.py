# -*- coding: utf-8 -*-
"""
Created on Fri Mar 20 08:57:42 2020

@author: veenstra
"""

import pytest
#import inspect
import os

if 'TEAMCITY_VERSION' in os.environ.keys(): #teamcity path
    dir_testinput = r'/opt/testdata/dfm_tools'
else: #default to this path
    dir_testinput = os.path.join(r'c:/DATA','dfm_tools_testdata')

#from dfm_tools.testutils import getmakeoutputdir



@pytest.mark.unittest    
def test_gettimesfromnc():

    from dfm_tools.get_nc import get_ncmodeldata
    
    #this file seems to have a reconstructable array, but contains gaps, so the time reader falls back to reading all
    file_nc = r'p:\1204257-dcsmzuno\data\waterbase2\sea_surface_height\id1-LICHTELGRE.nc'
    data_fromhis_time = get_ncmodeldata(file_nc,varname='time',timestep='all')

    
    
    
    

@pytest.mark.unittest    
def test_getvarnamemapnc():
    """
    this test tests if a netcdf varname can be retrieved from the 'dictionary' and if the variable can be retrieved from de netcdf
    """
    
    from netCDF4 import Dataset
    
    from dfm_tools.get_nc_helpers import get_varname_fromnc
    
    file_map = os.path.join(dir_testinput,r'DFM_3D_z_Grevelingen\computations\run01\DFM_OUTPUT_Grevelingen-FM\Grevelingen-FM_0000_map.nc')
    #file_net = os.path.join(dir_testinput,r'DFM_3D_z_Grevelingen\computations\run01\Grevelingen_FM_grid_20190603_net.nc')
    
    file_nc = file_map
    data_nc = Dataset(file_nc)
    varname_requested = 'NetNode_y' #is actually in file, so this is not a good test
    
    varname = get_varname_fromnc(data_nc,varname_requested, vardim='var')
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
    vars_pd_matching = vars_pd[vars_pd.loc[:,'long_name'].str.match(pattern)] #does not have to stop after pattern
    #vars_pd_matching = vars_pd[vars_pd.loc[:,'long_name'].str.startswith('Flow') & vars_pd.loc[:,'long_name'].str.endswith('component')]
    varkeys_list_matching = list(vars_pd_matching['nc_varkeys'])
    
    assert varkeys_list_matching == ['mesh2d_ucx', 'mesh2d_ucy', 'mesh2d_ucz', 'mesh2d_ucxa', 'mesh2d_ucya']




