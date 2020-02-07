#!/usr/bin/env python

"""Tests for `dfm_tools` package."""

import pytest


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
    
    
    # 4. Vefiry final expectations
    assert 1==1


def Test_UGrid(self):
    #from netCDF4 import Dataset
    from dfm_tools.grid import UGrid
    
    #file_map = r'n:\My Documents\werkmap\vanJulien_shortmodelfiles\DFM_3D_z_Grevelingen\computations\run01\DFM_OUTPUT_Grevelingen-FM\Grevelingen-FM_0000_map.nc'
    file_net = r'n:\My Documents\werkmap\vanJulien_shortmodelfiles\DFM_3D_z_Grevelingen\computations\run01\Grevelingen_FM_grid_20190603_net.nc'
    #file_net = r'p:\11205258-006-kpp2020_rmm-g6\C_Work\01_Rooster\final_totaalmodel\rooster_rmm_v1p5_net.nc'
    
    data_ncUG = UGrid.fromfile(file_net)
    test = data_ncUG.cellcoords()
    assert 1==1


def Test_mdu(self):
    #from netCDF4 import Dataset
    from dfm_tools.mdu import read_deltares_ini, write_deltares_ini
    
    filename_mdu = r'n:\My Documents\werkmap\vanJulien_shortmodelfiles\DFM_3D_z_Grevelingen\computations\run01\Grevelingen-FM.mdu'
    filename_mdu_out = r'n:\My Documents\werkmap\vanJulien_shortmodelfiles\DFM_3D_z_Grevelingen\computations\run01\Grevelingen-FM_out.mdu'
    
    data_mdu = read_deltares_ini(filename_mdu)
    write_deltares_ini(data_mdu, filename_mdu_out)
    
    
    assert 1==1


def Test_dflowutil_mesh(self):
    from dfm_tools.grid import plot_net
    
    #plt.close('all')
    file_net = r'n:\My Documents\werkmap\vanJulien_shortmodelfiles\DFM_3D_z_Grevelingen\computations\run01\Grevelingen_FM_grid_20190603_net.nc'
    file_map = r'n:\My Documents\werkmap\vanJulien_shortmodelfiles\DFM_3D_z_Grevelingen\computations\run01\DFM_OUTPUT_Grevelingen-FM\Grevelingen-FM_0000_map.nc'
    
    file_nc = file_net
    plot_net(file_nc)
    file_nc = file_map
    plot_net(file_nc,var_values='mesh2d_s1')
    
    
    
    
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
