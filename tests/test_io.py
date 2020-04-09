# -*- coding: utf-8 -*-
"""
Created on Wed Mar 18 16:33:52 2020

@author: veenstra
"""

import pytest
import inspect
import os

dir_testinput = os.path.join(r'c:/DATA','dfm_tools_testdata')
from dfm_tools.testutils import getmakeoutputdir




@pytest.mark.acceptance
def test_mdu():
    """
    tests whether mdu file can be imported and exported
    """

    dir_output = getmakeoutputdir(__file__,inspect.currentframe().f_code.co_name)
    #dir_output = './test_output'
    if not os.path.exists(dir_output):
        os.mkdir(dir_output)

    from dfm_tools.io import mdu
        
    try:
        filename_mdu = os.path.join(dir_testinput, r'DFM_3D_z_Grevelingen\computations\run01\Grevelingen-FM.mdu')
        data_mdu = mdu.read_deltares_ini(filename_mdu)
        print(data_mdu)
        import_success = True
    except:
        import_success = False

    try:
        filename_mdu_out = os.path.join(dir_output, 'Grevelingen-FM_out.mdu')
        mdu.write_deltares_ini(data_mdu, filename_mdu_out)
        export_success = True
    except:
        export_success = False
    
    assert import_success == True
    assert export_success == True






@pytest.mark.parametrize("file_pol", [pytest.param(os.path.join(dir_testinput,'DFM_3D_z_Grevelingen\\geometry\\structures\\Grevelingen-FM_BL_fxw.pliz'), id='Grevelingen pliz'),
                                      pytest.param(os.path.join(dir_testinput,'DFM_3D_z_Grevelingen\\geometry\\structures\\Grevelingen-FM_BL_fxw.pli'), id='Grevelingen pli'),
                                      pytest.param(os.path.join(dir_testinput,'world.ldb'), id='world'),
                                      pytest.param(os.path.join(dir_testinput,'Maeslant.tek'), id='Maeslant'),
                                      pytest.param(r'p:\1220688-lake-kivu\3_modelling\1_FLOW\7_heatfluxinhis\063_netcdf\shorterperiod\ballenplot\0200a.tek', id='Kivu tek map sal'),
                                      pytest.param(r'p:\1220688-lake-kivu\3_modelling\1_FLOW\7_heatfluxinhis\063_netcdf\shorterperiod\ballenplot\SDS-zd003b5dec2-sal.tek', id='Kivu tek map more'),
                                      pytest.param(os.path.join(dir_testinput,'test_new.tek'), id='ts_Theo')])
@pytest.mark.unittest
def test_readpolygon(file_pol):
    dir_output = getmakeoutputdir(__file__,inspect.currentframe().f_code.co_name)
    """
    this test tests if a netcdf varname can be retrieved from the 'dictionary' and if the variable can be retrieved from de netcdf
    file_pol = os.path.join(dir_testinput,'DFM_3D_z_Grevelingen\\geometry\\structures\\Grevelingen-FM_BL_fxw.pliz')
    file_pol = os.path.join(dir_testinput,'DFM_3D_z_Grevelingen\\geometry\\structures\\Grevelingen-FM_BL_fxw.pli')
    file_pol = os.path.join(dir_testinput,'world.ldb')
    file_pol = 'p:\\11205258-006-kpp2020_rmm-g6\\C_Work\\08_RMM_FMmodel\\geometry_j13_6-w3\\rmm_v1p3_fixed_weirs.pliz'
    file_pol = 'p:\\11205258-006-kpp2020_rmm-g6\\C_Work\\08_RMM_FMmodel\\geometry_j13_6-w3\\structures\\rmm_v1p3_structures.pli'
    file_pol = 'p:\\11205258-006-kpp2020_rmm-g6\\C_Work\\04_randvoorwaarden\\keringen\\Maeslantkering\\Maeslant.tek'
    file_pol = os.path.join(dir_testinput,'Maeslant.tek')
    file_pol = os.path.join(dir_testinput,'test_new.tek')
    file_pol = r'p:\1220688-lake-kivu\3_modelling\1_FLOW\7_heatfluxinhis\063_netcdf\shorterperiod\ballenplot\0200a.tek'
    file_pol = r'p:\1220688-lake-kivu\3_modelling\1_FLOW\7_heatfluxinhis\063_netcdf\shorterperiod\ballenplot\SDS-zd003b5dec2-sal.tek'

    dir_output = './test_output'
    """
    
    
    import matplotlib.pyplot as plt
    plt.close('all')
    
    from dfm_tools.io.polygon import Polygon
        
    pol_data_list, pol_name_list, pol_comment_list = Polygon.fromfile(file_pol, pd_output=False)
    pol_data_pd_list = Polygon.fromfile(file_pol, pd_output=True)
    pol_object_list = Polygon.fromfile(file_pol, obj_output=True)

    fig, ax = plt.subplots()
    for iP, pol_data in enumerate(pol_data_list):
        pd_collist = pol_data_pd_list[iP].columns.tolist()
        if 'datetime' in pd_collist:
            for iV in range(3,len(pd_collist)):
                ax.plot(pol_data_pd_list[iP].loc[:,'datetime'],pol_data_pd_list[iP].iloc[:,iV],'-',label=pd_collist[iV], linewidth=0.5)
        elif '*column 1 : X' in pd_collist or '*column 1 = x coordinate' in pd_collist:
            for iV in range(len(pd_collist)):
                ax.plot(pol_data_pd_list[iP].iloc[:,iV],'-',label=pd_collist[iV], linewidth=0.5)
        else: #eg ldb
            ax.plot(pol_data[:,0],pol_data[:,1],'-',linewidth=0.5)
    ax.legend()
    plt.savefig(os.path.join(dir_output,os.path.basename(file_pol).replace('.','')))




