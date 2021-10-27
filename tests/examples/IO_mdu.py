# -*- coding: utf-8 -*-
"""
Created on Wed Oct 27 22:09:09 2021

@author: veenstra

tests whether mdu file can be imported and exported
"""

import os
from dfm_tools.io import mdu

dir_testinput = r'c:\DATA\dfm_tools_testdata'
dir_output = '.'

    
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

