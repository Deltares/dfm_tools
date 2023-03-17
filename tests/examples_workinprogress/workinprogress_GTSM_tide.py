# -*- coding: utf-8 -*-
"""
Created on Fri Mar 17 16:03:56 2023

@author: veenstra
"""

import os
import datetime as dt
import contextily as ctx
from pathlib import Path
import matplotlib.pyplot as plt
plt.close('all')
import dfm_tools as dfmt
import hydrolib.core.dflowfm as hcdfm

#TODO: add coordinate conversion of pli-coordinates? (for nesting RD models in oceanmodels)
#TODO: additional models/sources for download/interpolate (evt xESMF for CMCC, climate forcing cmip6 procedure (=calendarconversion) and others)

nPoints = 3# None #amount of Points to process per PolyObject in the plifile (use int for testing, use None for all Points)
refdate_str = 'minutes since 2011-12-22 00:00:00 +00:00' # if None, xarray uses ds.time.encoding['units'] as refdate_str
dir_output = './test_interpolate_nc_to_bc_GTSMtide'

list_plifiles = [Path(r'p:\11208054-004-dcsm-fm\models\model_input\bnd_cond\pli\DCSM-FM_OB_all_20181108_nocomments.pli')] #TODO: reading this file without '_nocomments' results in empty Polyfile, should raise an error. https://github.com/Deltares/HYDROLIB-core/issues/320

# start of interpolation process
dtstart = dt.datetime.now()
ext_bnd = hcdfm.ExtModel()
if not os.path.isdir(dir_output):
    os.mkdir(dir_output)

for file_pli in list_plifiles:
    file_bc_basename = file_pli.name.replace('.pli','')
    
    quantity = 'tide'
    print(f'processing quantity: {quantity}')
    tidemodel = 'GTSM4.1preliminary' #FES2014, FES2012, EOT20
    component_list = None
    ForcingModel_object = dfmt.interpolate_tide_to_bc(tidemodel=tidemodel, file_pli=file_pli, component_list=component_list, nPoints=nPoints)
    for forcingobject in ForcingModel_object.forcing: #add A0 component
        forcingobject.datablock.append(['A0',0.0,0.0])
    
    file_bc_basename = file_pli.name.replace('.pli','')
    file_bc_out = Path(dir_output,f'{quantity}_{file_bc_basename}_{tidemodel}.bc')
    
    print(f'writing ForcingModel to bc file with hydrolib ({file_bc_out.name})')
    bc_type = 'bc' #TODO: add netcdf bc support. https://github.com/Deltares/HYDROLIB-core/issues/318
    if bc_type=='bc':
        #ForcingModel_object.serializer_config.float_format = '.3f' #TODO SOLVED: improve formatting of bc file: https://github.com/Deltares/HYDROLIB-core/issues/308
        #ForcingModel_object.serializer_config.float_format_datablock = '.5f' #maybe move this to interp_regularnc_to_plipoints/interpolate_tide_to_bc?
        ForcingModel_object.save(filepath=file_bc_out)
    
    #TODO: support for relative paths?
    #generate boundary object for the ext file (quantity, pli-filename, bc-filename)
    boundary_object = hcdfm.Boundary(quantity=quantity.replace('tide','waterlevelbnd'), #the FM quantity for tide is also waterlevelbnd
                                     locationfile=file_pli,
                                     forcingfile=ForcingModel_object)
    ext_bnd.boundary.append(boundary_object)

file_ext_out = Path(dir_output,'example_bnd.ext')
ext_bnd.save(filepath=file_ext_out)

time_passed = (dt.datetime.now()-dtstart).total_seconds()
print(f'>>total script time passed: {time_passed:.2f} sec')

