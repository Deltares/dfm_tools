# -*- coding: utf-8 -*-
"""
Created on Thu Aug 18 16:37:38 2022

@author: veenstra
"""

import datetime as dt
#import os
from pathlib import Path
import matplotlib.pyplot as plt
plt.close('all')
from dfm_tools.CMEMS_interpolate import get_conversions_dicts, interpolate_FES, interpolate_nc_to_bc
from hydrolib.core.io.ext.models import Boundary, ExtModel
#TODO: add other issues in other scripts
#TODO: add uxuy support (merged arrays) >> or avoid by using separate blocks if accepted by dflowfm. https://github.com/Deltares/HYDROLIB-core/issues/316

convert_180to360 = False

dir_sourcefiles_hydro = r'p:\1204257-dcsmzuno\data\CMEMS\nc\DCSM_allAvailableTimes'
dir_sourcefiles_waq = r'p:\11206304-futuremares\python_scripts\ocean_boundaryCMEMS\data_monthly' #CMEMS
#dir_sourcefiles_waq, convert_180to360 = r'p:\11206304-futuremares\data\CMIP6_BC\GFDL-ESM4', True #GFDL
#dir_sourcefiles_waq, convert_180to360 = r'p:\11206304-futuremares\data\CMIP6_BC\CMCC-ESM2', True #CMCC

#copied plifile from DCSM folder: r'p:\1204257-dcsmzuno\data\CMEMS\bnd\NorthSeaAndBaltic_1993-2019_20210510'
#list_plifiles = [Path(r'n:\My Documents\werkmap\hydrolib_test\DCSM\DCSM-FM_OB_all_20181108.pli')] #TODO: reading this file results in empty Polyfile, should raise an error. https://github.com/Deltares/HYDROLIB-core/issues/320
list_plifiles = [Path(r'n:\My Documents\werkmap\hydrolib_test\DCSM\DCSM-FM_OB_all_20181108_nocomments.pli')]

dir_out = r'n:\My Documents\werkmap\hydrolib_test\DCSM'
bc_type = 'bc' #currently only 'bc' supported #TODO: add netcdf bc support. https://github.com/Deltares/HYDROLIB-core/issues/318

refdate_str = 'minutes since 2011-12-22 00:00:00 +00:00' # this is copied from the reference bc file, but can be changed by the user
tstart = dt.datetime(1993, 1, 1, 12, 0) #CMEMS phys has daily values at 12:00 (not at midnight), so make sure to include a day extra if necessary
tstop = dt.datetime(1993, 2, 1, 12, 0)
#tstart = dt.datetime(2011, 12, 16, 12, 0)
#tstop = dt.datetime(2012, 12, 1, 12, 0)
#tstart = dt.datetime(2015, 6, 16, 12, 0)
#tstop = dt.datetime(2015, 12, 1, 12, 0)
nPoints = 2 #amount of Points to process per PolyObject in the plifile (for testing, use None for all Points)
debug = False

usefor, constituent_boundary_type = get_conversions_dicts()
list_modelvarnames = ['NO3']
list_modelvarnames = ['steric','salinity','tide']#,'tide']#,['salinity','temperature','steric'] #should be in varnames_dict.keys()

ext_bnd = ExtModel()

for file_pli in list_plifiles:
    for modelvarname in list_modelvarnames:
        print(f'processing modelvarname: {modelvarname}')
        if modelvarname in ['tide']: #TODO: tide compares not too well, 2cm M2 difference. Why?
            dir_pattern = Path(r'p:\1230882-emodnet_hrsm\FES2014\fes2014_linux64_gnu\share\data\fes\2014\ocean_tide_extrapolated','*.nc') #TODO: or ocean_tide_extrapolated folder?
            ForcingModel_object = interpolate_FES(dir_pattern, file_pli, convert_180to360=True, nPoints=nPoints, debug=debug)
        elif modelvarname in ['NO3']:
            varname_file = usefor[modelvarname]['substance'][0] #TODO: [1] is also necessary for uxuy
            dir_pattern = Path(dir_sourcefiles_waq,f'cmems_mod_glo_bgc_my_0.25_P1M-m_{varname_file}_*.nc') # CMEMS waq
            #dir_pattern = Path(dir_sourcefiles_waq,f'{varname_file}_esm-hist.nc') # GFDL. crashes because of NoLeap calendar
            #dir_pattern = Path(dir_sourcefiles_waq,f'{varname_file}_Omon_CMCC-ESM2_ssp126_r1i1p1f1_gn_*.nc') #CMCC, crashes because of missing lat coords
            ForcingModel_object = interpolate_nc_to_bc(dir_pattern=dir_pattern, file_pli=file_pli, modelvarname=modelvarname,
                                                       convert_180to360=convert_180to360,
                                                       tstart=tstart, tstop=tstop, refdate_str=refdate_str,
                                                       nPoints=nPoints, debug=debug)
        else: #['steric','salinity','temperature'] ['uxuy']
            varname_file = usefor[modelvarname]['substance'][0] #TODO: [1] is also necessary for uxuy
            dir_pattern = Path(dir_sourcefiles_hydro,f'{varname_file}_1993*.nc') # later remove 1993 from string, but this is faster for testing
            ForcingModel_object = interpolate_nc_to_bc(dir_pattern=dir_pattern, file_pli=file_pli, modelvarname=modelvarname, 
                                                       convert_180to360=convert_180to360,
                                                       tstart=tstart, tstop=tstop, refdate_str=refdate_str,
                                                       nPoints=nPoints, debug=debug)
            #ForcingModel_object.filepath = Path(str(ForcingModel_object.filepath).replace(dir_out,'')) #TODO QUESTION: add relative paths in ext file (already possible?)
        
        file_bc_basename = file_pli.name.replace('.pli','.bc')
        file_bc_out = Path(dir_out,f'{modelvarname}_{file_bc_basename}')
        print(f'writing ForcingModel to bc file with hydrolib ({file_bc_out.name})')
        dtstart = dt.datetime.now()
        if bc_type=='bc':
            ForcingModel_object.save(filepath=file_bc_out)
            #TODO: improve performance of bc file writing with hydrolib (numpy.savetxt() is way faster because of formatting): https://github.com/Deltares/HYDROLIB-core/issues/313
            #TODO: improve formatting of bc file to make nicer, save diskspace and maybe write faster: https://github.com/Deltares/HYDROLIB-core/issues/308
        else:
            raise Exception(f'invalid bc_type: {bc_type}')
        time_passed = (dt.datetime.now()-dtstart).total_seconds()
        #print(f'>>time passed: {time_passed:.2f} sec')
        
        boundary_object = Boundary(quantity=modelvarname, #TODO: nodeId / bndWidth1D / bndBlDepth are written as empty values, but they should not be written if not supplied. https://github.com/Deltares/HYDROLIB-core/issues/319
                                   locationfile=Path(dir_out,file_pli.name),
                                   forcingfile=ForcingModel_object,
                                   )
        ext_bnd.boundary.append(boundary_object)

file_ext_out = Path(dir_out,'example_bnd.ext')
ext_bnd.save(filepath=file_ext_out)




