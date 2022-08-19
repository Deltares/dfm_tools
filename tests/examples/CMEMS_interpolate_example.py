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
from dfm_tools.CMEMS_interpolate import get_varnames_dict, interpolate_FES, interpolate_nc_to_bc
from hydrolib.core.io.ext.models import Boundary, ExtModel
#TODO: improve formatting of bcfile (and other issues in other scripts)
#TODO REPORT: add uxuy functionality. How to write bc file with merged arrays?

dir_sourcefiles = r'p:\1204257-dcsmzuno\data\CMEMS\nc\DCSM_allAvailableTimes'

#copied plifile from DCSM folder: r'p:\1204257-dcsmzuno\data\CMEMS\bnd\NorthSeaAndBaltic_1993-2019_20210510'
#list_plifiles = [Path(r'n:\My Documents\werkmap\hydrolib_test\DCSM\DCSM-FM_OB_all_20181108.pli')] #TODO: reading this file results in empty Polyfile, should raise an error (report issue)
list_plifiles = [Path(r'n:\My Documents\werkmap\hydrolib_test\DCSM\DCSM-FM_OB_all_20181108_nocomments.pli')]

dir_out = r'n:\My Documents\werkmap\hydrolib_test\DCSM'
bc_type = 'bc' #currently only 'bc' supported #TODO: add netcdf bc support (is this already in hydrolib?)

refdate_str = 'minutes since 2011-12-22 00:00:00 +00:00' # this is copied from the reference bc file, but can be changed by the user
tstart = dt.datetime(1993, 1, 1, 12, 0) #CMEMS has daily values at 12:00 (not at midnight), so make sure to include a day extra if necessary
#tstop = dt.datetime(2019, 12, 31, 12, 0)
tstop = dt.datetime(1993, 2, 1, 12, 0)
nPoints = 2 #amount of Points to process per PolyObject in the plifile (for testing, use None for all Points)

varnames_dict = get_varnames_dict('cmems')
list_modelvarnames = ['tide','salinity','steric']#,['salinity','temperature','steric'] #should be in varnames_dict.keys()

ext_bnd = ExtModel()

for file_pli in list_plifiles:
    for modelvarname in list_modelvarnames:
        print(f'processing modelvarname: {modelvarname}')
        if modelvarname == 'tide': #TODO: tide compares not too well, 2cm M2 difference. Why?
            dir_pattern = Path(r'p:\1230882-emodnet_hrsm\FES2014\fes2014_linux64_gnu\share\data\fes\2014\ocean_tide','*.nc') #TODO: or ocean_tide_extrapolated folder?
            ForcingModel_object = interpolate_FES(dir_pattern, file_pli, nPoints=nPoints, debug=True)
            
        else:
            file_pattern = f'{varnames_dict[modelvarname]}_1993*.nc' # later remove 1993 from string, but this is faster for testing
            dir_pattern = Path(dir_sourcefiles,file_pattern)
                    
            ForcingModel_object = interpolate_nc_to_bc(dir_pattern=dir_pattern, file_pli=file_pli, 
                                                       modelvarname=modelvarname, varnames_dict=varnames_dict,
                                                       tstart=tstart, tstop=tstop, refdate_str=refdate_str,
                                                       nPoints=nPoints, debug=True)#, ForcingModel_object=None)
        
        file_bc_basename = file_pli.name.replace('.pli','.bc')
        file_bc_out = Path(dir_out,f'{modelvarname}_{file_bc_basename}')
        print(f'writing ForcingModel to bc file with hydrolib ({file_bc_out.name})')
        dtstart = dt.datetime.now()
        if bc_type=='bc':
            ForcingModel_object.save(filepath=file_bc_out)
        else:
            raise Exception(f'invalid bc_type: {bc_type}')
        time_passed = (dt.datetime.now()-dtstart).total_seconds()
        #print(f'>>time passed: {time_passed:.2f} sec')
        
        
        """
        quantity: str = Field(alias="quantity")
        nodeid: Optional[str] = Field(alias="nodeId")
        locationfile: Optional[Path] = Field(alias="locationFile")
        forcingfile: ForcingModel = Field(alias="forcingFile")
        bndwidth1d: Optional[float] = Field(alias="bndWidth1D")
        bndbldepth: Optional[float] = Field(alias="bndBlDepth")
        """
        boundary_object = Boundary(quantity=modelvarname, #TODO REPORT: nodeId / bndWidth1D / bndBlDepth are written as empty values, but they should not be written if not supplied
                                   locationfile=Path(dir_out,file_pli.name),
                                   forcingfile=ForcingModel_object,
                                   )
        ext_bnd.boundary.append(boundary_object)

file_ext_out = Path(dir_out,'example_bnd.ext')
ext_bnd.save(filepath=file_ext_out)




