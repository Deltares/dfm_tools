# -*- coding: utf-8 -*-
"""
Created on Thu Aug 18 16:37:38 2022

@author: veenstra
"""

import datetime as dt
from pathlib import Path
import matplotlib.pyplot as plt
plt.close('all')
from dfm_tools.interpolate_grid2bnd import get_conversion_dict, interpolate_FES, interpolate_nc_to_bc
from hydrolib.core.io.ext.models import Boundary, ExtModel

dir_sourcefiles_hydro = r'p:\1204257-dcsmzuno\data\CMEMS\nc\DCSM_allAvailableTimes'
dir_sourcefiles_waq = r'p:\11206304-futuremares\python_scripts\ocean_boundaryCMEMS\data_monthly' #CMEMS
#dir_sourcefiles_waq = r'p:\11206304-futuremares\data\CMIP6_BC\GFDL-ESM4' #GFDL
#dir_sourcefiles_waq = r'p:\11206304-futuremares\data\CMIP6_BC\CMCC-ESM2' #CMCC

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

conversion_dict = get_conversion_dict()
list_modelvarnames = ['NO3']
list_modelvarnames = ['steric','salinity','tide']#,'tide']#,['salinity','temperature','steric'] #should be in varnames_dict.keys()
list_modelvarnames = ['tide']

ext_bnd = ExtModel()

for file_pli in list_plifiles:
    for modelvarname in list_modelvarnames:
        print(f'processing modelvarname: {modelvarname}')
        if modelvarname in ['tide']: #TODO: tide compares not too well, 2cm M2 difference. Why?
            dir_pattern,convert_180to360 = Path(r'p:\1230882-emodnet_hrsm\FES2014\fes2014_linux64_gnu\share\data\fes\2014\ocean_tide_extrapolated','*.nc'),True #TODO: or ocean_tide_extrapolated folder? (extrapolated to the coast)
            if 1:
                from hydrolib.core.io.polyfile.parser import read_polyfile
                import os
                import glob
                import xarray as xr
                file_list_nc = glob.glob(str(dir_pattern))
                component_list = [os.path.basename(x).replace('.nc','') for x in file_list_nc] #TODO: add sorting, manually? Add A0? translate dict for component names?
                #TODO: resulting amplitudes are slightly different, also imaginary numbers in orig code, why? c:\DATA\hydro_tools\FES\PreProcessing_FES_TideModel_imaginary.m
                #TODO: MV: "Cornelis gebruikt complexe getallen. Los van de interpolatie hoort dat hetzelfde te doen. Bij de interpolatie in de ruimte is het interpoleren van de complexe getallen hetzelfde als interpolatie van de coeffiecienten voor de cos en sin componenten. Die levert wel iets andere waarden dan wanneer je de amplitude en fase interpoleert."
                
                #load boundary file
                polyfile_object = read_polyfile(file_pli,has_z_values=False)
                
                pli_PolyObjects = polyfile_object['objects']
                for iPO, pli_PolyObject_sel in enumerate(pli_PolyObjects):
                    print(f'processing PolyObject {iPO+1} of {len(pli_PolyObjects)}: name={pli_PolyObject_sel.metadata.name}')
                    
                    for iP, pli_Point_sel in enumerate(pli_PolyObject_sel.points[:nPoints]):
                        print(f'processing Point {iP+1} of {len(pli_PolyObject_sel.points)}: ',end='')
                        
                        lonx, laty = pli_Point_sel.x, pli_Point_sel.y
                        if convert_180to360:
                            lonx = lonx%360 #for FES since it ranges from 0 to 360 instead of -180 to 180
                        print(f'(x={lonx}, y={laty})')
                        pli_PolyObject_name_num = f'{pli_PolyObject_sel.metadata.name}_{iP+1:04d}'
                        
                        datablock_list = []
                        for iC,component in enumerate(component_list):
                            file_nc = os.path.join(os.path.dirname(dir_pattern),f'{component}.nc')
                            data_xr = xr.open_dataset(file_nc)
                            """
                            def extract_component(ds):
                                comp_name = os.path.basename(x).replace('.nc','') #x here is the filename to extract the componentname from, but that is not available in the dataset
                                ds.assign(component=comp_name)
                                return ds
                            data_xrall = xr.open_mfdataset(file_list_nc, concat_dim='component',preprocess=extract_component)
                            """
                            lonvar_vals = data_xr['lon'].to_numpy()
                            latvar_vals = data_xr['lat'].to_numpy()
                            data_xr_amp = data_xr['amplitude']
                            data_xr_phs = data_xr['phase']
                            breakit
                
            ForcingModel_object = interpolate_FES(dir_pattern, file_pli, convert_180to360=convert_180to360, nPoints=nPoints, debug=debug)
        elif modelvarname in ['NO3']:
            varname_file = conversion_dict[modelvarname]['substance'][0] #TODO: [1] is also necessary for uxuy
            dir_pattern,convert_180to360 = Path(dir_sourcefiles_waq,f'cmems_mod_glo_bgc_my_0.25_P1M-m_{varname_file}_*.nc'),False # CMEMS waq
            #dir_pattern,convert_180to360 = Path(dir_sourcefiles_waq,f'{varname_file}_esm-hist.nc'),True # GFDL
            #dir_pattern,convert_180to360 = Path(dir_sourcefiles_waq,f'{varname_file}_Omon_CMCC-ESM2_ssp126_r1i1p1f1_gn_*.nc'),True #CMCC, TODO: crashes because of missing lat coords
            
            if 0:
                import glob
                import xarray as xr
                file_list_nc = glob.glob(str(dir_pattern))
                data_xr = xr.open_mfdataset(file_list_nc)#,decode_times=False)#, combine='by_coords', decode_times=False)
                breakit
            ForcingModel_object = interpolate_nc_to_bc(dir_pattern=dir_pattern, file_pli=file_pli, modelvarname=modelvarname,
                                                       convert_180to360=convert_180to360,
                                                       tstart=tstart, tstop=tstop, refdate_str=refdate_str,
                                                       nPoints=nPoints, debug=debug)
        else: #['steric','salinity','temperature'] ['uxuy']
            varname_file = conversion_dict[modelvarname]['substance'][0] #TODO: [1] is also necessary for uxuy
            dir_pattern,convert_180to360 = Path(dir_sourcefiles_hydro,f'{varname_file}_1993*.nc'),False # later remove 1993 from string, but this is faster for testing
            ForcingModel_object = interpolate_nc_to_bc(dir_pattern=dir_pattern, file_pli=file_pli, modelvarname=modelvarname, 
                                                       convert_180to360=convert_180to360,
                                                       tstart=tstart, tstop=tstop, refdate_str=refdate_str,
                                                       nPoints=nPoints, debug=debug)
            #ForcingModel_object.filepath = Path(str(ForcingModel_object.filepath).replace(dir_out,'')) #TODO: convert to relative paths in ext file possible?
        
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
        
        #generate boundary object for the ext file (quantity, pli-filename, bc-filename)
        boundary_object = Boundary(quantity=modelvarname, #TODO: nodeId / bndWidth1D / bndBlDepth are written as empty values, but they should not be written if not supplied. https://github.com/Deltares/HYDROLIB-core/issues/319
                                   locationfile=Path(dir_out,file_pli.name),
                                   forcingfile=ForcingModel_object,
                                   )
        ext_bnd.boundary.append(boundary_object)

file_ext_out = Path(dir_out,'example_bnd.ext')
ext_bnd.save(filepath=file_ext_out)




