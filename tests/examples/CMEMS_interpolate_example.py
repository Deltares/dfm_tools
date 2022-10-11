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
from dfm_tools.hydrolib_helpers import forcingobject_to_dataframe
from hydrolib.core.io.ext.models import Boundary, ExtModel

model = 'CMEMS' #CMEMS GFDL CMCC HYCOM

#copied plifile from DCSM folder: r'p:\1204257-dcsmzuno\data\CMEMS\bnd\NorthSeaAndBaltic_1993-2019_20210510'
#list_plifiles = [Path(r'n:\My Documents\werkmap\hydrolib_test\DCSM\DCSM-FM_OB_all_20181108.pli')] #TODO: reading this file results in empty Polyfile, should raise an error. https://github.com/Deltares/HYDROLIB-core/issues/320
list_plifiles = [Path(r'n:\My Documents\werkmap\hydrolib_test\DCSM\DCSM-FM_OB_all_20181108_nocomments.pli')]

dir_out = r'n:\My Documents\werkmap\hydrolib_test\DCSM'
bc_type = 'bc' #currently only 'bc' supported #TODO: add netcdf bc support. https://github.com/Deltares/HYDROLIB-core/issues/318

refdate_str = 'minutes since 2011-12-22 00:00:00 +00:00' # this is copied from the reference bc file, but can be changed by the user
tstart = dt.datetime(1993, 1, 1, 12, 0) #CMEMS phys has daily values at 12:00 (not at midnight), so make sure to include a day extra if necessary. also NO3_GFDL
tstop = dt.datetime(1993, 3, 1, 12, 0)
#tstart = dt.datetime(2011, 12, 16, 12, 0) #NO3_CMEMS
#tstop = dt.datetime(2012, 12, 1, 12, 0)
#tstart = dt.datetime(2015, 6, 16, 12, 0)
#tstop = dt.datetime(2015, 12, 1, 12, 0)
#tstart = dt.datetime(2016, 4, 20, 0, 0) #HYCOM
#tstop = dt.datetime(2016, 5, 3, 0, 0)

nPoints = 3 #amount of Points to process per PolyObject in the plifile (for testing, use None for all Points)
debug = False

#list_quantities = ['NO3']
list_quantities = ['steric','salinity','tide']#,['salinity','temperature','steric'] #should be in conversion_dict.keys()
#list_quantities = ['salinity','temperature']

dtstart = dt.datetime.now()
ext_bnd = ExtModel()

for file_pli in list_plifiles:
    for quantity in list_quantities:
        conversion_dict = get_conversion_dict()
        ncvarname = conversion_dict[quantity]['ncvarname']
        bcvarname = conversion_dict[quantity]['bcvarname']
        print(f'processing quantity: {quantity}/{ncvarname}/{bcvarname}')
        if quantity in ['tide']:
            dir_pattern,convert_360to180 = Path(r'P:\metocean-data\licensed\FES2014','*.nc'),True #source: p:\1230882-emodnet_hrsm\FES2014\fes2014_linux64_gnu\share\data\fes\2014\ocean_tide_extrapolated
            #dir_pattern,convert_360to180 = Path(r'P:\metocean-data\open\FES2012\data','*_FES2012_SLEV.nc'),True #is eigenlijk ook licensed
            #dir_pattern,convert_360to180 = Path(r'P:\metocean-data\open\EOT20\ocean_tides','*_ocean_eot20.nc'),True
            component_list = ['2n2','mf','p1','m2','mks2','mu2','q1','t2','j1','m3','mm','n2','r2','k1','m4','mn4','s1','k2','m6','ms4','nu2','s2','l2','m8','msf','o1','s4'] #None results in all FES components
            ForcingModel_object = interpolate_FES(dir_pattern, file_pli, component_list=component_list, convert_360to180=convert_360to180, nPoints=nPoints, debug=debug)
            for forcingobject in ForcingModel_object.forcing: #add A0 component
                forcingobject.datablock.append(['A0',0.0,0.0])
        elif quantity in ['NO3']:
            if model=='CMEMS':
                dir_sourcefiles_waq = r'p:\11206304-futuremares\python_scripts\ocean_boundaryCMEMS\data_monthly' #CMEMS
                dir_pattern,convert_360to180 = Path(dir_sourcefiles_waq,f'cmems_mod_glo_bgc_my_0.25_P1M-m_{ncvarname}_*.nc'),False # CMEMS waq
                #dir_pattern,convert_360to180 = Path(dir_sourcefiles_waq,f'{ncvarname}_*.nc'),False # CMEMS waq
            elif model=='GFDL':
                dir_sourcefiles_waq = r'p:\11206304-futuremares\data\CMIP6_BC\GFDL-ESM4' #GFDL
                dir_pattern,convert_360to180 = Path(dir_sourcefiles_waq,f'{ncvarname}_esm-hist.nc'),True # GFDL
            elif model=='CMCC':
                dir_sourcefiles_waq = r'p:\11206304-futuremares\data\CMIP6_BC\CMCC-ESM2' #CMCC
                dir_pattern,convert_360to180 = Path(dir_sourcefiles_waq,f'{ncvarname}_Omon_CMCC-ESM2_ssp126_r1i1p1f1_gn_*.nc'),True #CMCC, TODO: crashes because of missing lat coords
            ForcingModel_object = interpolate_nc_to_bc(dir_pattern=dir_pattern, file_pli=file_pli, quantity=quantity,
                                                       convert_360to180=convert_360to180,
                                                       tstart=tstart, tstop=tstop, refdate_str=refdate_str,
                                                       reverse_depth=True, #to compare with coastserv files, this argument will be phased out
                                                       nPoints=nPoints, debug=debug)
        else: #['steric','salinity','temperature'] ['uxuy']
            if model=='CMEMS':
                dir_sourcefiles_hydro = r'p:\1204257-dcsmzuno\data\CMEMS\nc\DCSM_allAvailableTimes' #CMEMS
                dir_pattern,convert_360to180 = Path(dir_sourcefiles_hydro,f'{ncvarname}_1993*.nc'),False # CMEMS. later remove 1993 from string, but this is faster for testing
            elif model=='HYCOM':
                dir_sourcefiles_hydro = 'c:\\DATA\\dfm_tools_testdata\\GLBu0.08_expt_91.2' #HYCOM
                dir_pattern,convert_360to180 = Path(dir_sourcefiles_hydro,'HYCOM_ST_GoO_*.nc'),False # HYCOM
                file_pli = Path(r'c:\DATA\dfm_tools_testdata\GLBu0.08_expt_91.2\bcline.pli')
                conversion_dict = get_conversion_dict('HYCOM') #rename_variables = {'salinity':'so', 'water_temp':'thetao'}
            ForcingModel_object = interpolate_nc_to_bc(dir_pattern=dir_pattern, file_pli=file_pli, quantity=quantity, 
                                                       convert_360to180=convert_360to180, conversion_dict=conversion_dict,
                                                       tstart=tstart, tstop=tstop, refdate_str=refdate_str,
                                                       reverse_depth=True, #to compare with coastserv files, this argument will be phased out
                                                       nPoints=nPoints, debug=debug)
        
        if 1: #plotting example data point
            for iF in [2]:#range(nPoints):
                forcingobject_one = ForcingModel_object.forcing[iF]
                forcingobject_one_df = forcingobject_to_dataframe(forcingobject_one)
                fig,ax1 = plt.subplots()
                if hasattr(forcingobject_one,'verticalpositionspecification'):
                    pc = ax1.pcolormesh(forcingobject_one_df.index,forcingobject_one.verticalpositionspecification,forcingobject_one_df.T)
                    fig.colorbar(pc,ax=ax1)
                else:
                    forcingobject_one_df.plot(ax=ax1)
                ax1.set_title(forcingobject_one.quantityunitpair[1].quantity)
        
        file_bc_basename = file_pli.name.replace('.pli','.bc')
        file_bc_out = Path(dir_out,f'{quantity}_{file_bc_basename}')
        print(f'writing ForcingModel to bc file with hydrolib ({file_bc_out.name})')
        if bc_type=='bc':
            ForcingModel_object.save(filepath=file_bc_out)
            #TODO: improve formatting of bc file to make nicer, save diskspace and maybe write faster: https://github.com/Deltares/HYDROLIB-core/issues/308 (and https://github.com/Deltares/HYDROLIB-core/issues/313)
        else:
            raise Exception(f'invalid bc_type: {bc_type}')
        
        #make paths relative (sort of) (also necessary for locationfile) /../ should also be supported? 
        #ForcingModel_object.filepath = Path(str(ForcingModel_object.filepath).replace(dir_out,'')) #TODO: convert to relative paths in ext file possible? This path is the same as file_bc_out
        
        #generate boundary object for the ext file (quantity, pli-filename, bc-filename)
        boundary_object = Boundary(quantity=bcvarname, #TODO: nodeId / bndWidth1D / bndBlDepth are written as empty values, but they should not be written if not supplied. https://github.com/Deltares/HYDROLIB-core/issues/319
                                   locationfile=Path(dir_out,file_pli.name),
                                   forcingfile=ForcingModel_object,
                                   )
        ext_bnd.boundary.append(boundary_object)

file_ext_out = Path(dir_out,'example_bnd.ext')
ext_bnd.save(filepath=file_ext_out)

time_passed = (dt.datetime.now()-dtstart).total_seconds()
print(f'>>time passed: {time_passed:.2f} sec')

