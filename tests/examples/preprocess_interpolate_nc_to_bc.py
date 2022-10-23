# -*- coding: utf-8 -*-
"""
Created on Thu Aug 18 16:37:38 2022

@author: veenstra

This script can be used to interpolate pre-downloaded CMEMS data (or other netcdf files) to a boundary locations in a pli-file
"""

import datetime as dt
from pathlib import Path
import matplotlib.pyplot as plt
plt.close('all')
from dfm_tools.interpolate_grid2bnd import get_conversion_dict, interpolate_tide_to_bc, interpolate_nc_to_bc
from dfm_tools.hydrolib_helpers import forcinglike_to_DataFrame, forcinglike_to_DataArray, T3Dvector_to_T3Dtuple
from hydrolib.core.io.ext.models import Boundary, ExtModel

#TODO: add coordinate conversion of pli-coordinates (for nesting RD models)

model = 'CMCC' #CMEMS GFDL CMCC HYCOM

if model == 'HYCOM': #not available in dcsm area
    list_plifiles = [Path(r'c:\DATA\dfm_tools_testdata\GLBu0.08_expt_91.2\bcline.pli')]
else:
    #copied plifile from DCSM folder: r'p:\1204257-dcsmzuno\data\CMEMS\bnd\NorthSeaAndBaltic_1993-2019_20210510'
    #list_plifiles = [Path(r'n:\My Documents\werkmap\hydrolib_test\DCSM\DCSM-FM_OB_all_20181108.pli')] #TODO: reading this file results in empty Polyfile, should raise an error. https://github.com/Deltares/HYDROLIB-core/issues/320
    list_plifiles = [Path(r'n:\My Documents\werkmap\hydrolib_test\DCSM\DCSM-FM_OB_all_20181108_nocomments.pli')]

dir_out = r'n:\My Documents\werkmap\hydrolib_test\DCSM'
bc_type = 'bc' #currently only 'bc' supported #TODO: add netcdf bc support. https://github.com/Deltares/HYDROLIB-core/issues/318

refdate_str = 'minutes since 2011-12-22 00:00:00 +00:00' # if None, xarray uses ds.time.encoding['units'] as refdate_str

nPoints = 3 #amount of Points to process per PolyObject in the plifile (for testing, use None for all Points)

#quantities should be in conversion_dict.keys(). waterlevelbnd is steric/zos, tide is tidal components from FES/EOT
list_quantities = ['waterlevelbnd','salinitybnd','tide','ux,uy','temperaturebnd','tracerbndNO3']
#list_quantities = ['salinitybnd','tracerbndNO3']

dtstart = dt.datetime.now()
ext_bnd = ExtModel()


for file_pli in list_plifiles:
    for quantity in list_quantities:
        
        #TODO: take out of loop (but quantity is required)
        if model=='CMEMS': #2012-01-06 12:00:00 to 2013-01-03 12:00:00
            conversion_dict = get_conversion_dict()
            tstart = dt.datetime(2012, 1, 16, 12, 0)
            tstop = dt.datetime(2012, 4, 1, 12, 0)
            dir_sourcefiles_hydro = r'p:\1204257-dcsmzuno\data\CMEMS\nc\DCSM_allAvailableTimes' #CMEMS hydro: bottomT, so, thetao, uo, vo, zos (2012-01-06 12:00:00 to 2013-01-03 12:00:00) (daily values at noon, not at midnight)
            dir_pattern_hydro = Path(dir_sourcefiles_hydro,f'{conversion_dict[quantity]["ncvarname"]}_2012*.nc') # later remove 2012 from string, but this is faster for testing #TODO: it is quite slow, maybe speed up possible?
            dir_sourcefiles_waq = r'p:\11206304-futuremares\python_scripts\ocean_boundaryCMEMS\data_monthly' #CMEMS waq: no3, o2, phyc, so4, si (2011-12-16 12:00:00 to 2019-01-16 12:00:00)
            dir_pattern_waq = Path(dir_sourcefiles_waq,f'cmems_mod_glo_bgc_my_0.25_P1M-m_{conversion_dict[quantity]["ncvarname"]}_*.nc') 
        elif model=='GFDL':
            conversion_dict = get_conversion_dict()
            tstart = dt.datetime(2012, 1, 16, 12, 0)
            tstop = dt.datetime(2012, 4, 1, 12, 0)
            dir_sourcefiles_hydro = None
            dir_pattern_hydro = None
            dir_sourcefiles_waq = r'p:\11206304-futuremares\data\CMIP6_BC\GFDL-ESM4' #GFDL waq: no3 (1850-01-16 12:00:00 to 2014-12-16 12:00:00)
            dir_pattern_waq = Path(dir_sourcefiles_waq,f'{conversion_dict[quantity]["ncvarname"]}_esm-hist.nc')
        elif model=='CMCC': #TODO: check method, now finding nearest points (so always has values)
            conversion_dict = get_conversion_dict(ncvarname_updates={'salinitybnd':'sos', 'temperaturebnd':'tos'})
            tstart = dt.datetime(2015, 6, 16, 12, 0)
            tstop = dt.datetime(2016, 12, 1, 12, 0)
            dir_sourcefiles_hydro = r'p:\11206304-futuremares\data\CMIP6_BC\CMCC-ESM2'
            dir_pattern_hydro = Path(dir_sourcefiles_hydro,f'{conversion_dict[quantity]["ncvarname"]}_Omon_CMCC-ESM2_ssp126_r1i1p1f1_gn_*.nc')
            dir_sourcefiles_waq = dir_sourcefiles_hydro #CMCC waq: (2015-01-16 12:00:00 to 2100-12-16 12:00:00)
            dir_pattern_waq = dir_pattern_hydro
        elif model=='HYCOM':
            if quantity not in ['salinitybnd','temperaturebnd']: #only contains two quantities
                continue
            conversion_dict = get_conversion_dict(ncvarname_updates={'salinitybnd':'salinity', 'temperaturebnd':'water_temp'})
            tstart = dt.datetime(2016, 4, 20, 0, 0) #HYCOM
            tstop = dt.datetime(2016, 5, 3, 0, 0)
            dir_sourcefiles_hydro = 'c:\\DATA\\dfm_tools_testdata\\GLBu0.08_expt_91.2' #HYCOM hydro: salinity/so, water_temp/thetao (2016-04-19 00:00:00 to 2016-05-06 00:00:00)
            dir_pattern_hydro = Path(dir_sourcefiles_hydro,'HYCOM_ST_GoO_*.nc')
            dir_sourcefiles_waq = None
            dir_pattern_waq = None
        else:
            raise Exception(f'invalid model: {model}')
        
        print(f'processing quantity: {quantity}/{conversion_dict[quantity]["ncvarname"]}')
        if quantity in ['tide']: #tide #TODO: choose flexible/generic component notation
            tidemodel = 'FES2014' #FES2014, FES2012, EOT20
            component_list = ['2n2','mf','p1','m2','mks2','mu2','q1','t2','j1','m3','mm','n2','r2','k1','m4','mn4','s1','k2','m6','ms4','nu2','s2','l2','m8','msf','o1','s4'] #None results in all FES components
            ForcingModel_object = interpolate_tide_to_bc(tidemodel=tidemodel, file_pli=file_pli, component_list=component_list, nPoints=nPoints)
            for forcingobject in ForcingModel_object.forcing: #add A0 component
                forcingobject.datablock.append(['A0',0.0,0.0])
        elif quantity in ['waterlevelbnd','salinitybnd','temperaturebnd','ux,uy']: #hydro
            if dir_sourcefiles_hydro is None:
                continue
            ForcingModel_object = interpolate_nc_to_bc(dir_pattern=dir_pattern_hydro, file_pli=file_pli,
                                                       quantity=quantity, conversion_dict=conversion_dict,
                                                       tstart=tstart, tstop=tstop, refdate_str=refdate_str,
                                                       #reverse_depth=True, #to compare with coastserv files, this argument will be phased out
                                                       nPoints=nPoints)
        else: #waq
            if dir_pattern_waq is None:
                continue
            ForcingModel_object = interpolate_nc_to_bc(dir_pattern=dir_pattern_waq, file_pli=file_pli,
                                                       quantity=quantity, conversion_dict=conversion_dict,
                                                       tstart=tstart, tstop=tstop, refdate_str=refdate_str,
                                                       #reverse_depth=True, #to compare with coastserv files, this argument will be phased out
                                                       nPoints=nPoints)
        
        if 1: #plotting example data point
            for iF in [2]:#range(nPoints):
                forcingobject_one = ForcingModel_object.forcing[iF]
                if hasattr(forcingobject_one.quantityunitpair[1],'elementname'): #T3Dvector, take only first one
                    forcingobject_one, dummy = T3Dvector_to_T3Dtuple(forcingobject_one) #TODO: maybe make DataSet with two DataArrays instead
                forcingobject_one_xr = forcinglike_to_DataArray(forcingobject_one)
                fig,ax1 = plt.subplots()
                if hasattr(forcingobject_one,'vertpositions'):
                    pc = forcingobject_one_xr.T.plot.pcolormesh(ax=ax1)
                elif quantity=='tide': #TODO: improve tide/astronomic DataArray so this exception is not necessary anymore (e.g. DataSet with two DataArrays)
                    forcingobject_one_xr.isel(quantity=0).plot(ax=ax1, label='amplitude')
                    forcingobject_one_xr.isel(quantity=1).plot(ax=ax1, label='phase')
                    ax1.legend(loc=1)
                else:
                    forcingobject_one_xr.plot(ax=ax1)
        
        file_bc_basename = file_pli.name.replace('.pli','')
        if quantity=='tide':
            file_bc_out = Path(dir_out,f'{quantity}_{file_bc_basename}_{tidemodel}.bc')
        else:
            file_bc_out = Path(dir_out,f'{quantity}_{file_bc_basename}_{model}.bc')
        print(f'writing ForcingModel to bc file with hydrolib ({file_bc_out.name})')
        if bc_type=='bc':
            ForcingModel_object.save(filepath=file_bc_out)
            #TODO: improve formatting of bc file to make nicer, save diskspace and maybe write faster: https://github.com/Deltares/HYDROLIB-core/issues/308 (and https://github.com/Deltares/HYDROLIB-core/issues/313)
        else:
            raise Exception(f'invalid bc_type: {bc_type}')
        
        #make paths relative (sort of) (also necessary for locationfile) /../ should also be supported? 
        #ForcingModel_object.filepath = Path(str(ForcingModel_object.filepath).replace(dir_out,'')) #TODO: convert to relative paths in ext file possible? This path is the same as file_bc_out
        
        #generate boundary object for the ext file (quantity, pli-filename, bc-filename)
        boundary_object = Boundary(quantity=quantity,
                                   locationfile=Path(dir_out,file_pli.name),
                                   forcingfile=ForcingModel_object,
                                   )
        ext_bnd.boundary.append(boundary_object)

file_ext_out = Path(dir_out,'example_bnd.ext')
ext_bnd.save(filepath=file_ext_out)

time_passed = (dt.datetime.now()-dtstart).total_seconds()
print(f'>>time passed: {time_passed:.2f} sec')

