# -*- coding: utf-8 -*-
"""
Created on Thu Aug 18 16:37:38 2022

@author: veenstra

This script can be used to interpolate pre-downloaded CMEMS data (or other netcdf files) to points of a pli-file, resulting in a boundary forcingfile (bc file)
"""

import os
import datetime as dt
import contextily as ctx
from pathlib import Path
import matplotlib.pyplot as plt
plt.close('all')
import dfm_tools as dfmt
try: #0.3.1 release
    from hydrolib.core.io.ext.models import Boundary, ExtModel
except: #main branch and next release #TODO: move to easy imports after https://github.com/Deltares/HYDROLIB-core/issues/410
    from hydrolib.core.dflowfm.ext.models import Boundary, ExtModel

#TODO: add coordinate conversion of pli-coordinates? (for nesting RD models in oceanmodels)
#TODO: additional models/sources for download/interpolate (evt xESMF for CMCC, climate forcing cmip6 procedure (=calendarconversion) and others)

nPoints = 3# None #amount of Points to process per PolyObject in the plifile (use int for testing, use None for all Points)
refdate_str = 'minutes since 2011-12-22 00:00:00 +00:00' # if None, xarray uses ds.time.encoding['units'] as refdate_str
dir_output = './test_interpolate_nc_to_bc'

#quantities should be in conversion_dict.keys(). waterlevelbnd is steric/zos, tide is tidal components from FES/EOT
list_quantities = ['waterlevelbnd','salinitybnd','temperaturebnd','uxuy','tracerbndNO3','tide']
#list_quantities = ['waterlevelbnd','salinitybnd','temperaturebnd','tracerbndNO3']
list_quantities = ['salinitybnd','tracerbndNO3']

model = 'CMEMS' #CMEMS GFDL CMCC HYCOM

#The {ncvarname} wildcard in dir_pattern_hydro/dir_patern_waq is used to replace it with conversion_dict[quantity]['ncvarname'] by using str(dir_pattern).format(ncvarname)
if model=='CMEMS': #2012-01-06 12:00:00 to 2013-01-03 12:00:00
    conversion_dict = dfmt.get_conversion_dict()
    tstart = '2012-01-16 12:00'
    tstop = '2012-04-01 12:00'
    #tstop = '2013-01-01 12:00'
    list_plifiles = [Path(r'p:\11208054-004-dcsm-fm\models\model_input\bnd_cond\pli\DCSM-FM_OB_all_20181108_nocomments.pli')] #TODO: reading this file without '_nocomments' results in empty Polyfile, should raise an error. https://github.com/Deltares/HYDROLIB-core/issues/320
    dir_sourcefiles_hydro = r'p:\1204257-dcsmzuno\data\CMEMS\nc\DCSM_allAvailableTimes' #CMEMS hydro: bottomT, so, thetao, uo, vo, zos (2012-01-06 12:00:00 to 2013-01-03 12:00:00) (daily values at noon, not at midnight)
    dir_pattern_hydro = Path(dir_sourcefiles_hydro,'{ncvarname}_2012*.nc') # later remove 2012 from string, but this is faster for testing #TODO: it is quite slow, maybe speed up possible?
    dir_sourcefiles_waq = r'p:\11206304-futuremares\python_scripts\ocean_boundaryCMEMS\data_monthly' #CMEMS waq: no3, o2, phyc, so4, si (2011-12-16 12:00:00 to 2019-01-16 12:00:00)
    dir_pattern_waq = Path(dir_sourcefiles_waq,'cmems_mod_glo_bgc_my_0.25_P1M-m_{ncvarname}_*.nc') 
    #to reproduce old CMEMS data (icw reverse_depth=True) (from p:\1204257-dcsmzuno\data\CMEMS\bnd\NorthSeaAndBaltic_1993-2019_20210510)
    #tstart = dt.datetime(1993,1,1,12,0)
    #tstop = tstart + dt.timedelta(days=5)
    #dir_pattern_hydro = Path(dir_sourcefiles_hydro,'{ncvarname}_1993*.nc')
elif model=='GFDL':
    conversion_dict = dfmt.get_conversion_dict()
    tstart = '2012-01-16 12:00'
    tstop = '2012-04-01 12:00'
    list_plifiles = [Path(r'p:\11208054-004-dcsm-fm\models\model_input\bnd_cond\pli\DCSM-FM_OB_all_20181108_nocomments.pli')] #TODO: reading this file without '_nocomments' results in empty Polyfile, should raise an error. https://github.com/Deltares/HYDROLIB-core/issues/320
    dir_sourcefiles_hydro = None
    dir_pattern_hydro = None
    dir_sourcefiles_waq = r'p:\11206304-futuremares\data\CMIP6_BC\GFDL-ESM4' #GFDL waq: no3 (1850-01-16 12:00:00 to 2014-12-16 12:00:00)
    dir_pattern_waq = Path(dir_sourcefiles_waq,'{ncvarname}_esm-hist.nc')
elif model=='CMCC': #TODO: check method, now finding nearest points (so always has values)
    #TODO: time_bnds/lev_bnds are available, take into account in bc file?
    conversion_dict = dfmt.get_conversion_dict(ncvarname_updates={'salinitybnd':'sos', 'temperaturebnd':'tos'})
    conversion_dict['tracerbndNO3'] = {'ncvarname':'no3', 'unit':'g/m3', 'conversion':14.0} #other vars also have different conversion than cmems
    tstart = '2015-06-16 12:00'
    tstop = '2015-12-01 12:00'
    list_plifiles = [Path(r'p:\11208054-004-dcsm-fm\models\model_input\bnd_cond\pli\DCSM-FM_OB_all_20181108_nocomments.pli')] #TODO: reading this file without '_nocomments' results in empty Polyfile, should raise an error. https://github.com/Deltares/HYDROLIB-core/issues/320
    dir_sourcefiles_hydro = r'p:\11206304-futuremares\data\CMIP6_BC\CMCC-ESM2'
    dir_pattern_hydro = Path(dir_sourcefiles_hydro,'{ncvarname}_Omon_CMCC-ESM2_ssp126_r1i1p1f1_gn_*.nc')
    dir_sourcefiles_waq = dir_sourcefiles_hydro #CMCC waq: (2015-01-16 12:00:00 to 2100-12-16 12:00:00)
    dir_pattern_waq = dir_pattern_hydro
elif model=='HYCOM':
    conversion_dict = dfmt.get_conversion_dict(ncvarname_updates={'salinitybnd':'salinity', 'temperaturebnd':'water_temp'})
    tstart = '2016-04-20'
    tstop = '2016-05-03'
    list_plifiles = [Path(r'c:\DATA\dfm_tools_testdata\GLBu0.08_expt_91.2\bcline.pli')] #HYCOM not available in DCSM area, so use other pli-file
    dir_sourcefiles_hydro = 'c:\\DATA\\dfm_tools_testdata\\GLBu0.08_expt_91.2' #HYCOM hydro: salinity/so, water_temp/thetao (2016-04-19 00:00:00 to 2016-05-06 00:00:00)
    dir_pattern_hydro = Path(dir_sourcefiles_hydro,'HYCOM_ST_GoO_*.nc')
    dir_sourcefiles_waq = None
    dir_pattern_waq = None
else:
    raise Exception(f'invalid model: {model}')


# start of interpolation process
dtstart = dt.datetime.now()
ext_bnd = ExtModel()
if not os.path.isdir(dir_output):
    os.mkdir(dir_output)

for file_pli in list_plifiles:
    file_bc_basename = file_pli.name.replace('.pli','')
    for quantity in list_quantities:
        print(f'processing quantity: {quantity}')
        if quantity=='tide': #TODO: choose flexible/generic component notation
            tidemodel = 'FES2014' #FES2014, FES2012, EOT20
            component_list = ['2n2','mf','p1','m2','mks2','mu2','q1','t2','j1','m3','mm','n2','r2','k1','m4','mn4','s1','k2','m6','ms4','nu2','s2','l2','m8','msf','o1','s4'] #None results in all FES components
            ForcingModel_object = dfmt.interpolate_tide_to_bc(tidemodel=tidemodel, file_pli=file_pli, component_list=component_list, nPoints=nPoints)
            for forcingobject in ForcingModel_object.forcing: #add A0 component
                forcingobject.datablock.append(['A0',0.0,0.0])
        else:
            if quantity in ['waterlevelbnd','salinitybnd','temperaturebnd','uxuy']: #hydro
                if dir_sourcefiles_hydro is None:
                    continue
                if (model=='HYCOM') & (quantity not in ['salinitybnd','temperaturebnd']): #only contains quantities salinity and water_temp, so crashes on others
                    continue
                dir_pattern = dir_pattern_hydro
            else: #waq
                if dir_pattern_waq is None:
                    continue
                dir_pattern = dir_pattern_waq
            
            #open regulargridDataset and do some basic stuff (time selection, renaming depth/lat/lon/varname, converting units, etc)
            data_xr_vars = dfmt.open_dataset_extra(dir_pattern=dir_pattern, quantity=quantity, #TODO: maybe replace renaming part with package CMCC/Lisa?
                                                   tstart=tstart, tstop=tstop,
                                                   conversion_dict=conversion_dict,
                                                   refdate_str=refdate_str)
            #interpolate regulargridDataset to plipointsDataset
            data_interp = dfmt.interp_regularnc_to_plipoints(data_xr_reg=data_xr_vars, file_pli=file_pli,
                                                             nPoints=nPoints) #argument for testing
            #data_interp = data_interp.ffill(dim="plipoints").bfill(dim="plipoints") #to fill allnan plipoints with values from the neighbour point
            
            #convert plipointsDataset to hydrolib ForcingModel
            ForcingModel_object = dfmt.plipointsDataset_to_ForcingModel(plipointsDataset=data_interp)
                    
        file_bc_basename = file_pli.name.replace('.pli','')
        if quantity=='tide':
            file_bc_out = Path(dir_output,f'{quantity}_{file_bc_basename}_{tidemodel}.bc')
        else:
            file_bc_out = Path(dir_output,f'{quantity}_{file_bc_basename}_{model}.bc')
        
        print(f'writing ForcingModel to bc file with hydrolib ({file_bc_out.name})')
        bc_type = 'bc' #TODO: add netcdf bc support. https://github.com/Deltares/HYDROLIB-core/issues/318
        if bc_type=='bc':
            #ForcingModel_object.serializer_config.float_format = '.3f' #TODO SOLVED: improve formatting of bc file: https://github.com/Deltares/HYDROLIB-core/issues/308
            #ForcingModel_object.serializer_config.float_format_datablock = '.5f' #maybe move this to interp_regularnc_to_plipoints/interpolate_tide_to_bc?
            ForcingModel_object.save(filepath=file_bc_out)
        
        #TODO: support for relative paths?
        #generate boundary object for the ext file (quantity, pli-filename, bc-filename)
        boundary_object = Boundary(quantity=quantity.replace('tide','waterlevelbnd'), #the FM quantity for tide is also waterlevelbnd
                                   locationfile=file_pli,
                                   forcingfile=ForcingModel_object)
        ext_bnd.boundary.append(boundary_object)

        if 1 and quantity!='tide': #TODO: data_xr_vars/data_interp does not exist for tide yet
            #plotting dataset and polyline (is wrong for CMCC)
            varname0 = list(data_xr_vars.data_vars)[0] 
            fig,ax = plt.subplots()
            if 'depth' in data_xr_vars[varname0].dims:
                data_xr_vars[varname0].isel(time=0,depth=0).plot(ax=ax)
            else:
                data_xr_vars[varname0].isel(time=0).plot(ax=ax)
            plipoint_coords = data_interp.plipoints.to_dataframe()
            ax.plot(plipoint_coords['plipoint_x'],plipoint_coords['plipoint_y'],'r-')
            ctx.add_basemap(ax=ax,crs="EPSG:4326",attribution=False)
            fig.tight_layout()
            fig.savefig(str(file_bc_out).replace('.bc','_polyline'))
            
            #plotting example data point
            for iF in [2]:#range(nPoints): 
                data_vars = list(data_interp.data_vars)
                fig,ax1 = plt.subplots(figsize=(10, 6))
                data_interp[data_vars[0]].isel(plipoints=iF).T.plot()
                fig.tight_layout()
                fig.savefig(str(file_bc_out).replace('.bc',''))

file_ext_out = Path(dir_output,'example_bnd.ext')
ext_bnd.save(filepath=file_ext_out)

time_passed = (dt.datetime.now()-dtstart).total_seconds()
print(f'>>total script time passed: {time_passed:.2f} sec')


