# -*- coding: utf-8 -*-
"""
Created on Thu Aug 18 16:37:38 2022

@author: veenstra

This script can be used to interpolate pre-downloaded CMEMS data (or other netcdf files) to a boundary locations in a pli-file
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
    from hydrolib.core.io.dflowfm.ext.models import Boundary, ExtModel

#TODO: put all functions under dfm_tools init for more convenience
#TODO: discuss structure with arthur/stendert (new modularity, naming of functions)
#TODO: make both interpolate functions accept tstart/tstop as datestrings
#TODO: add coordinate conversion of pli-coordinates (for nesting RD models)
#copied plifile from DCSM folder: r'p:\1204257-dcsmzuno\data\CMEMS\bnd\NorthSeaAndBaltic_1993-2019_20210510'
#list_plifiles = [Path(r'c:\DATA\dfm_tools_testdata\hydrolib_bc\DCSM\DCSM-FM_OB_all_20181108.pli')] #TODO: reading this file results in empty Polyfile, should raise an error. https://github.com/Deltares/HYDROLIB-core/issues/320
list_plifiles = [Path(r'c:\DATA\dfm_tools_testdata\hydrolib_bc\DCSM\DCSM-FM_OB_all_20181108_nocomments.pli')]
nPoints = None #amount of Points to process per PolyObject in the plifile (use int for testing, use None for all Points)

dir_output = './test_interpolate_nc_to_bc'
bc_type = 'bc' #currently only 'bc' supported #TODO: add netcdf bc support. https://github.com/Deltares/HYDROLIB-core/issues/318

refdate_str = 'minutes since 2011-12-22 00:00:00 +00:00' # if None, xarray uses ds.time.encoding['units'] as refdate_str

#quantities should be in conversion_dict.keys(). waterlevelbnd is steric/zos, tide is tidal components from FES/EOT
list_quantities = ['waterlevelbnd','salinitybnd','temperaturebnd','uxuy','tracerbndNO3']#,'tide']
#list_quantities = ['waterlevelbnd','salinitybnd','temperaturebnd','tracerbndNO3']
list_quantities = ['tracerbndNO3']

model = 'CMEMS' #CMEMS GFDL CMCC HYCOM

#The {ncvarname} wildcard in dir_pattern_hydro/dir_patern_waq is used to replace it with conversion_dict[quantity]['ncvarname'] by using str(dir_pattern).format(ncvarname)
reverse_depth = False #to compare with coastserv files, this argument will be phased out
kdtree_k = 3
if model=='CMEMS': #2012-01-06 12:00:00 to 2013-01-03 12:00:00
    conversion_dict = dfmt.get_conversion_dict()
    tstart = dt.datetime(2012, 1, 16, 12, 0)
    tstop = dt.datetime(2012, 4, 1, 12, 0)
    #tstop = dt.datetime(2013, 1, 1, 12, 0)
    dir_sourcefiles_hydro = r'p:\1204257-dcsmzuno\data\CMEMS\nc\DCSM_allAvailableTimes' #CMEMS hydro: bottomT, so, thetao, uo, vo, zos (2012-01-06 12:00:00 to 2013-01-03 12:00:00) (daily values at noon, not at midnight)
    dir_pattern_hydro = Path(dir_sourcefiles_hydro,'{ncvarname}_2012*.nc') # later remove 2012 from string, but this is faster for testing #TODO: it is quite slow, maybe speed up possible?
    dir_sourcefiles_waq = r'p:\11206304-futuremares\python_scripts\ocean_boundaryCMEMS\data_monthly' #CMEMS waq: no3, o2, phyc, so4, si (2011-12-16 12:00:00 to 2019-01-16 12:00:00)
    dir_pattern_waq = Path(dir_sourcefiles_waq,'cmems_mod_glo_bgc_my_0.25_P1M-m_{ncvarname}_*.nc') 
    #to reproduce old CMEMS data (icw reverse_depth=True)
    #tstart = dt.datetime(1993,1,1,12,0)
    #tstop = tstart + dt.timedelta(days=5)
    #dir_pattern_hydro = Path(dir_sourcefiles_hydro,'{ncvarname}_1993*.nc')
elif model=='GFDL':
    conversion_dict = dfmt.get_conversion_dict()
    tstart = dt.datetime(2012, 1, 16, 12, 0)
    tstop = dt.datetime(2012, 4, 1, 12, 0)
    dir_sourcefiles_hydro = None
    dir_pattern_hydro = None
    dir_sourcefiles_waq = r'p:\11206304-futuremares\data\CMIP6_BC\GFDL-ESM4' #GFDL waq: no3 (1850-01-16 12:00:00 to 2014-12-16 12:00:00)
    dir_pattern_waq = Path(dir_sourcefiles_waq,'{ncvarname}_esm-hist.nc')
elif model=='CMCC': #TODO: check method, now finding nearest points (so always has values)
    #TODO: time_bnds/lev_bnds are available, take into account in bc file?
    conversion_dict = dfmt.get_conversion_dict(ncvarname_updates={'salinitybnd':'sos', 'temperaturebnd':'tos'})
    conversion_dict['tracerbndNO3'] = {'ncvarname':'no3', 'unit':'g/m3', 'conversion':14.0} #other vars also have different conversion than cmems
    tstart = dt.datetime(2015, 6, 16, 12, 0)
    tstop = dt.datetime(2016, 12, 1, 12, 0)
    dir_sourcefiles_hydro = r'p:\11206304-futuremares\data\CMIP6_BC\CMCC-ESM2'
    dir_pattern_hydro = Path(dir_sourcefiles_hydro,'{ncvarname}_Omon_CMCC-ESM2_ssp126_r1i1p1f1_gn_*.nc')
    dir_sourcefiles_waq = dir_sourcefiles_hydro #CMCC waq: (2015-01-16 12:00:00 to 2100-12-16 12:00:00)
    dir_pattern_waq = dir_pattern_hydro
    KDTree_invdistweigth = True #only relevant for CMCC
elif model=='HYCOM':
    conversion_dict = dfmt.get_conversion_dict(ncvarname_updates={'salinitybnd':'salinity', 'temperaturebnd':'water_temp'})
    tstart = dt.datetime(2016, 4, 20, 0, 0) #HYCOM
    tstop = dt.datetime(2016, 5, 3, 0, 0)
    dir_sourcefiles_hydro = 'c:\\DATA\\dfm_tools_testdata\\GLBu0.08_expt_91.2' #HYCOM hydro: salinity/so, water_temp/thetao (2016-04-19 00:00:00 to 2016-05-06 00:00:00)
    dir_pattern_hydro = Path(dir_sourcefiles_hydro,'HYCOM_ST_GoO_*.nc')
    dir_sourcefiles_waq = None
    dir_pattern_waq = None
    list_plifiles = [Path(r'c:\DATA\dfm_tools_testdata\GLBu0.08_expt_91.2\bcline.pli')] #HYCOM not available in DCSM area, so use other pli-file
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
                                                   refdate_str=refdate_str, 
                                                   reverse_depth=reverse_depth) #temporary argument to compare easier with old coastserv files            
            #interpolate regulargridDataset to plipointsDataset
            data_interp = dfmt.interp_regularnc_to_plipoints(data_xr_reg=data_xr_vars, file_pli=file_pli,
                                                             nPoints=nPoints, #argument for testing
                                                             kdtree_k=kdtree_k) #argument for curvi grids like CMCC
            #data_interp = data_interp.ffill(dim="plipoints").bfill(dim="plipoints") #to fill allnan plipoints with values from the neighbour point
            
            #convert plipointsDataset to hydrolib ForcingModel
            ForcingModel_object = dfmt.plipointsDataset_to_ForcingModel(plipointsDataset=data_interp)
                    
        file_bc_basename = file_pli.name.replace('.pli','')
        if quantity=='tide':
            file_bc_out = Path(dir_output,f'{quantity}_{file_bc_basename}_{tidemodel}.bc')
        else:
            file_bc_out = Path(dir_output,f'{quantity}_{file_bc_basename}_{model}.bc')
        print(f'writing ForcingModel to bc file with hydrolib ({file_bc_out.name})')
        if bc_type=='bc':
            #ForcingModel_object.serializer_config.float_format = '.3f'
            #ForcingModel_object.serializer_config.float_format_datablock = '.5f'
            ForcingModel_object.save(filepath=file_bc_out)
            #TODO SOLVED: improve formatting of bc file: https://github.com/Deltares/HYDROLIB-core/issues/308 (became quite slow: https://github.com/Deltares/HYDROLIB-core/issues/313)
        else:
            raise Exception(f'invalid bc_type: {bc_type}')
        
        #make paths relative (sort of) (also necessary for locationfile) /../ should also be supported? 
        #ForcingModel_object.filepath = Path(str(ForcingModel_object.filepath).replace(dir_out,'')) #TODO: convert to relative paths in ext file possible? This path is the same as file_bc_out
        
        #generate boundary object for the ext file (quantity, pli-filename, bc-filename)
        boundary_object = Boundary(quantity=quantity.replace('tide','waterlevelbnd'), #the FM quantity for tide is also waterlevelbnd
                                   locationfile=file_pli,
                                   forcingfile=ForcingModel_object,
                                   )
        ext_bnd.boundary.append(boundary_object)

        if 1 and quantity!='tide': #plotting dataset and polyline (is wrong for CMCC)
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
            
        if 1: #plotting example data point
            # for iF in [2]:#range(nPoints): #from data_interp TODO: this is slightly easier, but negative depth is not done yet (is in Dataset_to_T3D/Dataset_to_TimeSeries functions)
            #     data_vars = list(data_interp.data_vars)
            #     fig,ax1 = plt.subplots(figsize=(10, 6))
            #     data_interp[data_vars[0]].isel(plipoints=iF).T.plot()
            #     fig.tight_layout()
            for iF in [2]:#range(nPoints): #from ForcingModel_object
                forcingobject_one = ForcingModel_object.forcing[iF]
                forcingobject_one_xr = dfmt.forcinglike_to_Dataset(forcingobject_one,convertnan=True)
                data_vars = list(forcingobject_one_xr.data_vars)
                fig,ax1 = plt.subplots(figsize=(10, 6))
                if hasattr(forcingobject_one,'vertpositions'): #3D quantity (time/depth dimensions)
                    if hasattr(forcingobject_one.quantityunitpair[1],'elementname'): #uxuy vector
                        plt.close(fig)
                        fig, axes = plt.subplots(2,1,figsize=(10, 6),sharex=True,sharey=True)
                        forcingobject_one_xr[data_vars[0]].T.plot(ax=axes[0])
                        forcingobject_one_xr[data_vars[1]].T.plot(ax=axes[1])
                    else:
                        forcingobject_one_xr[data_vars[0]].T.plot(ax=ax1)
                elif quantity=='tide':
                    ax2 = ax1.twinx()
                    forcingobject_one_xr[data_vars[0]].plot(ax=ax1)
                    forcingobject_one_xr[data_vars[1]].plot(ax=ax2)
                else:
                    forcingobject_one_xr[data_vars[0]].plot(ax=ax1)
                fig.tight_layout()
                fig.savefig(str(file_bc_out).replace('.bc',''))


file_ext_out = Path(dir_output,'example_bnd.ext')
ext_bnd.save(filepath=file_ext_out)

time_passed = (dt.datetime.now()-dtstart).total_seconds()
print(f'>>total script time passed: {time_passed:.2f} sec')


