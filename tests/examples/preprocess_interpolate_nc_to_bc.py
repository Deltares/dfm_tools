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
import hydrolib.core.dflowfm as hcdfm

#TODO: add coordinate conversion of pli-coordinates? (for nesting RD models in oceanmodels)
#TODO: additional models/sources for download/interpolate (evt xESMF for CMCC, climate forcing cmip6 procedure (=calendarconversion) and others)

nPoints = 3# None #amount of Points to process per PolyObject in the plifile (use int for testing, use None for all Points)
refdate_str = 'minutes since 2011-12-22 00:00:00 +00:00' # if None, xarray uses ds.time.encoding['units'] as refdate_str
dir_output = './test_interpolate_nc_to_bc_TEMP'

#quantities should be in conversion_dict.keys(). waterlevelbnd is steric/zos, tide is tidal components from FES/EOT
list_quantities = ['waterlevelbnd','salinitybnd','temperaturebnd','uxuyadvectionvelocitybnd','tracerbndNO3','tide']
#list_quantities = ['waterlevelbnd','salinitybnd','temperaturebnd','tracerbndNO3']
list_quantities = ['salinitybnd','tracerbndNO3','tide']
list_quantities = ['salinitybnd','tracerbndNO3']
#list_quantities = ['waterlevelbnd','salinitybnd','temperaturebnd','uxuyadvectionvelocitybnd','tracerbndNO3','tracerbndOpal','tracerbndDON','tide'] #also waq vars with same ncvarname, opal not available for GFDL and CMCC

model = 'CMEMS' #CMEMS GFDL

#The {ncvarname} wildcard in dir_pattern_hydro/dir_patern_waq is used to replace it with conversion_dict[quantity]['ncvarname'] by using str(dir_pattern).format(ncvarname)
list_plifiles = [Path(r'p:\archivedprojects\11208054-004-dcsm-fm\models\model_input\bnd_cond\pli\DCSM-FM_OB_all_20181108.pli')]
if model=='CMEMS': #2012-01-06 12:00:00 to 2013-01-03 12:00:00
    conversion_dict = dfmt.get_conversion_dict()
    tstart = '2012-01-16 12:00'
    tstop = '2012-04-01 12:00'
    #tstop = '2013-01-01 12:00'
    dir_sourcefiles_hydro = r'p:\1204257-dcsmzuno\data\CMEMS\nc\DCSM_allAvailableTimes' #CMEMS hydro: bottomT, so, thetao, uo, vo, zos (2012-01-06 12:00:00 to 2013-01-03 12:00:00) (daily values at noon, not at midnight)
    dir_pattern_hydro = Path(dir_sourcefiles_hydro,'{ncvarname}_2012*.nc') # later remove 2012 from string, but this is faster for testing #TODO: it is quite slow, maybe speed up possible?
    dir_sourcefiles_waq = r'p:\archivedprojects\11206304-futuremares-rawdata-preps\python_scripts\ocean_boundaryCMEMS\data_monthly' #CMEMS waq: no3, o2, phyc, so4, si (2011-12-16 12:00:00 to 2019-01-16 12:00:00)
    dir_pattern_waq = Path(dir_sourcefiles_waq,'cmems_mod_glo_bgc_my_0.25_P1M-m_{ncvarname}_*.nc') 
    #to reproduce old CMEMS data (icw reverse_depth=True) (from p:\1204257-dcsmzuno\data\CMEMS\bnd\NorthSeaAndBaltic_1993-2019_20210510)
    #tstart = dt.datetime(1993,1,1,12,0)
    #tstop = tstart + dt.timedelta(days=5)
    #dir_pattern_hydro = Path(dir_sourcefiles_hydro,'{ncvarname}_1993*.nc')
elif model=='GFDL': # TODO: data is not available anymore
    conversion_dict = dfmt.get_conversion_dict()
    tstart = '2012-01-16 12:00'
    tstop = '2012-04-01 12:00'
    dir_sourcefiles_hydro = None
    dir_pattern_hydro = None
    dir_sourcefiles_waq = r'p:\11206304-futuremares\data\CMIP6_BC\GFDL-ESM4' #GFDL waq: no3 (1850-01-16 12:00:00 to 2014-12-16 12:00:00)
    dir_pattern_waq = Path(dir_sourcefiles_waq,'{ncvarname}_esm-hist.nc')
else:
    raise KeyError(f'invalid model: {model}')


# start of interpolation process
dtstart = dt.datetime.now()
ext_bnd = hcdfm.ExtModel()
os.makedirs(dir_output, exist_ok=True)

for file_pli in list_plifiles:
    file_bc_basename = file_pli.name.replace('.pli','')
    for quantity in list_quantities:
        print(f'processing quantity: {quantity}')
        if quantity=='tide': 
            tidemodel = 'FES2014' #FES2014, FES2012, EOT20, GTSMv4.1, GTSMv4.1_opendap
            if tidemodel == 'FES2014': #for comparing to older FES bc-files #TODO: choose flexible/generic component notation
                component_list = ['2n2','mf','p1','m2','mks2','mu2','q1','t2','j1','m3','mm','n2','r2','k1','m4','mn4','s1','k2','m6','ms4','nu2','s2','l2','m8','msf','o1','s4']
            else:
                component_list = None #None results in all tidemodel components
            ForcingModel_object = dfmt.interpolate_tide_to_bc(tidemodel=tidemodel, file_pli=file_pli, component_list=component_list, nPoints=nPoints)
            for forcingobject in ForcingModel_object.forcing: #add A0 component
                forcingobject.datablock.append(['A0',0.0,0.0])
        else:
            if quantity in ['waterlevelbnd','salinitybnd','temperaturebnd','uxuyadvectionvelocitybnd']: #hydro
                if dir_sourcefiles_hydro is None:
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
            # read polyfile as geodataframe
            polyfile_object = hcdfm.PolyFile(file_pli)
            gdf_points_all = dfmt.PolyFile_to_geodataframe_points(polyfile_object)
            gdf_points = gdf_points_all.iloc[:nPoints]
            # interpolate regulargridDataset to plipointsDataset
            data_interp = dfmt.interp_regularnc_to_plipointsDataset(data_xr_reg=data_xr_vars, gdf_points=gdf_points, load=True)
            # data_interp = data_interp.ffill(dim="node").bfill(dim="node") #to fill allnan plipoints with values from the neighbour point #TODO: this also fills the belowbed layers from one point onto another, so should be done after ffill/bfill in depth dimension. Currently all-nan arrays are replaced with .fillna(0)
            
            #convert plipointsDataset to hydrolib ForcingModel
            ForcingModel_object = dfmt.plipointsDataset_to_ForcingModel(plipointsDataset=data_interp)
                    
        bctype = 'bc' #TODO: add netcdf bc support to hcdfm.Boundary: https://github.com/Deltares/HYDROLIB-core/issues/318
        file_bc_basename = file_pli.name.replace('.pli','')
        if quantity=='tide':
            file_bc_out = Path(dir_output,f'{quantity}_{file_bc_basename}_{tidemodel}.{bctype}')
        else:
            file_bc_out = Path(dir_output,f'{quantity}_{file_bc_basename}_{model}.{bctype}')
        if bctype=='bc':
            print(f'writing ForcingModel to bc file with hydrolib ({file_bc_out.name})')
            #ForcingModel_object.serializer_config.float_format = '.3f' #TODO: improve formatting of bc file, maybe move this to interp_regularnc_to_plipoints/interpolate_tide_to_bc?
            #ForcingModel_object.serializer_config.float_format_datablock = '.5f'
            ForcingModel_object.save(filepath=file_bc_out)
        elif bctype=='nc':
            data_interp.to_netcdf(file_bc_out)
        
        #TODO: support for relative paths?
        #generate boundary object for the ext file (quantity, pli-filename, bc-filename)
        boundary_object = hcdfm.Boundary(quantity=quantity.replace('tide','waterlevelbnd'), #the FM quantity for tide is also waterlevelbnd
                                         locationfile=file_pli,
                                         forcingfile=file_bc_out)
        ext_bnd.boundary.append(boundary_object)

        if quantity!='tide': #TODO: data_xr_vars/data_interp does not exist for tide yet
            #plotting dataset and polyline (is wrong for CMCC)
            varname0 = list(data_xr_vars.data_vars)[0] 
            fig,ax = plt.subplots()
            if 'z' in data_xr_vars[varname0].dims:
                data_xr_vars[varname0].isel(time=0,z=0).plot(ax=ax)
            else:
                data_xr_vars[varname0].isel(time=0).plot(ax=ax)
            plipoint_coords = data_interp.node.to_dataframe()
            ax.plot(plipoint_coords['lon'],plipoint_coords['lat'],'r-')
            ctx.add_basemap(ax=ax,crs="EPSG:4326",attribution=False)
            fig.tight_layout()
            fig.savefig(str(file_bc_out).replace('.bc','_polyline'))
            
            #plotting example data point
            for iF in [2]:#range(nPoints):
                data_vars = list(data_interp.data_vars)
                fig,ax1 = plt.subplots(figsize=(10, 6))
                data_interp[data_vars[0]].isel(node=iF).T.plot()
                fig.tight_layout()
                fig.savefig(str(file_bc_out).replace('.bc',''))

file_ext_out = Path(dir_output,'example_bnd.ext')
ext_bnd.save(filepath=file_ext_out)

time_passed = (dt.datetime.now()-dtstart).total_seconds()
print(f'>>total script time passed: {time_passed:.2f} sec')


