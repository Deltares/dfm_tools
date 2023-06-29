"""
Created on Tue Apr  4 16:12:56 2023

@author: groenenb, sclaan
edited by: veenstra, zijlker

This is a proof of concept that will be properly coded in the near future and probably end up under hydromt-delft3dfm
Since the functions in this script contain hardcoded parameters, it is not exposed to public and you need to import like dfmt.modelbuilder.[function]
"""

import os
import xarray as xr
import pandas as pd
import dfm_tools as dfmt
import hydrolib.core.dflowfm as hcdfm
import datetime as dt
import glob


def preprocess_interpolate_nc_to_bc(ext_bnd,
                                    list_quantities,#quantities should be in conversion_dict.keys(). waterlevelbnd is steric/zos, tide is tidal components from FES/EOT
                                    tstart,
                                    tstop,
                                    list_plifiles,
                                    dir_pattern,
                                    dir_output='.',
                                    refdate_str=None,
                                    conversion_dict=dfmt.get_conversion_dict(),
                                    tidemodel='FES2014'): #FES2014, FES2012, EOT20, GTSM4.1preliminary
    #input examples in https://github.com/Deltares/dfm_tools/blob/main/tests/examples/preprocess_interpolate_nc_to_bc.py
    model = 'CMEMS'
    
    for file_pli in list_plifiles:
        file_bc_basename = os.path.basename(file_pli).replace('.pli','')
        for quantity in list_quantities:
            print(f'processing quantity: {quantity}')
            if quantity=='tide': 
                if tidemodel == 'FES2014': #for comparing to older FES bc-files
                    component_list = ['2n2','la2','mf','mtm','p1','ssa','eps2','m2','mks2','mu2','q1','t2','j1','m3','mm','n2','r2','k1','m4','mn4','n4','s1','k2','m6','ms4','nu2','s2','l2','m8','msf','o1','s4','sa','msqm']
                else:
                    component_list = None #None results in all tidemodel components
                ForcingModel_object = dfmt.interpolate_tide_to_bc(tidemodel=tidemodel, file_pli=file_pli, component_list=component_list, nPoints=None)
                ForcingModel_object.serializer_config.float_format_datablock = "21.9e"
            else:
                tstart, tstop = dfmt.round_timestamp_to_outer_noon(tstart,tstop)
                #open regulargridDataset and do some basic stuff (time selection, renaming depth/lat/lon/varname, converting units, etc)
                data_xr_vars = dfmt.open_dataset_extra(dir_pattern=dir_pattern, quantity=quantity,
                                                       tstart=tstart, tstop=tstop,
                                                       conversion_dict=conversion_dict,
                                                       refdate_str=refdate_str)
                #interpolate regulargridDataset to plipointsDataset
                data_interp = dfmt.interp_regularnc_to_plipoints(data_xr_reg=data_xr_vars, file_pli=file_pli)
                
                #convert plipointsDataset to hydrolib ForcingModel
                ForcingModel_object = dfmt.plipointsDataset_to_ForcingModel(plipointsDataset=data_interp)
                        
            if quantity=='tide':
                file_bc_out = os.path.join(dir_output,f'{quantity}_{file_bc_basename}_{tidemodel}.bc')
            else:
                file_bc_out = os.path.join(dir_output,f'{quantity}_{file_bc_basename}_{model}.bc')
            
            ForcingModel_object.save(filepath=file_bc_out)
            
            #generate boundary object for the ext file (quantity, pli-filename, bc-filename)
            boundary_object = hcdfm.Boundary(quantity=quantity.replace('tide','waterlevelbnd'), #the FM quantity for tide is also waterlevelbnd
                                             locationfile=file_pli,
                                             forcingfile=ForcingModel_object)
            ext_bnd.boundary.append(boundary_object)
    
    return ext_bnd

    
def preprocess_ini_cmems_to_nc(ext_old, tstart='1998-01-01',
        dir_data  = r'p:\i1000668-tcoms\03_newModel\01_input\02_bnd\data_opendap', #folder containing CMEMS so and thetao netcdf files
        dir_out = '.'):
    
    file_nc_list_so = glob.glob(f'{dir_data}\\cmems_so_*.nc')
    file_nc_list_thetao = glob.glob(f'{dir_data}\\cmems_thetao_*.nc')
    file_nc_list = file_nc_list_so + file_nc_list_thetao
    
    print(f'opening {len(file_nc_list)} datasets')
    data_xr = xr.open_mfdataset(file_nc_list)
    
    tSimStart = pd.Timestamp(tstart)
    data_xr_ontime = data_xr.sel(time=slice(tSimStart-dt.timedelta(days=1),tSimStart+dt.timedelta(days=1)))
    
    print('writing file')
    outFile = os.path.join(dir_out,f'InitialField_{tSimStart.strftime("%Y-%m-%d_%H-%M-%S")}.nc')
    data_xr_ontime.to_netcdf(outFile,format="NETCDF4_CLASSIC") #TODO: why the format setting?
    
    #append forcings to ext
    # 3D initialsalinity/initialtemperature fields are silently ignored, initial 3D conditions are only possible via nudging 1st timestep via quantity=nudge_salinity_temperature
    forcing_saltem = hcdfm.ExtOldForcing(quantity='nudge_salinity_temperature',
                                         filename=outFile,
                                         filetype=hcdfm.ExtOldFileType.NetCDFGridData,
                                         method=hcdfm.ExtOldMethod.InterpolateTimeAndSpaceSaveWeights, #3
                                         operand=hcdfm.Operand.override, #O
                                         )
    ext_old.forcing.append(forcing_saltem)
    
    return ext_old
    
    
def preprocess_merge_meteofiles(ext_old,
        mode = 'HYCOM', # 'HIRLAM_meteo' 'HIRLAM_meteo-heatflux' 'HARMONIE' 'HYCOM' 'ERA5_wind_pressure' 'ERA5_heat_model' 'ERA5_radiation' 'ERA5_rainfall' 'WOA'
        varkey_list = [],
        dir_data = [],
        dir_output = '.',
        time_slice = slice('2013-12-30','2014-01-01')):

    if isinstance(varkey_list[0], list):
        varkey_lists = varkey_list
    else:
        varkey_lists = [varkey_list]
    
    for varkey_list in varkey_lists:
        if 'HIRLAM' in mode:
            if mode == 'HIRLAM_meteo': #1year voor meteo crasht (HIRLAM72_*\\h72_*) door conflicting dimension sizes, sourcefolders opruimen? meteo_heatflux folders zijn schoner dus daar werkt het wel
                dir_data = 'P:\\1204257-dcsmzuno\\2014\\data\\meteo\\HIRLAM72_*' #files contain: ['air_pressure_fixed_height','northward_wind','eastward_wind']
            elif mode == 'HIRLAM_meteo-heatflux':
                dir_data = 'P:\\1204257-dcsmzuno\\2014\\data\\meteo-heatflux\\HIRLAM72_*' # files contain: ['dew_point_temperature','air_temperature','cloud_area_fraction']
            fn_match_pattern = 'h72_20131*.nc'
            file_out_prefix = 'h72_'
            preprocess = dfmt.preprocess_hirlam #temporary(?) fix for >1D-vars with same name as its dim
        elif mode == 'HARMONIE':
            dir_data = 'p:\\1204257-dcsmzuno\\data\\meteo\\HARMONIE\\nc\\air_*' #many invalid files, so subsetting here
            fn_match_pattern = 'HARMONIE_*_2020_*.nc'
            file_out_prefix = 'HARMONIE_'
            preprocess = None
        elif mode == 'HYCOM':
            dir_data = 'c:\\DATA\\dfm_tools_testdata\\GLBu0.08_expt_91.2'
            fn_match_pattern = 'HYCOM_ST_GoO_*.nc'
            file_out_prefix = 'HYCOM_ST_GoO_'
            preprocess = None
            #rename_variables = {'salinity':'so', 'water_temp':'thetao'}
        elif 'ERA5' in mode:
            fn_match_pattern = f'era5_.*({"|".join(varkey_list)})_.*.nc' #simpler but selects more files: 'era5_*.nc'
            file_out_prefix = f'era5_{"_".join(varkey_list)}_'
            preprocess = dfmt.preprocess_ERA5 #reduce expver dimension if present
        elif mode == 'WOA':
            dir_data = r'p:\1204257-dcsmzuno\data\WOA13'
            fn_match_pattern = 'woa13_decav_s*.nc'
            file_out_prefix = 'woa13_decav_s_'
            preprocess = dfmt.preprocess_woa #add 360-day calendar unit to time attrs before decode_cf
        else:
            raise Exception('ERROR: wrong mode %s'%(mode))
        
        if not os.path.exists(dir_output):
            os.makedirs(dir_output)
        
        file_nc = os.path.join(dir_data,fn_match_pattern)
        
        data_xr_tsel = dfmt.merge_meteofiles(file_nc=file_nc, time_slice=time_slice, 
                                             preprocess=preprocess,
                                             add_global_overlap=False, #GTSM specific: extend data beyond -180 to 180 longitude
                                             zerostart=False) #GTSM specific: extend data with 0-value fields 1 and 2 days before all_tstart
        
        #write to netcdf file
        print('>> writing file (can take a while): ',end='')
        dtstart = dt.datetime.now()
        try:
            times_np = data_xr_tsel['time'].to_series()
        except:
            times_np = data_xr_tsel['time'].to_numpy() #.to_series() does not work for woa 360_day data. .to_numpy() results in numpy.datetime64 for other datasets, which has no attribute strftime
        time_start_str = times_np[0].strftime("%Y%m%d")
        time_stop_str = times_np[-1].strftime("%Y%m%d")
        file_out = os.path.join(dir_output, f'{file_out_prefix}{time_start_str}to{time_stop_str}_{mode}.nc')
        data_xr_tsel.to_netcdf(file_out)
        print(f'{(dt.datetime.now()-dtstart).total_seconds():.2f} sec')
        
        #append to ext model
        if varkey_list == ['msl','u10n','v10n','chnk']:
            forcing_meteo = hcdfm.ExtOldForcing(quantity='airpressure_windx_windy_charnock',
                                                filename=file_out,
                                                varname='msl u10n v10n chnk',
                                                filetype=hcdfm.ExtOldFileType.NetCDFGridData, #11
                                                method=hcdfm.ExtOldMethod.InterpolateTimeAndSpaceSaveWeights, #3
                                                operand=hcdfm.Operand.override, #O
                                                )
            ext_old.forcing.append(forcing_meteo)
        elif varkey_list == ['d2m','t2m','tcc']:
            forcing_meteo = hcdfm.ExtOldForcing(quantity='dewpoint_airtemperature_cloudiness',
                                                filename=file_out,
                                                varname='d2m t2m tcc',
                                                filetype=hcdfm.ExtOldFileType.NetCDFGridData, #11
                                                method=hcdfm.ExtOldMethod.InterpolateTimeAndSpaceSaveWeights, #3
                                                operand=hcdfm.Operand.override, #O
                                                )
            ext_old.forcing.append(forcing_meteo)
        elif varkey_list == ['ssr','strd']:
            forcing_meteo = hcdfm.ExtOldForcing(quantity='solarradiation',
                                                filename=file_out,
                                                varname='ssr',
                                                filetype=hcdfm.ExtOldFileType.NetCDFGridData, #11
                                                method=hcdfm.ExtOldMethod.InterpolateTimeAndSpaceSaveWeights, #3
                                                operand=hcdfm.Operand.override, #O
                                                )
            ext_old.forcing.append(forcing_meteo)
            forcing_meteo = hcdfm.ExtOldForcing(quantity='longwaveradiation',
                                                filename=file_out,
                                                varname='strd',
                                                filetype=hcdfm.ExtOldFileType.NetCDFGridData, #11
                                                method=hcdfm.ExtOldMethod.InterpolateTimeAndSpaceSaveWeights, #3
                                                operand=hcdfm.Operand.override, #O
                                                )
            ext_old.forcing.append(forcing_meteo)
        elif varkey_list == ['mer','mtpr']:
            forcing_meteo = hcdfm.ExtOldForcing(quantity='rainfall_rate',
                                                filename=file_out,
                                                varname='mtpr',
                                                filetype=hcdfm.ExtOldFileType.NetCDFGridData, #11
                                                method=hcdfm.ExtOldMethod.InterpolateTimeAndSpaceSaveWeights, #3
                                                operand=hcdfm.Operand.override, #O
                                                )
            ext_old.forcing.append(forcing_meteo)
            forcing_meteo = hcdfm.ExtOldForcing(quantity='rainfall_rate',
                                                filename=file_out,
                                                varname='mer',
                                                filetype=hcdfm.ExtOldFileType.NetCDFGridData, #11
                                                method=hcdfm.ExtOldMethod.InterpolateTimeAndSpaceSaveWeights, #3
                                                operand=hcdfm.Operand.add, #+
                                                )
            ext_old.forcing.append(forcing_meteo)
        
    return ext_old

