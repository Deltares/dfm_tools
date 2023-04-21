"""
Created on Tue Apr  4 16:12:56 2023

@author: groenenb, sclaan
edited by: veenstra

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
import numpy as np

#gridgen
import meshkernel


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
                    component_list = ['2n2','mf','p1','m2','mks2','mu2','q1','t2','j1','m3','mm','n2','r2','k1','m4','mn4','s1','k2','m6','ms4','nu2','s2','l2','m8','msf','o1','s4']
                else:
                    component_list = None #None results in all tidemodel components
                ForcingModel_object = dfmt.interpolate_tide_to_bc(tidemodel=tidemodel, file_pli=file_pli, component_list=component_list, nPoints=None)
            else:
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


#GRIDGENERATION WITH MESHKERNEL

def make_basegrid(lon_min,lon_max,lat_min,lat_max,dx=0.05,dy=0.05,angle=0):
    print('modelbuilder.make_basegrid()')
    # create base grid
    nox = int(np.round((lon_max-lon_min)/dx))
    noy = int(np.round((lat_max-lat_min)/dy))
    
    make_grid_parameters = meshkernel.MakeGridParameters(num_columns=nox,
                                                         num_rows=noy,
                                                         angle=angle,
                                                         origin_x=lon_min,
                                                         origin_y=lat_min,
                                                         block_size_x=dx,
                                                         block_size_y=dy)
    
    geometry_list = meshkernel.GeometryList(np.empty(0, dtype=np.double), np.empty(0, dtype=np.double)) # A polygon must to be provided. If empty it will not be used. If a polygon is provided it will be used in the generation of the curvilinear grid. The polygon must be closed
    mk = meshkernel.MeshKernel() #TODO: is_geographic=True was used in modelbuilder, but refinement super slow and raises "MeshKernelError: MeshRefinement::connect_hanging_nodes: The number of non-hanging nodes is neither 3 nor 4."
    mk.curvilinear_make_uniform(make_grid_parameters, geometry_list) #TODO: make geometry_list argument optional: https://github.com/Deltares/MeshKernelPy/issues/30
    mk.curvilinear_convert_to_mesh2d() #convert to ugrid/mesh2d
    
    return mk


def refine_basegrid(mk, data_bathy_sel,min_face_size=0.1):
    print('modelbuilder.refine_basegrid()')
    samp_x,samp_y = np.meshgrid(data_bathy_sel.lon.to_numpy(),data_bathy_sel.lat.to_numpy())
    samp_z = data_bathy_sel.elevation.to_numpy().astype(float) #TODO: without .astype(float), meshkernelpy generates "TypeError: incompatible types, c_short_Array_27120 instance instead of LP_c_double instance": https://github.com/Deltares/MeshKernelPy/issues/31
    samp_x = samp_x.ravel()
    samp_y = samp_y.ravel()
    samp_z = samp_z.ravel()
    geomlist = meshkernel.GeometryList(x_coordinates=samp_x, y_coordinates=samp_y, values=samp_z) #TODO: does not check if lenghts of input array is equal (samp_z[1:]) https://github.com/Deltares/MeshKernelPy/issues/32
    
    #refinement
    mesh_refinement_parameters = meshkernel.MeshRefinementParameters(refine_intersected=False, #TODO: provide defaults for several arguments, so less arguments are required
                                                                     use_mass_center_when_refining=False, #TODO: what does this do?
                                                                     min_face_size=min_face_size, #TODO: size in meters would be more convenient: https://github.com/Deltares/MeshKernelPy/issues/33
                                                                     refinement_type=meshkernel.RefinementType(1), #Wavecourant/1,
                                                                     connect_hanging_nodes=True, #set to False to do multiple refinement steps (e.g. for multiple regions)
                                                                     account_for_samples_outside_face=True, #outsidecell argument for --refine?
                                                                     max_refinement_iterations=5,
                                                                     ) #TODO: missing the arguments dtmax (necessary?), hmin (min_face_size but then in meters instead of degrees), smoothiters (currently refinement is patchy along coastlines, goes good in dflowfm exec after additional implementation of HK), spherical 1/0 (necessary?)
    
    mk.mesh2d_refine_based_on_samples(samples=geomlist,
                                       relative_search_radius=0.5, #TODO: bilin interp is preferred, but this is currently not supported (samples have to be ravelled): https://github.com/Deltares/MeshKernelPy/issues/34
                                       minimum_num_samples=3,
                                       mesh_refinement_params=mesh_refinement_parameters,
                                       )
    
    return mk


def generate_bndpli(lon_min, lon_max, lat_min, lat_max, dlon, dlat, name='bnd'):

    vals_lon_ar = np.arange(lon_min, lon_max, dlon)
    vals_lon = np.linspace(lon_min, lon_max,len(vals_lon_ar))
    vals_lat_ar = np.arange(lat_min, lat_max, dlat)
    vals_lat = np.linspace(lat_min, lat_max,len(vals_lat_ar))
    pli_p1 = np.c_[np.repeat(lon_min,len(vals_lat)),vals_lat]
    pli_p2 = np.c_[vals_lon,np.repeat(lat_max,len(vals_lon))]
    pli_p3 = np.c_[np.repeat(lon_max,len(vals_lat)),vals_lat[::-1]]
    pli_p4 = np.c_[vals_lon[::-1],np.repeat(lat_min,len(vals_lon))]
    
    pli_all = np.concatenate([pli_p1[:-1],pli_p2[:-1],pli_p3[:-1],pli_p4[:-1]],axis=0)
    
    pli_polyobject = hcdfm.PolyObject(metadata=hcdfm.Metadata(name=name, n_rows=pli_all.shape[0], n_columns=pli_all.shape[1]),
                                      points=[hcdfm.Point(x=x,y=y,data=[]) for x,y in pli_all])
    pli_polyfile = hcdfm.PolyFile(objects=[pli_polyobject])
    return pli_polyfile



