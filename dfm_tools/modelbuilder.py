import os
import logging
import pandas as pd
import dfm_tools as dfmt
import datetime as dt
import xarray as xr
import hydrolib.core.dflowfm as hcdfm
from hydrolib.core.dimr.models import DIMR, FMComponent, Start
from hydrolib.core.utils import get_path_style_for_current_operating_system
from dfm_tools.interpolate_grid2bnd import (ext_add_boundary_object_per_polyline,
                                            open_prepare_dataset,
                                            ds_apply_conversion_dict,
                                            _ds_sel_time_outside,
                                            )
from dfm_tools.xarray_helpers import interpolate_na_multidim

__all__ = [
    "constant_to_bc",
    "cmems_nc_to_bc",
    "cmems_nc_to_ini",
    "preprocess_ini_cmems_to_nc",
    "preprocess_merge_meteofiles_era5",
    "create_model_exec_files",
    "make_paths_relative",
    ]

logger = logging.getLogger(__name__)


def get_quantity_list(quantity):
    if quantity=='uxuyadvectionvelocitybnd': #T3Dvector
        quantity_list = ['ux','uy']
    elif isinstance(quantity, list):
        quantity_list = quantity
    else:
        quantity_list = [quantity]
    return quantity_list


def get_ncvarname(quantity, conversion_dict):
    # check existence of requested keys in conversion_dict
    if quantity not in conversion_dict.keys():
        raise KeyError(f"quantity '{quantity}' not in conversion_dict, (case sensitive) options are"
                       f": {str(list(conversion_dict.keys()))}")
    
    ncvarname = conversion_dict[quantity]['ncvarname']
    return ncvarname


def constant_to_bc(ext_new: hcdfm.ExtModel, file_pli, constant=0):
    """
    Generate a boundary conditions file (.bc) with a constant waterlevel.
    This can be used to enforce a know offset from zero, for instance to
    account for sea level rise.
    """
    quantity = "waterlevelbnd"
    
    # read polyfile as geodataframe
    polyfile_obj = hcdfm.PolyFile(file_pli)
    plinames_list = [x.metadata.name for x in polyfile_obj.objects]
    
    # generate constant forcingmodel object
    ForcingModel_object = hcdfm.ForcingModel()
    for pliname in plinames_list:
        locationname = f"{pliname}_0001"
        qup = [hcdfm.QuantityUnitPair(quantity=quantity, unit="m")]
        Constant_object = hcdfm.Constant(
            name=locationname,
            quantityunitpair=qup,
            datablock=[[constant]], 
            )
        ForcingModel_object.forcing.append(Constant_object)
    
    dir_output = os.path.dirname(file_pli)
    file_pli_name = polyfile_obj.filepath.stem
    file_bc_out = os.path.join(dir_output,f'{quantity}_constant_{file_pli_name}.bc')
    ForcingModel_object.save(filepath=file_bc_out)
    
    # generate hydrolib-core Boundary object to be appended to the ext file
    boundary_object = hcdfm.Boundary(quantity=quantity,
                                      locationfile=file_pli,
                                      forcingfile=ForcingModel_object)
    
    # add the boundary object to the ext file for each polyline in the polyfile
    ext_add_boundary_object_per_polyline(ext_new=ext_new, boundary_object=boundary_object)
    return ext_new


def cmems_nc_to_bc(ext_new, list_quantities, tstart, tstop, file_pli, dir_pattern, dir_output, conversion_dict=None, refdate_str=None):
    if conversion_dict is None:
        conversion_dict = dfmt.get_conversion_dict()
    
    for quantity in list_quantities: # loop over salinitybnd/uxuyadvectionvelocitybnd/etc
        print(f'processing quantity: {quantity}')
        
        quantity_list = get_quantity_list(quantity=quantity)
        
        for quantity_key in quantity_list: # loop over ux/uy
            ncvarname = get_ncvarname(quantity=quantity_key, conversion_dict=conversion_dict)
            dir_pattern_one = str(dir_pattern).format(ncvarname=ncvarname)
            #open regulargridDataset and do some basic stuff (time selection, renaming depth/lat/lon/varname, converting units, etc)
            data_xr_onevar = open_prepare_dataset(dir_pattern=dir_pattern_one, 
                                                  quantity=quantity_key,
                                                  tstart=tstart, tstop=tstop,
                                                  conversion_dict=conversion_dict,
                                                  refdate_str=refdate_str)
            if quantity_key == quantity_list[0]:
                data_xr_vars = data_xr_onevar
            else: # only relevant in case of ux/uy, others all have only one quantity
                data_xr_vars[quantity_key] = data_xr_onevar[quantity_key]
        
        # interpolate regulargridDataset to plipointsDataset
        polyfile_obj = hcdfm.PolyFile(file_pli)
        gdf_points = dfmt.PolyFile_to_geodataframe_points(polyfile_object=polyfile_obj)
        data_interp = dfmt.interp_regularnc_to_plipointsDataset(data_xr_reg=data_xr_vars, gdf_points=gdf_points, load=True)
        
        #convert plipointsDataset to hydrolib ForcingModel
        ForcingModel_object = dfmt.plipointsDataset_to_ForcingModel(plipointsDataset=data_interp)
        
        # generate boundary object for the ext file (quantity, pli-filename, bc-filename)
        file_pli_name = polyfile_obj.filepath.stem
        file_bc_out = os.path.join(dir_output,f'{quantity}_CMEMS_{file_pli_name}.bc')
        ForcingModel_object.save(filepath=file_bc_out)
        boundary_object = hcdfm.Boundary(quantity=quantity,
                                         # locationfile is updated if multiple polylines in polyfile
                                         locationfile=file_pli, 
                                         forcingfile=ForcingModel_object)
        
        # add the boundary object to the ext file for each polyline in the polyfile
        ext_add_boundary_object_per_polyline(ext_new=ext_new, boundary_object=boundary_object)

    return ext_new


def preprocess_ini_cmems_to_nc(**kwargs):
    raise DeprecationWarning("`dfmt.preprocess_ini_cmems_to_nc()` was "
                             "deprecated, use `cmems_nc_to_ini()` instead")


def cmems_nc_to_ini(
        ext_old,
        dir_output,
        list_quantities,
        tstart,
        dir_pattern,
        conversion_dict=None,
        ):
    """
    This makes quite specific 3D initial input files based on CMEMS data.
    delft3dfm is quite picky, so it works with CMEMS files because they have a
    depth variable with standard_name='depth'. If this is not present, or
    dimensions are ordered differently, there will be an error or incorrect
    model results.
    """
    
    if conversion_dict is None:
        conversion_dict = dfmt.get_conversion_dict()
    
    tstart_pd = pd.Timestamp(tstart)
    tstart_str = tstart_pd.strftime("%Y-%m-%d_%H-%M-%S")
    # tstop_pd is slightly higher than tstart_pd to ensure >1 timesteps in all
    # cases
    tstop_pd = tstart_pd + pd.Timedelta(hours=1)
    
    for quan_bnd in list_quantities:
        if (quan_bnd != "salinitybnd") and not quan_bnd.startswith("tracer"):
            # quantities other than salinitybnd, temperaturebnd and tracer* are
            # silently skipped since they are also not supported by delft3dfm
            # temperaturebnd is handled with salinity
            # uxuyadvectionvelocitybnd is not supported
            print(f"quantity {quan_bnd} is not supported by cmems_nc_to_ini, "
                   "silently skipped")
            continue
        
        ncvarname = get_ncvarname(
            quantity=quan_bnd,
            conversion_dict=conversion_dict,
            )
        dir_pattern_one = dir_pattern.format(ncvarname=ncvarname)
        
        if quan_bnd=="salinitybnd":
            # this also handles temperaturebnd
            # 3D initialsalinity/initialtemperature fields are silently ignored
            # initial 3D conditions are only possible via nudging 1st timestep
            # via quantity=nudge_salinity_temperature
            data_xr = xr.open_mfdataset(dir_pattern_one)
            ncvarname_tem = get_ncvarname(
                quantity="temperaturebnd",
                conversion_dict=conversion_dict,
                )
            dir_pattern_tem = dir_pattern.format(ncvarname=ncvarname_tem)
            data_xr_tem = xr.open_mfdataset(dir_pattern_tem)
            data_xr["thetao"] = data_xr_tem["thetao"]
            quantity = "nudge_salinity_temperature"
            varname = None
        elif quan_bnd.startswith("tracer"):
            data_xr = xr.open_mfdataset(dir_pattern_one)
            data_xr = ds_apply_conversion_dict(
                data_xr=data_xr,
                conversion_dict=conversion_dict,
                quantity=quan_bnd,
                )
            quantity = f'initial{quan_bnd.replace("bnd","")}'
            varname = quantity
            data_xr = data_xr.rename_vars({quan_bnd:quantity})
        
        # subset two times. interp to tstart would be the proper way to do it, 
        # but FM needs two timesteps for nudge_salinity_temperature and initial
        # waq vars
        data_xr = _ds_sel_time_outside(
            ds=data_xr, tstart=tstart_pd, tstop=tstop_pd,
            )
        
        # assert that there are at least two timesteps in the resulting dataset
        # delft3dfm will crash if there is only one timestep
        assert len(data_xr.time) >= 2
    
        # interpolate_na for all data_vars, first over lat/lon, then over depth
        for var in data_xr.data_vars:
            data_xr[var] = interpolate_na_multidim(data_xr[var], ["latitude", 
                                                                  "longitude"])
            data_xr[var] = interpolate_na_multidim(data_xr[var], ["depth"])
        
        print('writing file')
        file_output = os.path.join(dir_output,f"{quantity}_{tstart_str}.nc")
        data_xr.to_netcdf(file_output)
        
        #append forcings to ext
        forcing_saltem = hcdfm.ExtOldForcing(
            quantity=quantity,
            varname=varname,
            filename=file_output,
            filetype=hcdfm.ExtOldFileType.NetCDFGridData,
            method=hcdfm.ExtOldMethod.InterpolateTimeAndSpaceSaveWeights, #3
            operand=hcdfm.Operand.override, #O
            )
        ext_old.forcing.append(forcing_saltem)
    
    return ext_old


def preprocess_merge_meteofiles_era5(ext_old, varkey_list, dir_data, dir_output, time_slice):

    if isinstance(varkey_list[0], list):
        varkey_lists = varkey_list
    else:
        varkey_lists = [varkey_list]
    
    for varkey_list in varkey_lists:
        fn_match_pattern = f'era5_.*({"|".join(varkey_list)})_.*.nc' #simpler but selects more files: 'era5_*.nc'
        file_out_prefix = f'era5_{"_".join(varkey_list)}_'
        preprocess = dfmt.preprocess_ERA5 #reduce expver dimension if present
        
        if not os.path.exists(dir_output):
            os.makedirs(dir_output)
        
        file_nc = os.path.join(dir_data,fn_match_pattern)
        
        data_xr_tsel = dfmt.merge_meteofiles(file_nc=file_nc,
                                             time_slice=time_slice, 
                                             preprocess=preprocess,
                                             )
        
        #write to netcdf file
        print('>> writing file (can take a while): ',end='')
        dtstart = dt.datetime.now()
        times_pd = data_xr_tsel['time'].to_series()
        time_start_str = times_pd.iloc[0].strftime("%Y%m%d")
        time_stop_str = times_pd.iloc[-1].strftime("%Y%m%d")
        file_out = os.path.join(dir_output, f'{file_out_prefix}{time_start_str}to{time_stop_str}_ERA5.nc')
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
        else:
            # TODO: add support for other quantities: https://github.com/Deltares/dfm_tools/issues/887
            raise KeyError(f"'{varkey_list}' is not supported by dfmt.preprocess_merge_meteofiles_era5(), "
                           "the files were merged but the ExtModel Forcing cannot be appended.")
        
    return ext_old


def create_model_exec_files(file_mdu, nproc=1, dimrset_folder=None, path_style=None):
    """
    creates a dimr_config.xml and if desired a batfile to run the model
    """
    
    if not os.path.isfile(file_mdu):
        raise FileNotFoundError(f"file_mdu not found: {file_mdu}")
        
    dirname = os.path.dirname(file_mdu)
    mdu_name = os.path.basename(file_mdu)
    file_dimr = os.path.join(dirname,'dimr_config.xml')
    dimr_name = os.path.basename(file_dimr)
    
    # generate dimr_config.xml
    fm_modelname = "DFlowFM"
    control_comp = Start(name=fm_modelname)
    fm_comp = FMComponent(name=fm_modelname, workingDir='.', inputfile=mdu_name,
                          process=nproc, 
                          mpiCommunicator="DFM_COMM_DFMWORLD")
    dimr_model = DIMR(control=control_comp, component=fm_comp)
    print(f"writing {dimr_name}")
    dimr_model.save(file_dimr)
    
    if dimrset_folder is None:
        print("no dimrset_folder provided, cannot write bat/sh file")
        return
    elif dimrset_folder=='docker':
        generate_docker_file(dimr_model=dimr_model)
        return
    
    # continue with dimrset_folder which is not None or 'docker'
    if path_style is None:
        path_style = get_path_style_for_current_operating_system().value
    
    if path_style == 'windows':
        generate_bat_file(dimr_model=dimr_model, dimrset_folder=dimrset_folder)
    else:
        logger.warning(f"path_style/os {path_style} not yet supported by `dfmt.create_model_exec_files()`, no bat/sh file is written")


def generate_bat_file(dimr_model, dimrset_folder=None):
    """
    generate bat file for running on windows
    """
    
    if dimr_model.filepath is None:
        raise Exception('first save the dimr_model before passing it to generate_bat_file')
    
    dirname = os.path.dirname(dimr_model.filepath)
    file_bat = os.path.join(dirname, "run_parallel.bat")
    bat_name = os.path.basename(file_bat)
    
    dimr_name = os.path.basename(dimr_model.filepath)
    mdu_name = os.path.basename(dimr_model.component[0].inputFile)
    nproc = dimr_model.component[0].process
    
    if not os.path.exists(dimrset_folder):
        raise FileNotFoundError(f"dimrset_folder not found: {dimrset_folder}")
    
    bat_str = fr"""rem User input
set dimrset_folder="{dimrset_folder}"
set MDU_file="{mdu_name}"
set partitions={nproc}

rem Partition the network and mdu
call %dimrset_folder%\x64\bin\run_dflowfm.bat "--partition:ndomains=%partitions%:icgsolver=6" %MDU_file%

rem Execute the simulation
call %dimrset_folder%\x64\bin\run_dimr_parallel.bat %partitions% {dimr_name}

rem To prevent the DOS box from disappearing immediately: enable pause on the following line
pause
"""
    print(f"writing {bat_name}")
    with open(file_bat,'w') as f:
        f.write(bat_str)


def generate_docker_file(dimr_model):
    """
    generate run_docker.sh file for running on windows or unix with docker
    """
    
    if dimr_model.filepath is None:
        raise Exception('first save the dimr_model before passing it to generate_bat_file')
    
    dirname = os.path.dirname(dimr_model.filepath)
    file_docker = os.path.join(dirname, "run_docker.sh")
    docker_name = os.path.basename(file_docker)
    
    dimr_name = os.path.basename(dimr_model.filepath)
    mdu_name = os.path.basename(dimr_model.component[0].inputFile)
    nproc = dimr_model.component[0].process
    
    docker_str = fr"""#!/bin/bash
# export OMP_NUM_THREADS=1 # not sure what for
export I_MPI_FABRICS=shm # required on windows

# first pull or load a docker container
# docker pull deltares/delft3dfm
# docker load -i <file.tar>
# RUN THIS run_docker.sh FILE ON COMMAND LINE WITH (shm-size and ulimit seem optional):
# docker run -v /path/to/dimr:/data -t deltares/delft3dfm:latest /data/run_docker.sh --shm-size=4gb --ulimit stack=-1

# stop after an error occured
set -e

# set number of partitions
nPart={nproc}

# location of the binaries inside Docker image
delft3d=/opt/delft3dfm_latest/lnx64

# DIMR input-file; must already exist!
dimrFile={dimr_name}

# Replace number of processes in DIMR file
PROCESSSTR="$(seq -s " " 0 $((nPart-1)))"
sed -i "s/\(<process.*>\)[^<>]*\(<\/process.*\)/\1$PROCESSSTR\2/" $dimrFile

# Read MDU file from DIMR-file
# mduFile="$(sed -n 's/\r//; s/<inputFile>\(.*\).mdu<\/inputFile>/\1/p' $dimrFile)".mdu
mduFile={mdu_name}

if [ "$nPart" == "1" ]; then
    $delft3d/bin/run_dimr.sh -m $dimrFile
else
    $delft3d/bin/run_dflowfm.sh --partition:ndomains=$nPart:icgsolver=6 $mduFile
    $delft3d/bin/run_dimr.sh -c $nPart -m $dimrFile
fi
"""
    print(f"writing {docker_name}")
    # run_docker.sh requires unix file endings, so we use newline='\n'
    with open(file_docker, 'w', newline='\n') as f:
        f.write(docker_str)


def make_paths_relative(mdu_file:str):
    """
    Making paths on windows and unix relative by removing the dir_model prefix from paths in the mdu and ext files
    This is a temporary workaround until implemented in https://github.com/Deltares/HYDROLIB-core/issues/532    

    Parameters
    ----------
    mdu_file : str
        path to mdu_file.

    Returns
    -------
    None.

    """
    dir_model = os.path.dirname(mdu_file)
    mdu_existing = hcdfm.FMModel(mdu_file, recurse=False)
    file_list = [mdu_file]
    ext_old = mdu_existing.external_forcing.extforcefile
    if ext_old is not None:
        file_list.append(ext_old.filepath)
    ext_new = mdu_existing.external_forcing.extforcefilenew
    if ext_new is not None:
        file_list.append(ext_new.filepath)
    for filename in file_list:
        if filename is None:
            continue
        with open(filename, 'r') as file:
            filedata = file.read()
        filedata = filedata.replace(dir_model.replace('\\','/')+'/', '')
        with open(filename, 'w') as file:
            file.write(filedata)

