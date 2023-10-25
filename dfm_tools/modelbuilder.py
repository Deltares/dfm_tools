"""
Created on Tue Apr  4 16:12:56 2023

@author: groenenb, sclaan
edited by: veenstra, zijlker

This is a proof of concept that will be properly coded in the near future and probably end up under hydromt-delft3dfm
Since the functions in this script contain hardcoded parameters, it is not exposed to public and you need to import like dfmt.modelbuilder.[function]
"""

import os
from pathlib import Path
import xarray as xr
import pandas as pd
import dfm_tools as dfmt
import hydrolib.core.dflowfm as hcdfm
import datetime as dt
import glob
import geopandas as gpd
import numpy as np
import pyproj
import shapely
from shapely.geometry import LineString, MultiPoint, Point
from scipy.spatial import KDTree
from hydrolib.core.dimr.models import DIMR, FMComponent, Start
import warnings

def generate_coastline_obspoints(lon_min : float, # degrees
                                 lon_max : float, # degrees
                                 lat_min : float, # degrees
                                 lat_max : float, # degrees
                                 interval : float,# degrees
                                 resolution : str ,# degrees
                                 threshold_mindepth : float, # meters
                                 file_nc : Path, # path to netcdf file
                                 output_file : Path): # path to output file,):

    #%% create obs file and snap to grid 
    # regular  grid in coast
    lat =  np.linspace(lat_min,lat_max,int((lat_max-lat_min)/0.5))
    lon =  np.linspace(lon_min,lon_max,int((lon_max-lon_min)/0.5))

    # Points along coastline with dfm_tools
    coastline = dfmt.get_coastlines_gdb(res = resolution,
                             bbox = [lon_min, lat_min, lon_max, lat_max],
                             crs = "EPSG:4326")
    linestrings = [polygon.boundary for polygon in coastline['geometry']]
    linstrings_cropped = [shapely.ops.clip_by_rect(linestring, lon_min, lat_min, lon_max, lat_max) for linestring in linestrings]
    shp_clipped = gpd.GeoDataFrame(geometry = gpd.GeoSeries(linstrings_cropped))
    shp_clipped = shp_clipped[~shp_clipped.is_empty]
    
    #Build GeoSeries with points along coastline 
    cpoints = gpd.GeoSeries()
    for index, line in shp_clipped.iterrows():
        if 'MULTILINESTRING' in str(line.values):       
            shapes = str(line.values[0]).strip().split("\n")       
            gdf = gpd.GeoDataFrame({'geometry': shapes})
            gdf['geometry'] = gpd.GeoSeries.from_wkt(gdf['geometry'])
            gdf = gdf.set_geometry('geometry').explode(index_parts=True)
            for iindex, iline in  gdf.iterrows():
                distances = np.arange(0, iline.geometry.length, interval)
                points = MultiPoint([LineString(iline.geometry).interpolate(distance) for distance in distances])
                gs =gpd.GeoSeries(Point(pnt.x,pnt.y) for pnt in points.geoms)
                cpoints=pd.concat([cpoints, gs])   #.append(gs)
        else:
            distances = np.arange(0, line.geometry.length, interval)
            points = MultiPoint([LineString(line.geometry).interpolate(distance) for distance in distances])
            gs =gpd.GeoSeries(Point(pnt.x,pnt.y) for pnt in points.geoms)
            cpoints=pd.concat([cpoints, gs])# cpoints=cpoints.append(gs)
    locs_coast = gpd.GeoDataFrame(cpoints, geometry=cpoints.geometry, crs="EPSG:4326")

    # Snap to model points and -5 m depth ## TO DO
    obs=pd.DataFrame()
    obs['x'] = locs_coast.geometry.x.values
    obs['y'] = locs_coast.geometry.y.values

    #%% calculate face center x, y and z
    print('interpolating cell centers from nodes')
    ds_net = xr.open_dataset(file_nc)
    net_face_x = []
    net_face_y = []
    net_face_z = []
    for face in ds_net.mesh2d_nFaces:
        idx_node = ds_net.mesh2d_face_nodes\
            .sel(mesh2d_nFaces = face)\
            .dropna(dim='mesh2d_nMax_face_nodes').data - 1
        idx_node =  [int(idx_node[i]) for i in range(len(idx_node))]
        net_face_x_sel = ds_net.mesh2d_node_x.sel(mesh2d_nNodes = idx_node).data.mean()
        net_face_y_sel = ds_net.mesh2d_node_y.sel(mesh2d_nNodes = idx_node).data.mean()
        net_face_z_sel = ds_net.mesh2d_node_z.sel(mesh2d_nNodes = idx_node).data.mean()
        net_face_x.append(net_face_x_sel)
        net_face_y.append(net_face_y_sel)
        net_face_z.append(net_face_z_sel)

    #Make dictionary with cell centers
    net_faces = pd.DataFrame({'x':net_face_x, 'y':net_face_y, 'z':net_face_z})

    #Select only cells with z < threshold_mindepth
    bool_valid_cells = (net_faces['z']<-threshold_mindepth)

    #creating kdtree with valid cell centers (cartesian coordinates)
    def xlonylat2xyzcartesian(data):
        """
        necessary to calculate cartesian distances, otherwise nearest neigbour can fail.
        https://stackoverflow.com/questions/45127141/find-the-nearest-point-in-distance-for-all-the-points-in-the-dataset-python
        """
        R = 6367
        phi = np.deg2rad(data['y'])
        theta = np.deg2rad(data['x'])
        data = pd.DataFrame()
        data['x_cart'] = R * np.cos(phi) * np.cos(theta)
        data['y_cart'] = R * np.cos(phi) * np.sin(theta)
        data['z_cart'] = R * np.sin(phi)
        return data

    data_celcenxy_valid = net_faces.loc[bool_valid_cells, ['x', 'y']].reset_index() #pd.DataFrame({'x':net_node_x[bool_valid_cells],'y':net_node_y[bool_valid_cells]})#,'area':data_cellarea[bool_valid_cells]})
    data_celcenxy_valid_cart = xlonylat2xyzcartesian(data_celcenxy_valid)
    tree = KDTree(data_celcenxy_valid_cart[['x_cart','y_cart','z_cart']])

    def dist_to_arclength(chord_length):
        """
        https://stackoverflow.com/questions/45127141/find-the-nearest-point-in-distance-for-all-the-points-in-the-dataset-python
        """
        R = 6367 # earth radius
        central_angle = 2*np.arcsin(chord_length/(2.0*R)) 
        arclength = R*central_angle*1000
        return arclength

    #finding nearest cellcenter-neighbors of each obspoint in file
    data_obsorg_cart = xlonylat2xyzcartesian(obs)
    distance_cart, index = tree.query(data_obsorg_cart, k=1)
    data_celcenxy_validsel = data_celcenxy_valid.loc[index,:]

    # write 
    obs_snapped=pd.DataFrame()
    obs_snapped['x'] = data_celcenxy_validsel['x'].values
    obs_snapped['y'] = data_celcenxy_validsel['y'].values
    obs_snapped = obs_snapped.drop_duplicates()

    #Create obs file for dflowfm
    obs_sfincs = hcdfm.XYNModel()
    points = [hcdfm.XYNPoint(x=x, y=y, n=f'sfincs_bnd_{i:03d}') for i, (x, y) in enumerate(zip(obs_snapped['x'].values, obs_snapped['y'].values))]
    obs_sfincs.points = points
    return obs_sfincs

def cmems_nc_to_bc(ext_bnd, list_quantities, tstart, tstop, file_pli, dir_pattern, dir_output, refdate_str):
    #input examples in https://github.com/Deltares/dfm_tools/blob/main/tests/examples/preprocess_interpolate_nc_to_bc.py
    
    file_bc_basename = os.path.basename(file_pli).replace('.pli','')
    for quantity in list_quantities:
        print(f'processing quantity: {quantity}')
        
        tstart, tstop = dfmt.round_timestamp_to_outer_noon(tstart,tstop)
        #open regulargridDataset and do some basic stuff (time selection, renaming depth/lat/lon/varname, converting units, etc)
        data_xr_vars = dfmt.open_dataset_extra(dir_pattern=dir_pattern, quantity=quantity,
                                               tstart=tstart, tstop=tstop,
                                               conversion_dict=dfmt.get_conversion_dict(),
                                               refdate_str=refdate_str)
        #interpolate regulargridDataset to plipointsDataset
        data_interp = dfmt.interp_regularnc_to_plipoints(data_xr_reg=data_xr_vars, file_pli=file_pli)
        
        #convert plipointsDataset to hydrolib ForcingModel
        ForcingModel_object = dfmt.plipointsDataset_to_ForcingModel(plipointsDataset=data_interp)
        
        file_bc_out = os.path.join(dir_output,f'{quantity}_{file_bc_basename}_CMEMS.bc')
        
        ForcingModel_object.save(filepath=file_bc_out)
        
        #generate boundary object for the ext file (quantity, pli-filename, bc-filename)
        boundary_object = hcdfm.Boundary(quantity=quantity,
                                         locationfile=file_pli,
                                         forcingfile=ForcingModel_object)
        ext_bnd.boundary.append(boundary_object)
    
    return ext_bnd

    
def preprocess_ini_cmems_to_nc(ext_old, tstart, dir_data, dir_out):
    
    file_nc_list_so = glob.glob(f'{dir_data}\\cmems_so_*.nc')
    file_nc_list_thetao = glob.glob(f'{dir_data}\\cmems_thetao_*.nc')
    file_nc_list = file_nc_list_so + file_nc_list_thetao
    
    print(f'opening {len(file_nc_list)} datasets')
    data_xr = xr.open_mfdataset(file_nc_list)
    
    # fill nans
    # start with lat/lon to avoid values from shallow coastal areas in deep layers
    # first interpolate nans to get smooth filling of e.g. islands, this cannot fill nans at the edge of the dataset
    data_xr = data_xr.interpolate_na(dim='latitude').interpolate_na(dim='longitude')
    
    # then use bfill/ffill to fill nans at the edge for lat/lon/depth
    data_xr = data_xr.ffill(dim='latitude').bfill(dim='latitude')
    data_xr = data_xr.ffill(dim='longitude').bfill(dim='longitude')
    data_xr = data_xr.ffill(dim='depth').bfill(dim='depth')
    
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
        
        data_xr_tsel = dfmt.merge_meteofiles(file_nc=file_nc, time_slice=time_slice, 
                                             preprocess=preprocess,
                                             add_global_overlap=False, #GTSM specific: extend data beyond -180 to 180 longitude
                                             zerostart=False) #GTSM specific: extend data with 0-value fields 1 and 2 days before all_tstart
        
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
        
    return ext_old


def create_model_exec_files(file_dimr, file_mdu, model_name, nproc=1, dimrset_folder=None, path_style=None):
    """
    creates a dimr_config.xml and if desired a batfile to run the model
    """
    
    mdu_name = os.path.basename(file_mdu)
    
    # generate dimr_config.xml
    control_comp = Start(name=model_name)
    fm_comp = FMComponent(name=model_name, workingDir='.', inputfile=mdu_name,
                          process=nproc, 
                          mpiCommunicator="DFM_COMM_DFMWORLD")
    dimr_model = DIMR(control=control_comp, component=fm_comp)
    print(f"writing {file_dimr}")
    dimr_model.save(file_dimr)
    
    # TODO: hydrolib-core does not support multiple cores properly: https://github.com/Deltares/dfm_tools/issues/214
    # therefore we manually replace it in the file
    print(f"re-writing {file_dimr}")
    with open(file_dimr,'r') as f:
        lines = f.readlines()
    str_from = f"<process>{nproc}</process>"
    nproc_range_str = " ".join([str(x) for x in range(nproc)])
    str_to = f"<process>{nproc_range_str}</process>"
    lines_new = [line.replace(str_from,str_to) for line in lines]
    with open(file_dimr,'w') as f:
        for line in lines_new:
            f.write(line)
    
    if path_style is None:
        from hydrolib.core.utils import get_path_style_for_current_operating_system
        path_style = get_path_style_for_current_operating_system().value
    
    #TODO: currently only bat files are supported (for windows), but linux extension can easily be made
    if path_style == 'windows':
        _generate_bat_file(dimr_model=dimr_model, dimrset_folder=dimrset_folder)
    else:
        warnings.warn(UserWarning(f"path_style/os {path_style} not yet supported by `dfmt.create_model_exec_files()`, no bat/sh file is written"))


def _generate_bat_file(dimr_model, dimrset_folder=None):
    """
    generate bat file for running on windows
    """
    
    if dimr_model.filepath is None:
        raise Exception('first save the dimr_model before passing it to generate_bat_file')
    
    dirname = os.path.dirname(dimr_model.filepath)
    file_bat = os.path.join(dirname, "run_parallel.bat")
    
    dimr_name = os.path.basename(dimr_model.filepath)
    mdu_name = os.path.basename(dimr_model.component[0].inputFile)
    nproc = dimr_model.component[0].process
    if dimrset_folder is None:
        dimrset_folder = r"c:\Program Files\Deltares\Delft3D FM Suite 2023.02 HMWQ\plugins\DeltaShell.Dimr\kernels"
    
    if not os.path.exists(dimrset_folder):
        raise FileNotFoundError(f"dimrset_folder not found: {dimrset_folder}")
    
    bat_str = fr"""
rem User input
set dimrset_folder="{dimrset_folder}"
set MDU_file="{mdu_name}"
set partitions={nproc}

rem Partition the network and mdu
call %dimrset_folder%\x64\dflowfm\scripts\run_dflowfm.bat "--partition:ndomains=%partitions%:icgsolver=6" %MDU_file%

rem Execute the simulation
call %dimrset_folder%\x64\dimr\scripts\run_dimr_parallel.bat %partitions% {dimr_name}

rem To prevent the DOS box from disappearing immediately: enable pause on the following line
pause
"""
    print(f"writing {file_bat}")
    with open(file_bat,'w') as f:
        f.write(bat_str)

