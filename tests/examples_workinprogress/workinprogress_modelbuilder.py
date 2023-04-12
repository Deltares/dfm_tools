# -*- coding: utf-8 -*-
"""
Created on Tue Apr  4 16:12:56 2023

@author: groenenb, laan_st, veenstra
"""

import os
import matplotlib.pyplot as plt
plt.close('all')
import dfm_tools as dfmt
from dfm_tools import modelbuilder as mb #different import for modelbuilder since it is not exposed publicly
import hydrolib.core.dflowfm as hcdfm
import xarray as xr
import pandas as pd

## input
model_name = 'Bonaire'
dir_output = r'p:\11209231-003-bes-modellering\hydrodynamica\hackathon\preprocessing\ModelBuilderOutput_JV'
path_style = 'unix' # windows / unix
overwrite = False # used for downloading of forcing data. Always set to True when changing the domain

#TODO: files that are not created in this script: obsfiles, dimr.xml and the submit script
#TODO: reference run in: p:\11209231-003-bes-modellering\hydrodynamica\hackathon\simulations\run001_mapInterval_1800\
#TODO: also compare settings to p:\11208054-004-dcsm-fm\models\3D_DCSM-FM\2013-2017\B05_hydrolib_JV\DCSM-FM_0_5nm.mdu (e.g. tlfSmo)
# domain
lon_min, lon_max, lat_min, lat_max = -68.55, -67.9, 11.8, 12.6
bnd_dlon_dlat = 0.06

#dates as understood by pandas.period_range(). ERA5 has freq='M' (month) and CMEMS has freq='D' (day)
date_min = '2022-11-01'
date_max = '2022-11-03'
ref_date = '2022-01-01'

# meteo (ERA5)
# option 1 (2D-model) = [airpressure, windx, windy, charnock]
# option 2 (3D-model) = [airpressure, windx, windy, charnock, dewpoint, airtemperature, cloudiness, solarradiation, longwaveradiation, rainfall, evaporation]
ERA5_meteo_option = 2

## process
# make dirs
if not os.path.isdir(dir_output):
    os.mkdir(dir_output)
dir_output_data = os.path.join(dir_output, 'data')
if not os.path.isdir(dir_output_data):
    os.mkdir(dir_output_data)

#%% initialize mdu file
mdu_file = os.path.join(dir_output, f'{model_name}.mdu')
mdu = hcdfm.FMModel()


#%%bnd generation
pli_polyfile = mb.generate_bndpli(lon_min, lon_max, lat_min, lat_max, dlon=bnd_dlon_dlat, dlat=bnd_dlon_dlat, name=f'{model_name}_bnd')
poly_file = os.path.join(dir_output, f'{model_name}.pli')
pli_polyfile.save(poly_file)


#%% grid generation
# GEBCO bathymetry and Grid generation
#select and plot bathy
file_nc_bathy = r'p:\metocean-data\open\GEBCO\2021\GEBCO_2021.nc'
data_bathy = xr.open_dataset(file_nc_bathy)
data_bathy_sel = data_bathy.sel(lon=slice(lon_min-1,lon_max+1),lat=slice(lat_min-1,lat_max+1))

#TODO: grid generation and bathy-refinement is still to be improved in meshkernel (https://github.com/Deltares/dfm_tools/issues/234)
mk_object = mb.make_basegrid(lon_min, lon_max, lat_min, lat_max) #TODO: should be sperical, but is cartesian >> is_geographic keywork does not work yet

#refine
min_face_size = 200/(40075*1000/360) #convert meters to degrees
mb.refine_basegrid(mk=mk_object, data_bathy_sel=data_bathy_sel, min_face_size=min_face_size) #TODO: min_face_size is now in degrees instead of meters (maybe already works when is_geographic=True?)

#cutcells
file_ldb = r'p:\11209231-003-bes-modellering\hydrodynamica\hackathon\preprocessing\grid\coastline.pli'
dfmt.meshkernel_delete_withpol(mk=mk_object,file_ldb=file_ldb)
#TODO: illegalcells.pol necessary?

#TODO: cleanup grid necessary?
# print('mk_object.mesh2d_get_obtuse_triangles_mass_centers()')
# print(mk_object.mesh2d_get_obtuse_triangles_mass_centers())
# print('mk_object.mesh2d_get_hanging_edges()')
# print(mk_object.mesh2d_get_hanging_edges())

#convert to xugrid
xu_grid_uds = dfmt.meshkernel_to_UgridDataset(mk=mk_object) #TODO: after oldext is in main, remove face_coordinates from grid?

#interp bathy
data_bathy_interp = data_bathy_sel.interp(lon=xu_grid_uds.obj.mesh2d_node_x, lat=xu_grid_uds.obj.mesh2d_node_y).reset_coords(['lat','lon']) #interpolates lon/lat gebcodata to mesh2d_nNodes dimension #TODO: if these come from xu_grid_uds (without ojb), the mesh2d_node_z var has no ugrid accessor since the dims are lat/lon instead of mesh2d_nNodes
xu_grid_uds['mesh2d_node_z'] = data_bathy_interp.elevation.clip(max=10)

fig, ax = plt.subplots()
xu_grid_uds.grid.plot(ax=ax,linewidth=1)

fig, ax = plt.subplots()
xu_grid_uds.mesh2d_node_z.ugrid.plot(ax=ax,center=False,vmin=-500,vmax=10)
#ctx.add_basemap(ax=ax, crs='EPSG:4326', attribution=False)

#write xugrid grid to netcdf
netfile  = os.path.join(dir_output, f'{model_name}_net.nc')
xu_grid_uds.ugrid.to_netcdf(netfile)
mdu.geometry.netfile = netfile #TODO: path is windows/unix dependent #TODO: providing os.path.basename(netfile) raises "ValidationError: 1 validation error for Geometry - netfile:   File: `C:\SnapVolumesTemp\MountPoints\{45c63495-0000-0000-0000-100000000000}\{79DE0690-9470-4166-B9EE-4548DC416BBD}\SVROOT\DATA\dfm_tools\tests\examples_workinprogress\Bonaire_net.nc` not found, skipped parsing." (wrong current directory)
#breakit

dir_output_data_cmems = os.path.join(dir_output_data, 'cmems')
if not os.path.isdir(dir_output_data_cmems):
    os.mkdir(dir_output_data_cmems)

if 1:
    #%% new ext: initial and open boundary condition
    
    ext_file_new = os.path.join(dir_output, f'{model_name}_new.ext')
    ext_new = hcdfm.ExtModel()
    
    # CMEMS - download
        
    mb.download_meteodata_oceandata(
        longitude_min = lon_min, longitude_max = lon_max, latitude_min = lat_min, latitude_max = lat_max,
        model = 'CMEMS',
        overwrite = overwrite,
        date_min = date_min, date_max = date_max,
        varlist = ['so','thetao','uo','vo','zos'],
        dir_output = dir_output_data_cmems)
    
    # CMEMS - boundary conditions file (.bc) (and add to ext_bnd)
    ext_new = mb.preprocess_interpolate_nc_to_bc(ext_bnd=ext_new,
                                                 refdate_str = 'minutes since '+ref_date+' 00:00:00 +00:00',
                                                 dir_output = dir_output,
                                                 model = 'CMEMS',
                                                 tstart=pd.Timestamp(date_min)-pd.Timedelta(hours=12), tstop=pd.Timestamp(date_max)+pd.Timedelta(hours=12), #TODO: to account for noon-fields of CMEMS, build in safety?
                                                 list_plifiles = [poly_file],
                                                 list_quantities = ['salinitybnd','temperaturebnd','uxuy','waterlevelbnd','tide'],
                                                 #list_quantities = ['waterlevelbnd','salinitybnd','temperaturebnd','uxuy','tide'], #TODO: see comment block below
                                                 dir_sourcefiles_hydro = dir_output_data_cmems)
    #TODO: when adding both waterlevelbnd and tide as waterlevelbnd they should be consequetive. If there are other quantities in between, the model crashes with a update_ghostboundvals error, report this:
    """
    ** INFO   : partition_fixorientation_ghostlist: number of reversed flowlinks=              0
    ** INFO   : Done partitioning model.
    ** ERROR  : update_ghostboundvals: not all ghost boundary flowlinks are being updated
    ** INFO   : Closed file : /p/11209231-003-bes-modellering/hydrodynamica/hackathon/preprocessing/ModelBuilderOutput_JV/Bonaire_old.ext
    ** INFO   : partition_fixorientation_ghostlist: number of reversed flowlinks=              0
    ** INFO   : Done partitioning model.
    ** ERROR  : update_ghostboundvals: not all ghost boundary flowlinks are being updated
    ** INFO   : Closed file : /p/11209231-003-bes-modellering/hydrodynamica/hackathon/preprocessing/ModelBuilderOutput_JV/Bonaire_old.ext
    Abort(1) on node 0 (rank 0 in comm 0): application called MPI_Abort(MPI_COMM_WORLD, 1) - process 0
    Abort(1) on node 1 (rank 1 in comm 0): application called MPI_Abort(MPI_COMM_WORLD, 1) - process 1
    """
    #When doing a sequential run, the model initializes successfully, but first timestep cannot be solved.
    """
    ** WARNING: Comp. time step average below threshold:  0.5742E-03 <  0.1000E+00.
    ** INFO   : Performing direct write of solution state...
    ** INFO   : Simulation current time: nt = 100, time1 =  26265600.06s (2022-11-01T00:00:00Z).
    ** INFO   : Done writing solution state.
    #### ERROR: dimr update ABORT,: myNameDFlowFM update failed, with return value 22 
    """
    #TODO: now also happened once with "correct" ordering, but after ~80% of simulation with auto_bnd and bnd_dlon_dlat=0.08 and 0.06
    """
    ** WARNING: Comp. time step average below threshold:  0.9077E-01 <  0.1000E+00.
    ** INFO   : Performing direct write of solution state...
    ** INFO   : Simulation current time: nt = 10632, time1 =  26413083.42s (2022-11-02T16:58:03Z).
    ** INFO   : Done writing solution state.
    #### ERROR: dimr update ABORT,: myNameDFlowFM update failed, with return value 22 
    """
    #save new ext file
    ext_new.save(filepath=ext_file_new,path_style=path_style)
    
    
    #%% old ext
if 1:    
    # CMEMS - initial condition file
    ext_file_old = os.path.join(dir_output, f'{model_name}_old.ext')
    ext_old = hcdfm.ExtOldModel()
    
    ext_old = mb.preprocess_ini_cmems_to_nc(ext_old=ext_old,
                                            tstart=date_min,
                                            dir_data=dir_output_data_cmems,
                                            dir_out=dir_output)
    
    # ERA5 - download
    dir_output_data_era5 = os.path.join(dir_output_data,'ERA5')
    if not os.path.exists(dir_output_data_era5):
        os.mkdir(dir_output_data_era5)
        
    if ERA5_meteo_option == 1: #TODO: pass option instead of varlist to fuctions?
        varlist = [['msl','u10n','v10n','chnk']]
    elif ERA5_meteo_option == 2:
        varlist = [['msl','u10n','v10n','chnk'],['d2m','t2m','tcc'],['ssr','strd'],['mer','mtpr']]
    
    mb.download_meteodata_oceandata(
        longitude_min = lon_min, longitude_max = lon_max, latitude_min = lat_min, latitude_max = lat_max,
        model = 'ERA5',
        overwrite = overwrite,
        date_min = date_min, date_max = date_max,
        varlist = varlist, # check variables_dict in dfmt.download_ERA5() for valid names
        dir_output = dir_output_data_era5)
    
    # ERA5 meteo - convert to netCDF for usage in Delft3D FM
    ext_old = mb.preprocess_merge_meteofiles(ext_old=ext_old,
            mode = 'ERA5',
            varkey_list = varlist,
            dir_data = dir_output_data_era5,
            dir_output = dir_output,
            time_slice = slice(date_min, date_max))
    
    ext_old.save(filepath=ext_file_old,path_style=path_style)


#%% .mdu settings
mdu.geometry.bedlevuni = 5
mdu.geometry.kmx = 20
mdu.geometry.layertype = 1
mdu.geometry.numtopsig = 20
mdu.geometry.sigmagrowthfactor = 1.2
mdu.geometry.dxdoubleat1dendnodes = 1
mdu.geometry.changevelocityatstructures = 0
mdu.geometry.changestructuredimensions = 1
mdu.geometry.numtopsiguniform = 1
mdu.geometry.dztop = 5.0
mdu.geometry.floorlevtoplay = -5.0
mdu.geometry.dztopuniabovez = -100.0
mdu.geometry.keepzlayeringatbed = 2

mdu.numerics.tlfsmo = 86400
mdu.numerics.izbndpos = 1
mdu.numerics.mintimestepbreak = 0.1
mdu.numerics.keepstbndonoutflow = 1

mdu.physics.tidalforcing = 1
mdu.physics.salinity = 1
mdu.physics.temperature = 5
mdu.physics.initialsalinity = 33.8
mdu.physics.temperature = 5
mdu.physics.initialtemperature = 29.3
mdu.physics.rhomean = 1023
mdu.physics.secchidepth = 4
mdu.physics.salimax = 50
mdu.physics.tempmax = 50
#mdu.physics.iniwithnudge = 2 #TODO: commented in oldextfile in reference run, initial sal/tem profiles from deep layer were used instead (not yet derived, but 3D inifields also do not have an effect)

mdu.wind.icdtyp = 4
#mdu.wind.cdbreakpoints = [0.025] #TODO: overwritten by spacevarying charnock from ERA5
mdu.wind.rhoair = 1.2265
mdu.wind.relativewind = 0.5
mdu.wind.pavbnd = 101330

mdu.external_forcing.extforcefile = ext_file_old
#mdu.external_forcing.extforcefilenew = ext_file_new #TODO: extfile not found with path_style='unix': https://github.com/Deltares/HYDROLIB-core/issues/516

mdu.time.refdate = pd.Timestamp(ref_date).strftime('%Y%m%d')
mdu.time.tunit = 'S'
mdu.time.dtmax = 30
mdu.time.startdatetime = pd.Timestamp(date_min).strftime('%Y%m%d%H%M%S')
mdu.time.stopdatetime = pd.Timestamp(date_max).strftime('%Y%m%d%H%M%S')
mdu.time.autotimestep = 3

mdu.output.obsfile = [os.path.join(dir_output,x) for x in ['osm_beach_centroids_offset_bonaire.xyn','stations_obs.xyn']]
mdu.output.hisinterval = [60]
mdu.output.mapinterval = [1800]#[86400]
mdu.output.rstinterval = [0] #TODO: default is 0.0, but this translates to [0.0 tstart tstop] instead of [0]: https://github.com/Deltares/HYDROLIB-core/issues/525#issuecomment-1505111924
mdu.output.statsinterval = [3600]
#TODO: disable many outputs (preferrably change many default values): https://github.com/Deltares/HYDROLIB-core/issues/525#issuecomment-1505111924


#%% export model
mdu.save(mdu_file,path_style=path_style)

mdu_contents = open(str(mdu_file),'r').readlines()
for iL, line in enumerate(mdu_contents):
    if line.startswith('netFile'): #TODO: temporary workaround for non-meshkernel grid (converted to spherical with interacter)
        mdu_contents[iL] = 'netFile = bonaire_spherical_net.nc\n'
    if line.startswith('[External Forcing]'): #TODO: temporary workaround for extforcefilenew
        mdu_contents[iL] = line+'extForceFilenew = Bonaire_new.ext'+'\n'
    pass
with open(mdu_file,'w') as f:
    f.write(''.join(mdu_contents))

#TODO: if windows/unix newextfile validation is fixed, use relative paths in .ext and in .mdu files



#%% model run
#TODO: hydrolib-core written mdu-file contains unused keywords
"""
** WARNING: While reading 'Bonaire.mdu': keyword [numerics] qhrelax=0.01 was in file, but not used. Check possible typo.
** WARNING: While reading 'Bonaire.mdu': keyword [wind] windspeedbreakpoints=0.0 100.0 was in file, but not used. Check possible typo.
** WARNING: While reading 'Bonaire.mdu': keyword [output] wrishp_enc=0 was in file, but not used. Check possible typo.
** WARNING: While reading 'Bonaire.mdu': keyword [output] waterlevelclasses=0.0 was in file, but not used. Check possible typo.
** WARNING: While reading 'Bonaire.mdu': keyword [output] waterdepthclasses=0.0 was in file, but not used. Check possible typo.
"""

#TODO: salinity instable, also waterlevel and velocity magnitude are instable at northeast side of island
"""
** INFO   :  Min. salinity limited, number of cells Limmin =           20
** INFO   :  Min. salinity limited, min =  -1.037033177733807E-005
"""
