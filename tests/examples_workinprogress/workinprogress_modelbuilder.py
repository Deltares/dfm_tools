# -*- coding: utf-8 -*-
"""
Created on Tue Apr  4 16:12:56 2023

@author: groenenb, laan_st, veenstra
"""

import datetime as dt
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

# path to polyline: this is the only file that is not created in the script, besides: obsfiles, dimr.xml and the submit script
poly_file = os.path.join(dir_output, 'bnd.pli')
#TODO: automate plifile generation with lat/lon bnds and dxy (possible via meshkernelpy?)

#TODO: reference run in: p:\11209231-003-bes-modellering\hydrodynamica\hackathon\simulations\run001_mapInterval_1800_waq_newGrid

# domain
lon_min, lon_max, lat_min, lat_max = -68.55, -67.9, 11.8, 12.6

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

#%% initialize mdu/ext files
mdu_file = os.path.join(dir_output, f'{model_name}.mdu')
mdu = hcdfm.FMModel()
ext_file_new = os.path.join(dir_output, f'{model_name}_new.ext')
ext_new = hcdfm.ExtModel()
ext_file_old = os.path.join(dir_output, f'{model_name}_old.ext')
ext_old = hcdfm.ExtOldModel()


#%% grid generation
# GEBCO bathymetry and Grid generation
#select and plot bathy
file_nc_bathy = r'p:\metocean-data\open\GEBCO\2021\GEBCO_2021.nc'
data_bathy = xr.open_dataset(file_nc_bathy)
data_bathy_sel = data_bathy.sel(lon=slice(lon_min-1,lon_max+1),lat=slice(lat_min-1,lat_max+1))

#TODO: grid generation/refinement based on bathy still to be improved in meshkernel (https://github.com/Deltares/dfm_tools/issues/234)
mk_object = mb.make_basegrid(lon_min, lon_max, lat_min, lat_max) #TODO: should be sperical, but is cartesian >> is_geographic keywork does not work yet

#refine
min_face_size = 200/(40075*1000/360) #convert meters to degrees
mb.refine_basegrid(mk=mk_object, data_bathy_sel=data_bathy_sel, min_face_size=min_face_size) #TODO: min_face_size is now in degrees instead of meters (maybe already works when is_geographic=True?)

#cutcells
file_ldb = r'p:\11209231-003-bes-modellering\hydrodynamica\hackathon\preprocessing\grid\coastline.pli'
dfmt.meshkernel_delete_withpol(mk=mk_object,file_ldb=file_ldb)

#TODO: cleanup grid necessary?
# print('mk_object.mesh2d_get_obtuse_triangles_mass_centers()')
# print(mk_object.mesh2d_get_obtuse_triangles_mass_centers())
# print('mk_object.mesh2d_get_hanging_edges()')
# print(mk_object.mesh2d_get_hanging_edges())

#convert to xugrid
xu_grid_uds = dfmt.meshkernel_to_UgridDataset(mk=mk_object)

#TODO: when providing mk=net_base, the grid is still refined, so mk is changed in-place (how to properly code this)
#interp bathy
data_bathy_interp = data_bathy_sel.interp(lon=xu_grid_uds.obj.mesh2d_node_x, lat=xu_grid_uds.obj.mesh2d_node_y).reset_coords(['lat','lon']) #interpolates lon/lat gebcodata to mesh2d_nNodes dimension #TODO: if these come from xu_grid_uds (without ojb), the mesh2d_node_z var has no ugrid accessor since the dims are lat/lon instead of mesh2d_nNodes
xu_grid_uds['mesh2d_node_z'] = data_bathy_interp.elevation.clip(max=10)

fig, ax = plt.subplots()
xu_grid_uds.grid.plot(ax=ax,linewidth=1)

fig, ax = plt.subplots()
xu_grid_uds.mesh2d_node_z.ugrid.plot(ax=ax,center=False)
#ctx.add_basemap(ax=ax, crs='EPSG:4326', attribution=False)

#write xugrid grid to netcdf
netfile  = os.path.join(dir_output, f'{model_name}_net.nc')
xu_grid_uds.ugrid.to_netcdf(netfile)
mdu.geometry.netfile = netfile #TODO: path is windows/unix dependent #TODO: providing os.path.basename(netfile) raises "ValidationError: 1 validation error for Geometry - netfile:   File: `C:\SnapVolumesTemp\MountPoints\{45c63495-0000-0000-0000-100000000000}\{79DE0690-9470-4166-B9EE-4548DC416BBD}\SVROOT\DATA\dfm_tools\tests\examples_workinprogress\Bonaire_net.nc` not found, skipped parsing." (wrong current directory)
#breakit

if 1:
    #%% new ext: initial and open boundary condition
    
    # CMEMS - download
    dir_output_data_cmems = os.path.join(dir_output_data, 'cmems')
    if not os.path.isdir(dir_output_data_cmems):
        os.mkdir(dir_output_data_cmems)
        
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
                                                 tstart = date_min,
                                                 tstop = date_max,
                                                 list_plifiles = [poly_file],
                                                 dir_sourcefiles_hydro = dir_output_data_cmems)
    
    #save new ext file
    ext_new.save(filepath=ext_file_new,path_style=path_style)
    
    
    #%% old ext
    
    # CMEMS - initial condition file
    ext_old = mb.preprocess_ini_cmems_to_nc(ext_old=ext_old,
                                            tSimStart=dt.datetime.strptime(date_min, '%Y-%m-%d'),
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
        #TODO does it loop correctly over all lists?
    
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
#TODO: missing keywords to be added to hydrolib-core: https://github.com/Deltares/HYDROLIB-core/issues/486.
#TODO: reference mdu: p:\11209231-003-bes-modellering\hydrodynamica\hackathon\simulations\run001_mapInterval_1800_trac_newGrid\Bonaire.mdu
#TODO: also compare settings in diafiles again
mdu.geometry.bedlevuni = 5
if 0: #TODO: this raises a segmentation error in dflowfm, report this
    mdu.geometry.kmx = 2
    mdu.geometry.layertype = 2
    mdu.geometry.numtopsig = 20
    mdu.geometry.sigmagrowthfactor = 1.2
else:
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
#mdu.numerics.mintimestepbreak = 0.1 #TODO: add to hydrolib-core
#mdu.numerics.keepstbndonoutflow = 1 #TODO: add to hydrolib-core

mdu.physics.tidalforcing = 1 #TODO: is 0 in diafile (maybe because of non-spherical grid)
mdu.physics.salinity = 1
mdu.physics.temperature = 5 #TODO: is 1 in diafile
mdu.physics.initialsalinity = 33.8
mdu.physics.temperature = 5
mdu.physics.initialtemperature = 29.3
mdu.physics.rhomean = 1023
mdu.physics.secchidepth = 4
#mdu.physics.salimax = 50 #TODO: add to hydrolib-core
#mdu.physics.tempmax = 50 #TODO: add to hydrolib-core
#mdu.physics.iniwithnudge = 2 #TODO: add to hydrolib-core #TODO: not implemented in oldextfile, do we want this?
#mdu.physics.jasfer3d = 1 #TODO: is 0 in file but should be 1 automatically, might cause model instability (might be caused by non spherical grid)

mdu.wind.icdtyp = 4
#TODO: also add cdbreakpoints for ERA5?
mdu.wind.rhoair = 1.2265
mdu.wind.relativewind = 0.5
mdu.wind.pavbnd = 101330

mdu.external_forcing.extforcefile = ext_file_old #TODO: is not written, but should be (with meteo)
#mdu.external_forcing.extforcefilenew = ext_file_new #TODO: relative paths? #TODO: extfile not found with path_style='unix': https://github.com/Deltares/HYDROLIB-core/issues/516
#mdu.external_forcing.windext = 1 #TODO: is set in reference mdu, happens automatically?

mdu.time.refdate = pd.Timestamp(ref_date).strftime('%Y%m%d')
mdu.time.tunit   = 'S'
mdu.time.dtmax   = 300
mdu.time.tstart  = (dt.datetime.strptime(date_min,'%Y-%m-%d') - dt.datetime.strptime(ref_date,'%Y-%m-%d')).total_seconds()
mdu.time.tstop   = (dt.datetime.strptime(date_max,'%Y-%m-%d') - dt.datetime.strptime(ref_date,'%Y-%m-%d')).total_seconds()
#mdu.time.Startdatetime  = pd.Timestamp(date_min).strftime('%Y%m%d%H%M%S') #TODO: add to hydrolib-core and remove mdu.time.tstart
#mdu.time.Stopdatetime   = pd.Timestamp(date_max).strftime('%Y%m%d%H%M%S') #TODO: add to hydrolib-core and remove mdu.time.tstop
#mdu.time.autotimestep = 3 #TODO: add to hydrolib-core

mdu.output.obsfile = [os.path.join(dir_output,x) for x in ['osm_beach_centroids_offset_bonaire.xyn','stations_obs.xyn']]
mdu.output.hisinterval = [60]
mdu.output.mapinterval = [1800]#[86400]


#%% export model
mdu.save(mdu_file,path_style=path_style)

#TODO: temporary workaround for extforcefilenew
mdu_contents = open(str(mdu_file),'r').readlines()
for iL, line in enumerate(mdu_contents):
    if line.startswith('[External Forcing]'):
        mdu_contents[iL] = line+'extForceFilenew = Bonaire_new.ext'+'\n'
with open(mdu_file,'w') as f:
    f.write(''.join(mdu_contents))

#TODO: if windows/unix extfile validation is fixed, use relative paths in .ext and in .mdu files
#TODO: investigate recurse saving


#%% model run
#TODO: hydrolib-core written mdu-file contains old keywords
"""
** WARNING: While reading 'Bonaire.mdu': keyword [numerics] qhrelax=0.01 was in file, but not used. Check possible typo.
** WARNING: While reading 'Bonaire.mdu': keyword [wind] windspeedbreakpoints=0.0 100.0 was in file, but not used. Check possible typo.
** WARNING: While reading 'Bonaire.mdu': keyword [output] wrishp_enc=0 was in file, but not used. Check possible typo.
** WARNING: While reading 'Bonaire.mdu': keyword [output] waterlevelclasses=0.0 was in file, but not used. Check possible typo.
** WARNING: While reading 'Bonaire.mdu': keyword [output] waterdepthclasses=0.0 was in file, but not used. Check possible typo.
"""

#TODO: when manually enabling extforcefilenew, we get the messages below, but maybe because of non-spherical grid or mdu settings?
"""
** INFO   : added node-based ghostcells:              0
** INFO   : partition_fixorientation_ghostlist: number of reversed flowlinks=              0
** INFO   : partition_fixorientation_ghostlist: number of reversed flowlinks=              0
** INFO   : partition_fixorientation_ghostlist: number of reversed flowlinks=              0
** INFO   : partition_fixorientation_ghostlist: number of reversed flowlinks=              0
** INFO   : Done partitioning model.
** INFO   : Done partitioning model.
** INFO   : Done partitioning model.
** INFO   : Done partitioning model.
** ERROR  : update_ghostboundvals: not all ghost boundary flowlinks are being updated
** ERROR  : update_ghostboundvals: not all ghost boundary flowlinks are being updated
** ERROR  : update_ghostboundvals: not all ghost boundary flowlinks are being updated
** ERROR  : update_ghostboundvals: not all ghost boundary flowlinks are being updated
Abort(1) on node 0 (rank 0 in comm 0): application called MPI_Abort(MPI_COMM_WORLD, 1) - process 0
Abort(1) on node 3 (rank 3 in comm 0): application called MPI_Abort(MPI_COMM_WORLD, 1) - process 3
Abort(1) on node 2 (rank 2 in comm 0): application called MPI_Abort(MPI_COMM_WORLD, 1) - process 2
"""
#When doing a sequential run, the error is clearer (seems to be cartesian/spherical related)
"""
Dimr [2023-04-11 18:42:59.558644] #0 >> kernel: ...
** WARNING: Variable 'air_pressure' in NetCDF file '/p/11209231-003-bes-modellering/hydrodynamica/hackathon/preprocessing/ModelBuilderOutput_JV/era5_msl_u10n_v10n_chnk_20221101to20221103_ERA5.nc requires 'projection_x_coordinate' and 'projection_y_coordinate'.
Dimr [2023-04-11 18:42:59.558677] #0 >> kernel: Variable 'air_pressure' in NetCDF file '/p/11209231-003-bes-modellering/hydrodynamica/hackathon/preprocessing/ModelBuilderOutput_JV/era5_msl_u10n_v10n_chnk_20221101to20221103_ERA5.nc requires 'projection_x_coordinate' and 'projection_y_coordinate'.
** ERROR  : flow_initexternalforcings: Error while initializing quantity: airpressure_windx_windy_charnock Check preceding log lines for details.
"""


#TODO: missing faces in quickplot, maybe because instable model
