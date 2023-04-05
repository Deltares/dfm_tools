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
mb = dfmt #from dfm_tools import dfm_modelbuilder as mb
import hydrolib.core.dflowfm as hcdfm
#import shutil
import xarray as xr
import pandas as pd

## input
model_name = 'Bonaire'
dir_output = r'p:\11209231-003-bes-modellering\hydrodynamica\hackathon\preprocessing\ModelBuilderOutput_JV'
path_style = 'unix' # windows / unix

#TODO: reference run in: p:\11209231-003-bes-modellering\hydrodynamica\hackathon\simulations\run001_mapInterval_1800_waq_newGrid

# domain
lon_min, lon_max, lat_min, lat_max = -68.55, -67.9, 11.8, 12.6, 

#dates as understood by pandas.period_range(). ERA5 has freq='M' (month) and CMEMS has freq='D' (day)
date_min = '2022-11-01'
date_max = '2022-11-03'
ref_date = '2022-01-01'

# meteo (ERA5)
# option 0 (testing)  = [airpressure]
# option 1 (2D-model) = [airpressure, windx, windy, charnock]
# option 2 (3D-model) = [airpressure, windx, windy, charnock, dewpoint, airtemperature, cloudiness, solarradiation, longwaveradiation, rainfall, evaporation]
ERA5_meteo_option = 0

## process
# make dirs
if not os.path.isdir(dir_output):
    os.mkdir(dir_output)
dir_output_data = os.path.join(dir_output, 'data')
if not os.path.isdir(dir_output_data):
    os.mkdir(dir_output_data)

#%% mdu
mdu_file = os.path.join(dir_output, f'{model_name}.mdu')
ext_file_new = os.path.join(dir_output, f'{model_name}_bc.ext')

mdu = hcdfm.FMModel()
ext_bnd = hcdfm.ExtModel()
#TODO: initialize old extmodel

#%% GEBCO bathymetry and Grid generation

#select and plot bathy
file_nc_bathy = r'p:\metocean-data\open\GEBCO\2021\GEBCO_2021.nc'
data_bathy = xr.open_dataset(file_nc_bathy)
data_bathy_sel = data_bathy.sel(lon=slice(lon_min-1,lon_max+1),lat=slice(lat_min-1,lat_max+1))

#TODO: grid generation/refinement based on bathy still to be improved in meshkernel (https://github.com/Deltares/dfm_tools/issues/234), replace if fixed
mk_object = mb.make_basegrid(lon_min, lon_max, lat_min, lat_max) #TODO: should be sperical, but is cartesian >> is_geographic keywork does not work yet

#refine
min_face_size = 200/(40075*1000/360) #convert meters to degrees
mb.refine_basegrid(mk=mk_object, data_bathy_sel=data_bathy_sel, min_face_size=min_face_size) #TODO: min_face_size is now in degrees instead of meters

#cutcells
file_ldb = r'p:\11209231-003-bes-modellering\hydrodynamica\hackathon\preprocessing\grid\coastline.pli'
mb.remove_grid_withldb(mk=mk_object,file_ldb=file_ldb) #TODO: does not seem to have effect

#convert to xugrid and interp bathy #TODO: when providing mk=net_base, the grid is still refined, so mk is changed in-place (how to properly code this)
xu_grid_ds = mb.xugrid_interp_bathy(mk=mk_object, data_bathy_sel=data_bathy_sel) #TODO: now has topology_dimension=1 to avoid model crash, but grid is pink in interacter

#write xugrid grid to netcdf
netfile  = os.path.join(dir_output, f'{model_name}_net.nc')
xu_grid_ds.to_netcdf(netfile)

mdu.geometry.netfile = netfile #TODO: path is windows/unix dependent #TODO: providing os.path.basename(netfile) raises "ValidationError: 1 validation error for Geometry - netfile:   File: `C:\SnapVolumesTemp\MountPoints\{45c63495-0000-0000-0000-100000000000}\{79DE0690-9470-4166-B9EE-4548DC416BBD}\SVROOT\DATA\dfm_tools\tests\examples_workinprogress\Bonaire_net.nc` not found, skipped parsing." (wrong current directory)

#%% initial and open boundary condition
# TODO create polyline (currently this file is prerequisite)
poly_file = os.path.join(dir_output, 'bnd.pli')

#%% CMEMS
# CMEMS - download
dir_output_data_cmems = os.path.join(dir_output_data, 'cmems')
if not os.path.isdir(dir_output_data_cmems):
    os.mkdir(dir_output_data_cmems)
    
mb.download_meteodata_oceandata(
    longitude_min = lon_min, longitude_max = lon_max, latitude_min = lat_min, latitude_max = lat_max,
    model = 'CMEMS',
    date_min = date_min, date_max = date_max,
    varlist = ['so','thetao','uo','vo','zos'],
    dir_output = dir_output_data_cmems)

# CMEMS - boundary conditions file (.bc) (and add to ext_bnd)
ext_bnd = mb.preprocess_interpolate_nc_to_bc(ext_bnd=ext_bnd, #TODO: is 12h on day before also downloaded/selected?
                                             refdate_str = 'minutes since '+ref_date+' 00:00:00 +00:00',
                                             dir_output = dir_output,
                                             model = 'CMEMS',
                                             tstart = date_min,
                                             tstop = date_max,
                                             list_plifiles = [poly_file],
                                             dir_sourcefiles_hydro = dir_output_data_cmems)

#save new ext file
ext_bnd.save(filepath=ext_file_new,path_style=path_style)


#%% old ext

# CMEMS - initial condition file
mb.preprocess_ini_cmems_to_nc(tSimStart=dt.datetime.strptime(date_min, '%Y-%m-%d'),
                              dir_data=dir_output_data_cmems,
                              dir_out=dir_output)

# ERA5 - download
dir_output = os.path.join(dir_output_data, 'ERA5')
if not os.path.isdir(dir_output):
    os.mkdir(dir_output)
    
if ERA5_meteo_option == 0:
    varlist = [['msl']]
elif ERA5_meteo_option == 1:
    varlist = [['msl','u10n','v10n','chnk']]
elif ERA5_meteo_option == 2:
    varlist = [['msl','u10n','v10n','chnk'],['d2m','t2m','tcc'],['ssr','strd'],['mer','mtpr']]
    #TODO does it loop correctly over all lists?  
           
mb.download_meteodata_oceandata(
    longitude_min = lon_min, longitude_max = lon_max, latitude_min = lat_min, latitude_max = lat_max,
    model = 'ERA5',
    date_min = date_min, date_max = date_max,
    varlist = varlist, # check variables_dict in dfmt.download_ERA5() for valid names
    dir_output = dir_output_data)

# ERA5 meteo - convert to netCDF for usage in Delft3D FM
mb.preprocess_merge_meteofiles(
        mode = 'ERA5',
        varkey_list = varlist,
        dir_data = dir_output_data,
        dir_output = dir_output,
        time_slice = slice(date_min, date_max))

#%% obs points
#TODO: add obs points

#%% .mdu settings
#TODO: missing keywords to be added to hydrolib-core: https://github.com/Deltares/HYDROLIB-core/issues/486.
#TODO: reference mdu: p:\11209231-003-bes-modellering\hydrodynamica\hackathon\simulations\run001_mapInterval_1800_trac_newGrid\Bonaire.mdu
#TODO: also compare settings in diafiles again
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
#mdu.numerics.mintimestepbreak = 0.1 #TODO: add to hydrolib-core
#mdu.numerics.keepstbndonoutflow = 1 #TODO: add to hydrolib-core

mdu.physics.tidalforcing = 1 #TODO: is 0 in diafile
mdu.physics.salinity = 1
mdu.physics.temperature = 5 #TODO: is 1 in diafile
mdu.physics.initialsalinity = 33.8
mdu.physics.temperature = 5
mdu.physics.initialtemperature = 29.3
mdu.physics.rhomean = 1023
mdu.physics.secchidepth = 4
#mdu.physics.salimax = 50 #TODO: add to hydrolib-core
#mdu.physics.tempmax = 50 #TODO: add to hydrolib-core
#mdu.physics.iniwithnudge = 2 #TODO: add to hydrolib-core
#mdu.physics.jasfer3d = 1 #TODO: is 0 in file but should be 1 automatically, might cause model instability (caused by non spherical grid)

mdu.wind.icdtyp = 4
#TODO: also add cdbreakpoints?
mdu.wind.rhoair = 1.2265
mdu.wind.relativewind = 0.5
mdu.wind.pavbnd = 101330

#mdu.external_forcing.extforcefile = f'{model_name}.ext' #TODO: is not written, but should be (with meteo)
#mdu.external_forcing.extforcefilenew = ext_file_new #TODO: relative paths?
#mdu.external_forcing.windext = 1 #TODO: is set in reference mdu, automatically?

#TODO: ValidationError: 5 validation errors for ExternalForcing
#extforcefilenew -> Bonaire_bc.ext -> boundary -> 0 -> forcingFile
#  File: `\\directory.intra\Project\p\11209231-003-bes-modellering\hydrodynamica\hackathon\preprocessing\ModelBuilderOutput_JV\data\cmems\waterlevelbnd_bnd_CMEMS.bc` not found, skipped parsing. (type=value_error)
#TODO: this is since paths are unix and we are validating on a windows system

mdu.time.refdate = pd.Timestamp(ref_date).strftime('%Y%m%d')
mdu.time.tunit   = 'S'
mdu.time.dtmax   = 300
mdu.time.tstart  = (dt.datetime.strptime(date_min,'%Y-%m-%d') - dt.datetime.strptime(ref_date,'%Y-%m-%d')).total_seconds() #TODO: replace with timestring keyword
mdu.time.tstop   = (dt.datetime.strptime(date_max,'%Y-%m-%d') - dt.datetime.strptime(ref_date,'%Y-%m-%d')).total_seconds()
#mdu.time.Startdatetime  = pd.Timestamp(date_min).strftime('%Y%m%d%H%M%S') #TODO: add to hydrolib-core
#mdu.time.Stopdatetime   = pd.Timestamp(date_max).strftime('%Y%m%d%H%M%S') #TODO: add to hydrolib-core
#mdu.time.autotimestep = 3 #TODO: add to hydrolib-core

#mdu.output.obsfile = #TODO: add obsfile
mdu.output.hisinterval = [60]
mdu.output.mapinterval = [1800]#[86400]

#%% export model
mdu.save(mdu_file,path_style=path_style)
# TODO: relative paths in .ext
#TODO: model is instable (2000m waterlevel on second mapfile timestep)
"""
#TODO: hydrolib-core written mdu-file contains old keywords
** WARNING: While reading 'Bonaire.mdu': keyword [numerics] qhrelax=0.01 was in file, but not used. Check possible typo.
** WARNING: While reading 'Bonaire.mdu': keyword [wind] windspeedbreakpoints=0.0 100.0 was in file, but not used. Check possible typo.
** WARNING: While reading 'Bonaire.mdu': keyword [output] wrishp_enc=0 was in file, but not used. Check possible typo.
** WARNING: While reading 'Bonaire.mdu': keyword [output] waterlevelclasses=0.0 was in file, but not used. Check possible typo.
** WARNING: While reading 'Bonaire.mdu': keyword [output] waterdepthclasses=0.0 was in file, but not used. Check possible typo.
"""

#TODO: QP warning:
"""
Unable to merge the partitions
Too many outputs requested.  Most likely cause is missing [] around left hand side that has a comma separated list expansion.
In private\nc_interpret>nc_mapmerge at line 1878
In private\nc_interpret at line 91
"""