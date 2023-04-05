# -*- coding: utf-8 -*-
"""
Created on Tue Apr  4 16:12:56 2023

@author: groenenb, laan_st, veenstra
"""

import datetime as dt
import os
from pathlib import Path
import matplotlib.pyplot as plt
plt.close('all')
import dfm_tools as dfmt
mb = dfmt #from dfm_tools import dfm_modelbuilder as mb
import hydrolib.core.dflowfm as hcdfm
import shutil
import xarray as xr

## input
model_name = 'Bonaire'
dir_output = r'p:\11209231-003-bes-modellering\hydrodynamica\hackathon\preprocessing\ModelBuilderOutput_JV'
dir_output_main = dir_output
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
mdu_file = os.path.join(dir_output_main, f'{model_name}.mdu')
ext_file = os.path.join(dir_output_main, f'{model_name}_bc.ext')

mdu = hcdfm.FMModel()
ext = hcdfm.ExtModel()

#%% land-sea-boundary
#TODO get cutting with ldb from https://github.com/Deltares/dfm_tools/blob/main/tests/examples_workinprogress/workinprogress_meshkernel_creategrid.py

#%% GEBCO bathymetry and Grid generation
dir_output = os.path.join(dir_output_data, 'grid_bathymetry')
if not os.path.isdir(dir_output):
    os.mkdir(dir_output)

#select and plot bathy
file_nc_bathy = r'p:\metocean-data\open\GEBCO\2021\GEBCO_2021.nc'
data_bathy = xr.open_dataset(file_nc_bathy)
data_bathy_sel = data_bathy.sel(lon=slice(lon_min-1,lon_max+1),lat=slice(lat_min-1,lat_max+1))

#TODO: grid generation/refinement based on bathy still to be improved in meshkernel (https://github.com/Deltares/dfm_tools/issues/234), replace if fixed
net_base = mb.make_basegrid(lon_min, lon_max, lat_min, lat_max) #TODO: should be sperical, but is cartesian >> is_geographic keywork does not work yet

#refine
min_face_size = 200/(40075*1000/360) #convert meters to degrees
net_refined = mb.refine_basegrid(mk=net_base, data_bathy_sel=data_bathy_sel, min_face_size=min_face_size) #TODO: min_face_size is now in degrees instead of meters

#convert to xugrid and interp bathy
xu_grid_uds = mb.xugrid_interp_bathy(mk=net_refined, data_bathy_sel=data_bathy_sel)

#write xugrid grid to netcdf
netfile  = os.path.join(dir_output_main, f'{model_name}_net.nc')
xu_grid_uds.ugrid.to_netcdf(netfile)


mdu.geometry.netfile = netfile #TODO: path is windows/unix dependent #TODO: providing os.path.basename(netfile) raises "ValidationError: 1 validation error for Geometry - netfile:   File: `C:\SnapVolumesTemp\MountPoints\{45c63495-0000-0000-0000-100000000000}\{79DE0690-9470-4166-B9EE-4548DC416BBD}\SVROOT\DATA\dfm_tools\tests\examples_workinprogress\Bonaire_net.nc` not found, skipped parsing." (wrong current directory)

#%% initial and open boundary condition
# TODO create polyline (currently this file is prerequisite)
poly_file = os.path.join(dir_output_main,'bnd.pli')

#%% CMEMS
# CMEMS - download
dir_output = os.path.join(dir_output_data, 'cmems')
if not os.path.isdir(dir_output):
    os.mkdir(dir_output)
    
mb.download_meteodata_oceandata(
    longitude_min = lon_min, longitude_max = lon_max, latitude_min = lat_min, latitude_max = lat_max,
    model = 'CMEMS',
    date_min = date_min, date_max = date_max,
    varlist = ['so','thetao','uo','vo','zos'],
    dir_output = dir_output)

# CMEMS - initial condition file
mb.preprocess_ini_cmems_to_nc(tSimStart = dt.datetime.strptime(date_min, '%Y-%m-%d'),
    dir_data  = dir_output,
    dir_out = dir_output_main)

# CMEMS - boundary conditions file (.bc)
mb.preprocess_interpolate_nc_to_bc(
    refdate_str = 'minutes since '+ref_date+' 00:00:00 +00:00',
    dir_output = dir_output,
    model = 'CMEMS',
    tstart = date_min,
    tstop = date_max,
    list_plifiles = [poly_file],
    dir_sourcefiles_hydro = dir_output)

# copy .bc to main output dir
for filename in os.listdir(dir_output):
    if filename.endswith('.bc'):
        shutil.copy(os.path.join(dir_output, filename), os.path.join(dir_output_main, filename))

# add to .ext
cmems_ext_file = os.path.join(dir_output_data, 'cmems', 'example_bnd.ext')
boundary_object = hcdfm.ExtModel(Path(cmems_ext_file))
ext.boundary.append(boundary_object)

#save ext file
#ext.save(filepath=ext_file,path_style=path_style) #TODO: "AttributeError: 'ExtModel' object has no attribute '_to_section'"


#%% ERA5
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
    dir_output = dir_output)

# ERA5 meteo - convert to netCDF for usage in Delft3D FM
mb.preprocess_merge_meteofiles(
        mode = 'ERA5',
        varkey_list = varlist,
        dir_data = dir_output,
        dir_output = dir_output_main,
        time_slice = slice(date_min, date_max))

#%% obs points
#TODO: add obs points

#%% .mdu settings
mdu.geometry.bedlevuni = 5
mdu.geometry.kmx = 2
mdu.geometry.layertype = 2
mdu.geometry.numtopsig = 20
mdu.geometry.sigmagrowthfactor = 1.2

mdu.numerics.tlfsmo = 86400
mdu.numerics.izbndpos = 1

mdu.physics.tidalforcing = 1
mdu.physics.salinity = 1
mdu.physics.temperature = 5
mdu.physics.initialsalinity = 33.8
mdu.physics.temperature = 5
mdu.physics.initialtemperature = 29.3
mdu.physics.rhomean = 1023
mdu.physics.secchidepth = 4

mdu.wind.icdtyp = 4
mdu.wind.rhoair = 1.2265
mdu.wind.relativewind = 0.5
mdu.wind.pavbnd = 101330

#mdu.external_forcing.extforcefile = f'{model_name}.ext' #TODO: is not written, but should be (with meteo)
#mdu.external_forcing.extforcefilenew = ext_file #f'{model_name}_bc.ext' #TODO: relative paths?

tstart = dt.datetime.strptime(date_max,'%Y-%m-%d') 
tstop = dt.datetime.strptime(date_max,'%Y-%m-%d') 
mdu.time.refdate = dt.datetime.strptime(ref_date,'%Y-%m-%d').strftime('%Y%m%d')
mdu.time.tunit   = 'S'
mdu.time.tstart  = (dt.datetime.strptime(date_min,'%Y-%m-%d') - dt.datetime.strptime(ref_date,'%Y-%m-%d')).total_seconds() #TODO: replace with timestring keyword
mdu.time.tstop   = (dt.datetime.strptime(date_max,'%Y-%m-%d') - dt.datetime.strptime(ref_date,'%Y-%m-%d')).total_seconds()
#mdu.time.tstart  = (pd.Datetime(date_min) - pd.Datetime(ref_date)).total_seconds() #TODO: replace with timestring keyword
#mdu.time.tstop   = (pd.Datetime(date_max) - pd.Datetime(ref_date)).total_seconds()

mdu.output.hisinterval = [60]
mdu.output.mapinterval = [86400]


"""
#TODO: to be added to hydrolib-core: https://github.com/Deltares/HYDROLIB-core/issues/486. Now Manually add:
[geometry]
netFile                    = 5_bathy_net.nc

Numtopsiguniform                    = 1                                                                   # Number of sigma layers in top of z-layer model
Dztop                               = 5.0                                                                 # Z-layer thickness of layers above level Dztopuniabovez
Floorlevtoplay                      = -5.0                                                                # Floor level of top layer
Dztopuniabovez                      = -100.0                                                              # Above level Dztopuniabovez layers will have uniform Dztop, SigmaGrowthFactor below this level
Keepzlayeringatbed                  = 2                                                                   #  bedlayerthickness = zlayerthickness at bed 0 or 1

[numerics]
MinTimestepBreak = 0.1

[Physics]
salimax                       = 50
tempmax                       = 50
IniWithNudge                  = 2

[time]
autotimestep = 3

[External Forcing]
    extForceFile    = Bonaire.ext    # Old format for external forcings file *.ext, link with tim/cmp-format boundary conditions specification.
    extForceFileNew = Bonaire_bc.ext # New format for external forcings file *.ext, link with bcformat boundary conditions specification.

"""

#%% export model
mdu.save(mdu_file,path_style=path_style)
# TODO: relative paths in .ext

"""
#TODO: hydrolib-core written mdu-file contains old keywords
** WARNING: While reading 'Bonaire.mdu': keyword [numerics] qhrelax=0.01 was in file, but not used. Check possible typo.
** WARNING: While reading 'Bonaire.mdu': keyword [wind] windspeedbreakpoints=0.0 100.0 was in file, but not used. Check possible typo.
** WARNING: While reading 'Bonaire.mdu': keyword [output] wrishp_enc=0 was in file, but not used. Check possible typo.
** WARNING: While reading 'Bonaire.mdu': keyword [output] waterlevelclasses=0.0 was in file, but not used. Check possible typo.
** WARNING: While reading 'Bonaire.mdu': keyword [output] waterdepthclasses=0.0 was in file, but not used. Check possible typo.

"""

#TODO: fix dflowfm error when running this model (might be grid related, since it happens on partitioning process. opening grid in interacter shows it to be cartesian)
#ug_get_meshgeom, #12, ierr=0
#forrtl: severe (174): SIGSEGV, segmentation fault occurred

