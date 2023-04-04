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

## input
model_name = 'Bonaire'
dir_output = r'p:\11209231-003-bes-modellering\hydrodynamica\hackathon\preprocessing\ModelBuilderOutput_JV'
dir_output_main = dir_output

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
#TODO

#%% GEBCO bathymetry and Grid generation
dir_output = os.path.join(dir_output_data, 'grid_bathymetry')
if not os.path.isdir(dir_output):
    os.mkdir(dir_output)

#TODO: create netfile (currently this file is prerequisite)
#TODO: grid generation/refinement based on bathy still to be improved in meshkernel (https://github.com/Deltares/dfm_tools/issues/234), replace if fixed
mb.GEBCO2asc(lon_min-1,    lon_max+1,    lat_min-1,    lat_max+1,    suffix='full', dir_out=dir_output) #TODO: use inpolygon functionality from this function also in todevelop function
mb.GEBCO2asc(lon_min+0.05, lon_max-0.05, lat_min+0.05, lat_max-0.05, suffix='refine', dir_out=dir_output)
netfile  = os.path.join(dir_output_main, f'{model_name}_net.nc')
mb.make_basegrid(lon_min, lon_max, lat_min, lat_max, name=model_name, filepath=netfile.replace('_net.nc','_basegrid_net.nc'))
#TODO mb.refine_basegrid(bathydataset_refine='GEBCO2021_res1d240deg_res464m_refine.asc', dxmin_refine='0200', filepath_in=netfile, dir_out=dir_output_main)
#TODO mb.interp_bathy()

#mdu.geometry.netfile = netfile.replace('_net.nc','_basegrid_net.nc') #TODO: "KeyError: 'An attribute for "mesh2d_face_x" was not found in nc file. Expected "mesh2d_face_x"'" (happens with xugrid netfile, was ugrid netfile before)
mdu.geometry.netfile = netfile

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

# mdu.external_forcing.extforcefile = model_name+'.ext'
# mdu.external_forcing.extforcefilenew = model_name+'_bc.ext'

tstart = dt.datetime.strptime(date_max,'%Y-%m-%d') 
tstop = dt.datetime.strptime(date_max,'%Y-%m-%d') 
mdu.time.refdate = dt.datetime.strptime(ref_date,'%Y-%m-%d').strftime('%Y%m%d')
mdu.time.tunit   = 'S'
mdu.time.tstart  = (dt.datetime.strptime(date_min,'%Y-%m-%d') - dt.datetime.strptime(ref_date,'%Y-%m-%d')).total_seconds() #TODO: replace with timestring keyword
mdu.time.tstop   = (dt.datetime.strptime(date_max,'%Y-%m-%d') - dt.datetime.strptime(ref_date,'%Y-%m-%d')).total_seconds()

mdu.output.hisinterval = [60]
mdu.output.mapinterval = [86400]


"""
#TODO: to be added to hydrolib-core: https://github.com/Deltares/HYDROLIB-core/issues/486
#TODO now Manually add:
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
mdu.save(mdu_file)
#TODO ext.save(filepath=ext_file)
# TODO: relative paths in .ext
