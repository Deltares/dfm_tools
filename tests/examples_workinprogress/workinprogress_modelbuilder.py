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
import contextily as ctx

## input
model_name = 'Bonaire'
dir_output = r'p:\11209231-003-bes-modellering\hydrodynamica\hackathon\preprocessing\ModelBuilderOutput_JV'
path_style = 'unix' # windows / unix
overwrite = False # used for downloading of forcing data. Always set to True when changing the domain
paths_relative = False #TODO: currently only works with path_style='windows' (same OS as IDE)
is_geographic = True
crs = 'EPSG:4326'

inisaltem = True #initialsalinity/initialtemperature gives 33.8ppt uniform and sal instabilities right from the start of the model run. Proper way seems to be with iniwithnudge=2 and nudge_salinity_temperature, which gives ini sal/tem indeed but also instable. Providing nudge_salinity_temperature and iniwithnudge=0 gives more stable model but also inisal is 33.8 (not spatially varying) (is same as False maybe?)
#TODO: salinity instable, also waterlevel and velocity magnitude are instable at northeast side of island (latter is with incorrect ordering/selection in extfile)
"""
** INFO   :  Min. salinity limited, number of cells Limmin =           20
** INFO   :  Min. salinity limited, min =  -1.037033177733807E-005
"""


#TODO: files that are not created in this script: obsfiles, dimr.xml and the submit script (and GEBCO+GSHHS datasets)
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
os.makedirs(dir_output, exist_ok=True)
dir_output_data = os.path.join(dir_output, 'data')
os.makedirs(dir_output_data, exist_ok=True)


#%%plifile generation
pli_polyfile = dfmt.generate_bndpli(lon_min, lon_max, lat_min, lat_max, dlon=bnd_dlon_dlat, dlat=bnd_dlon_dlat, name=f'{model_name}_bnd')
#TODO: generate pli from mk with mk_object.mesh2d_get_mesh_boundaries_as_polygons()
poly_file = os.path.join(dir_output, f'{model_name}.pli')
pli_polyfile.save(poly_file)


#%% grid generation
# GEBCO bathymetry and Grid generation
#select and plot bathy
file_nc_bathy = r'p:\metocean-data\open\GEBCO\2021\GEBCO_2021.nc'
data_bathy = xr.open_dataset(file_nc_bathy)
data_bathy_sel = data_bathy.sel(lon=slice(lon_min-1,lon_max+1),lat=slice(lat_min-1,lat_max+1))

#TODO: grid generation and bathy-refinement is still to be improved in meshkernel (https://github.com/Deltares/dfm_tools/issues/234)
mk_object = dfmt.make_basegrid(lon_min, lon_max, lat_min, lat_max, is_geographic=is_geographic)

#refine
min_edge_size = 300 #in meters
dfmt.refine_basegrid(mk=mk_object, data_bathy_sel=data_bathy_sel, min_edge_size=min_edge_size)

#cutcells
dfmt.meshkernel_delete_withcoastlines(mk=mk_object, res='h') #TODO: write used coastline to ldbfile?
#TODO: illegalcells.pol necessary?

#TODO: cleanup grid necessary?
# print('mk_object.mesh2d_get_obtuse_triangles_mass_centers()')
# print(mk_object.mesh2d_get_obtuse_triangles_mass_centers().values)
# print('mk_object.mesh2d_get_orthogonality()')
# print(mk_object.mesh2d_get_orthogonality().values.max()) #TODO: couple back to uds, currently ordering mismatch: https://github.com/Deltares/MeshKernelPy/issues/72
# print('mk_object.mesh2d_get_hanging_edges()')
# print(mk_object.mesh2d_get_hanging_edges())
# mk_object.mesh2d_delete_hanging_edges()

#convert to xugrid
xu_grid_uds = dfmt.meshkernel_to_UgridDataset(mk=mk_object)

#TODO: temporary fix until code from issue-fix is in main branch: https://github.com/Deltares/dfm_tools/issues/421
from netCDF4 import default_fillvals
import numpy as np
attribute_dict = {
            'name': 'WGS84',
            'epsg': np.array(4326, dtype=int),
            'grid_mapping_name': 'latitude_longitude',
            'longitude_of_prime_meridian': np.array(0.0, dtype=float),
            'semi_major_axis': np.array(6378137.0, dtype=float),
            'semi_minor_axis': np.array(6356752.314245, dtype=float),
            'inverse_flattening': np.array(298.257223563, dtype=float),
            'EPSG_code': 'EPSG:4326',
            }
xu_grid_uds['wgs84'] = xr.DataArray(np.array(default_fillvals['i4'],dtype=int),dims=(),attrs=attribute_dict)

#interp bathy
data_bathy_interp = data_bathy_sel.interp(lon=xu_grid_uds.obj.mesh2d_node_x, lat=xu_grid_uds.obj.mesh2d_node_y).reset_coords(['lat','lon']) #interpolates lon/lat gebcodata to mesh2d_nNodes dimension #TODO: if these come from xu_grid_uds (without ojb), the mesh2d_node_z var has no ugrid accessor since the dims are lat/lon instead of mesh2d_nNodes
xu_grid_uds['mesh2d_node_z'] = data_bathy_interp.elevation.clip(max=10)

fig, ax = plt.subplots()
xu_grid_uds.grid.plot(ax=ax,linewidth=1)
ctx.add_basemap(ax=ax, crs=crs, attribution=False)
dfmt.plot_coastlines(ax=ax, crs=crs)

fig, ax = plt.subplots()
xu_grid_uds.mesh2d_node_z.ugrid.plot(ax=ax,center=False)
ctx.add_basemap(ax=ax, crs=crs, attribution=False)
dfmt.plot_coastlines(ax=ax, crs=crs)

#write xugrid grid to netcdf
netfile  = os.path.join(dir_output, f'{model_name}_net.nc')
xu_grid_uds.ugrid.to_netcdf(netfile)


#%% new ext: initial and open boundary condition
ext_file_new = os.path.join(dir_output, f'{model_name}_new.ext')
ext_new = hcdfm.ExtModel()

# FES2014 tidal components bc file
file_bc_basename = os.path.basename(poly_file).replace('.pli','')
ForcingModel_object = dfmt.interpolate_tide_to_bc(tidemodel='FES2014', file_pli=poly_file, component_list=None) # tidemodel: FES2014, FES2012, EOT20, GTSM4.1preliminary
file_bc_out = os.path.join(dir_output,f'tide_{file_bc_basename}_FES2014.bc')
ForcingModel_object.save(filepath=file_bc_out)
boundary_object = hcdfm.Boundary(quantity='waterlevelbnd', #the FM quantity for tide is also waterlevelbnd
                                 locationfile=poly_file,
                                 forcingfile=ForcingModel_object)
ext_new.boundary.append(boundary_object)

# CMEMS - download
dir_output_data_cmems = os.path.join(dir_output_data, 'cmems')
os.makedirs(dir_output_data_cmems, exist_ok=True)
for varkey in ['so','thetao','uo','vo','zos']:
    dfmt.download_CMEMS(credentials=None, #credentials=['username','password'], or create "%USERPROFILE%/CMEMS_credentials.txt" with username on line 1 and password on line 2. Register at: https://resources.marine.copernicus.eu/registration-form'
                        varkey=varkey,
                        longitude_min=lon_min, longitude_max=lon_max, latitude_min=lat_min, latitude_max=lat_max,
                        date_min=date_min, date_max=date_max,
                        dir_output=dir_output_data_cmems, file_prefix='cmems_', overwrite=overwrite)

# CMEMS - boundary conditions file (.bc) (and add to ext_bnd)
list_quantities = ['waterlevelbnd','salinitybnd','temperaturebnd','uxuyadvectionvelocitybnd'] # when supplying two waterlevelbnds to FM (tide and steric) with other quantities in between, dimrset>=2.24.00 is required or else "ERROR  : update_ghostboundvals: not all ghost boundary flowlinks are being updated" is raised (https://issuetracker.deltares.nl/browse/UNST-7011). Two waterlevelbnds need to share same physical plifile in order to be appended (https://issuetracker.deltares.nl/browse/UNST-5320).
ext_new = mb.cmems_nc_to_bc(ext_bnd=ext_new,
                            refdate_str=f'minutes since {ref_date} 00:00:00 +00:00',
                            dir_output=dir_output,
                            list_quantities=list_quantities,
                            tstart=date_min,
                            tstop=date_max, 
                            file_pli=poly_file,
                            dir_pattern=os.path.join(dir_output_data_cmems,'cmems_{ncvarname}_*.nc'))

#save new ext file
ext_new.save(filepath=ext_file_new,path_style=path_style)

    
#%% old ext

# CMEMS - initial condition file
ext_file_old = os.path.join(dir_output, f'{model_name}_old.ext')
ext_old = hcdfm.ExtOldModel()

if inisaltem:
    ext_old = mb.preprocess_ini_cmems_to_nc(ext_old=ext_old,
                                            tstart=date_min,
                                            dir_data=dir_output_data_cmems,
                                            dir_out=dir_output)

# ERA5 - download
dir_output_data_era5 = os.path.join(dir_output_data,'ERA5')
os.makedirs(dir_output_data_era5, exist_ok=True)
    
if ERA5_meteo_option == 1:
    varlist_list = [['msl','u10n','v10n','chnk']]
elif ERA5_meteo_option == 2:
    varlist_list = [['msl','u10n','v10n','chnk'],['d2m','t2m','tcc'],['ssr','strd'],['mer','mtpr']]

for varlist in varlist_list:
    for varkey in varlist:
        dfmt.download_ERA5(varkey, 
                           longitude_min=lon_min, longitude_max=lon_max, latitude_min=lat_min, latitude_max=lat_max,
                           date_min=date_min, date_max=date_max,
                           dir_output=dir_output_data_era5, overwrite=overwrite)

# ERA5 meteo - convert to netCDF for usage in Delft3D FM
ext_old = mb.preprocess_merge_meteofiles_era5(ext_old=ext_old,
                                              varkey_list = varlist_list,
                                              dir_data = dir_output_data_era5,
                                              dir_output = dir_output,
                                              time_slice = slice(date_min, date_max))

ext_old.save(filepath=ext_file_old,path_style=path_style)


#%% initialize mdu file and update settings
mdu_file = os.path.join(dir_output, f'{model_name}.mdu')
mdu = hcdfm.FMModel()

mdu.geometry.netfile = netfile #TODO: path is windows/unix dependent #TODO: providing os.path.basename(netfile) raises "ValidationError: 1 validation error for Geometry - netfile:   File: `C:\SnapVolumesTemp\MountPoints\{45c63495-0000-0000-0000-100000000000}\{79DE0690-9470-4166-B9EE-4548DC416BBD}\SVROOT\DATA\dfm_tools\tests\examples_workinprogress\Bonaire_net.nc` not found, skipped parsing." (wrong current directory)
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
#TODO: investigate keywords herman for stable model: keepzlayeringatbed=1 en baroczlaybed=1 en keepzlaybedvol=1 en zerozbndinflow/zerozbndadvectionvelocity=1

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
if inisaltem:
    mdu.physics.iniwithnudge = 2 #TODO: commented in oldextfile in reference run, initial sal/tem profiles from deep layer were used instead (not yet derived, but 3D inifields also do not have an effect)

mdu.wind.icdtyp = 4
#mdu.wind.cdbreakpoints = [0.025] #TODO: overwritten by spacevarying charnock from ERA5
mdu.wind.rhoair = 1.2265
mdu.wind.relativewind = 0.5
mdu.wind.pavbnd = 101330

mdu.external_forcing.extforcefile = ext_file_old
mdu.external_forcing.extforcefilenew = ext_new #TODO: extfile not found with path_style='unix': https://github.com/Deltares/HYDROLIB-core/issues/516, but workaround is with extobject instead of extpath.

mdu.time.refdate = pd.Timestamp(ref_date).strftime('%Y%m%d')
mdu.time.tunit = 'S'
mdu.time.dtmax = 30
mdu.time.startdatetime = pd.Timestamp(date_min).strftime('%Y%m%d%H%M%S')
mdu.time.stopdatetime = pd.Timestamp(date_max).strftime('%Y%m%d%H%M%S')
mdu.time.autotimestep = 3

mdu.output.obsfile = [os.path.join(dir_output,x) for x in ['osm_beach_centroids_offset_bonaire.xyn','stations_obs.xyn','test_obs.xyn']]
mdu.output.hisinterval = [60]
mdu.output.mapinterval = [1800]#[86400]
mdu.output.rstinterval = [0] #TODO: default is 0.0, but this translates to [0.0 tstart tstop] instead of [0]: https://github.com/Deltares/HYDROLIB-core/issues/525#issuecomment-1505111924
mdu.output.statsinterval = [3600]
#TODO: disable many outputs (preferrably change many default values): https://github.com/Deltares/HYDROLIB-core/issues/525#issuecomment-1505111924


#%% export model
mdu.save(mdu_file,path_style=path_style)

#TODO: if windows/unix newextfile validation is fixed, use relative paths in .ext and in .mdu files (this is a workaround that only works for same-OS so windows paths)
if paths_relative:
    for filename in [mdu_file,ext_file_old,ext_file_new]:
        with open(filename, 'r') as file :
            filedata = file.read()
        filedata = filedata.replace(dir_output.replace('\\','/')+'/', '') #dir_output or os.path.dirname(mdu_file)
        with open(filename, 'w') as file:
            file.write(filedata)

