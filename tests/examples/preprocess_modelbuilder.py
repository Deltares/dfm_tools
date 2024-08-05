# -*- coding: utf-8 -*-
"""
Created on Tue Apr  4 16:12:56 2023

@author: groenenb, laan_st, veenstra
"""

import os
import matplotlib.pyplot as plt
plt.close('all')
import dfm_tools as dfmt
import hydrolib.core.dflowfm as hcdfm
import xarray as xr
import pandas as pd
import contextily as ctx

## input
model_name = 'Bonaire'
dir_output = r'p:\11209231-003-bes-modellering\hydrodynamica\hackathon\preprocessing\ModelBuilderOutput_JV2'
path_style = 'windows' # windows / unix, making relative paths only works when path_style is equal to os
overwrite = False # used for downloading of forcing data. Always set to True when changing the domain
crs = 'EPSG:4326'

#TODO: salinity instable, also waterlevel and velocity magnitude are instable at northeast side of island (latter is with incorrect ordering/selection in extfile)
"""
** INFO   :  Min. salinity limited, number of cells Limmin =           20
** INFO   :  Min. salinity limited, min =  -1.037033177733807E-005
"""


#TODO: files that are not created in this script: obsfiles, dimr.xml and the submit script (and GEBCO+GSHHS datasets)
#TODO: reference run in: p:\11209231-003-bes-modellering\hydrodynamica\hackathon\simulations\run001_mapInterval_1800\
#TODO: also compare settings to p:\11208054-004-dcsm-fm\models\3D_DCSM-FM\2013-2017\B05_hydrolib_JV\DCSM-FM_0_5nm.mdu (e.g. tlfSmo)
# domain and resolution
lon_min, lon_max, lat_min, lat_max = -68.55, -67.9, 11.8, 12.6
dxy = 0.05

#dates as understood by pandas.period_range(). ERA5 has freq='M' (month) and CMEMS has freq='D' (day)
date_min = '2022-11-01'
date_max = '2022-11-03'
ref_date = '2022-01-01'

# meteo (ERA5)
# option 1 (2D-model) = [airpressure, windx, windy, charnock]
# option 2 (3D-model) = [airpressure, windx, windy, charnock, dewpoint, airtemperature, cloudiness, solarradiation, longwaveradiation, rainfall, evaporation]
ERA5_meteo_option = 2

# make dirs
os.makedirs(dir_output, exist_ok=True)
dir_output_data = os.path.join(dir_output, 'data')
os.makedirs(dir_output_data, exist_ok=True)


#%% grid generation and refinement with GEBCO bathymetry
# connect to bathymetry dataset
file_nc_bathy = r'p:\metocean-data\open\GEBCO\2021\GEBCO_2021.nc'
data_bathy = xr.open_dataset(file_nc_bathy).elevation
# alternatively you can connect to ETOPO, for which there is also a 15s (15 arcseconds) resolution dataset available
# file_nc_bathy = "https://www.ngdc.noaa.gov/thredds/dodsC/global/ETOPO2022/30s/30s_surface_elev_netcdf/ETOPO_2022_v1_30s_N90W180_surface.nc"
# data_bathy = xr.open_dataset(file_nc_bathy).z

# subset bathy to area of interest 
data_bathy_sel = data_bathy.sel(lon=slice(lon_min-1, lon_max+1), lat=slice(lat_min-1, lat_max+1))


#TODO: grid generation and bathy-refinement is still to be improved in meshkernel (https://github.com/Deltares/dfm_tools/issues/234)
mk_object = dfmt.make_basegrid(lon_min, lon_max, lat_min, lat_max, dx=dxy, dy=dxy, crs=crs)

# generate plifile from grid extent and coastlines
bnd_gdf = dfmt.generate_bndpli_cutland(mk=mk_object, res='h', buffer=0.01)
bnd_gdf_interp = dfmt.interpolate_bndpli(bnd_gdf,res=0.06)
pli_polyfile = dfmt.geodataframe_to_PolyFile(bnd_gdf_interp, name=f'{model_name}_bnd')
poly_file = os.path.join(dir_output, f'{model_name}.pli')
pli_polyfile.save(poly_file)

#refine
min_edge_size = 300 #in meters
dfmt.refine_basegrid(mk=mk_object, data_bathy_sel=data_bathy_sel, min_edge_size=min_edge_size)

#cutcells
dfmt.meshkernel_delete_withcoastlines(mk=mk_object, res='h') #TODO: write used coastline to ldbfile?
#TODO: illegalcells.pol can be acquired with dfmt.meshkernel_get_illegalcells()

#convert to xugrid
xu_grid_uds = dfmt.meshkernel_to_UgridDataset(mk=mk_object, crs=crs)

#TODO: cleanup grid necessary?
# print('mk_object.mesh2d_get_obtuse_triangles_mass_centers()')
# print(mk_object.mesh2d_get_obtuse_triangles_mass_centers().values)
# print('mk_object.mesh2d_get_orthogonality()')
# print(mk_object.mesh2d_get_orthogonality().values.max())
# print('mk_object.mesh2d_get_hanging_edges()')
# print(mk_object.mesh2d_get_hanging_edges())
# mk_object.mesh2d_delete_hanging_edges()
mk_ortho = mk_object.mesh2d_get_orthogonality()
mk_ortho_vals = mk_ortho.values
mk_ortho_vals[mk_ortho_vals==mk_ortho.geometry_separator] = 0 # or np.nan, but results in invisible edges (grey would be better)
xu_grid_uds['orthogonality'] = xr.DataArray(mk_ortho_vals, dims=xu_grid_uds.grid.edge_dimension)
xu_grid_uds['orthogonality'].ugrid.plot()

#interp bathy
data_bathy_interp = data_bathy_sel.interp(lon=xu_grid_uds.obj.mesh2d_node_x, lat=xu_grid_uds.obj.mesh2d_node_y)
xu_grid_uds['mesh2d_node_z'] = data_bathy_interp.clip(max=10)

fig, ax = plt.subplots()
xu_grid_uds.grid.plot(ax=ax,linewidth=1)
ctx.add_basemap(ax=ax, crs=crs, attribution=False)
dfmt.plot_coastlines(ax=ax, crs=crs)

fig, ax = plt.subplots()
xu_grid_uds.mesh2d_node_z.ugrid.plot(ax=ax,center=False)
ctx.add_basemap(ax=ax, crs=crs, attribution=False)
dfmt.plot_coastlines(ax=ax, crs=crs)
bnd_gdf_interp.plot(ax=ax,color='r')

#write xugrid grid to netcdf
netfile  = os.path.join(dir_output, f'{model_name}_net.nc')
xu_grid_uds.ugrid.to_netcdf(netfile)

#%% new ext: initial and open boundary condition
ext_file_new = os.path.join(dir_output, f'{model_name}_new.ext')
ext_new = hcdfm.ExtModel()

# FES2014 tidal components bc file
tidemodel = 'EOT20' # tidemodel: FES2014, FES2012, EOT20, GTSMv4.1, GTSMv4.1_opendap
dfmt.interpolate_tide_to_bc(ext_new=ext_new, tidemodel=tidemodel, file_pli=poly_file, component_list=None)

# CMEMS - download
# you can also add WAQ variables like 'no3' and 'phyc'
# check dfmt.get_conversion_dict() for an overview of parameter/quantity names
dir_output_data_cmems = os.path.join(dir_output_data, 'cmems')
os.makedirs(dir_output_data_cmems, exist_ok=True)
for varkey in ['zos','so','thetao','uo','vo','no3']:#,'phyc']: # TODO: phyc not available in reanalysis: https://github.com/Deltares/dfm_tools/issues/847
    dfmt.download_CMEMS(varkey=varkey,
                        longitude_min=lon_min, longitude_max=lon_max, latitude_min=lat_min, latitude_max=lat_max,
                        date_min=date_min, date_max=date_max,
                        dir_output=dir_output_data_cmems, file_prefix='cmems_', overwrite=overwrite)

# CMEMS - boundary conditions file (.bc) (and add to ext_bnd)
# you can also add WAQ variables like 'tracerbndNO3' and 'tracerbndPON1'
# check dfmt.get_conversion_dict() for an overview of parameter/quantity names
# when supplying two waterlevelbnds to FM (tide and steric) with other quantities in between, dimrset>=2.24.00 is required
# or else "ERROR  : update_ghostboundvals: not all ghost boundary flowlinks are being updated" is raised (https://issuetracker.deltares.nl/browse/UNST-7011).
# Two waterlevelbnds need to share same physical plifile in order to be appended (https://issuetracker.deltares.nl/browse/UNST-5320).
list_quantities = ['waterlevelbnd','salinitybnd','temperaturebnd','uxuyadvectionvelocitybnd','tracerbndNO3']#,'tracerbndPON1']
dir_pattern = os.path.join(dir_output_data_cmems,'cmems_{ncvarname}_*.nc')
ext_new = dfmt.cmems_nc_to_bc(ext_bnd=ext_new,
                              refdate_str=f'minutes since {ref_date} 00:00:00 +00:00',
                              dir_output=dir_output,
                              list_quantities=list_quantities,
                              tstart=date_min,
                              tstop=date_max, 
                              file_pli=poly_file,
                              dir_pattern=dir_pattern)

#save new ext file
ext_new.save(filepath=ext_file_new,path_style=path_style)

    
#%% old ext
ext_file_old = os.path.join(dir_output, f'{model_name}_old.ext')
ext_old = hcdfm.ExtOldModel()

# CMEMS - initial conditions
ext_old = dfmt.cmems_nc_to_ini(ext_old=ext_old,
                               dir_output=dir_output,
                               list_quantities=list_quantities,
                               tstart=date_min,
                               dir_pattern=dir_pattern)

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
ext_old = dfmt.preprocess_merge_meteofiles_era5(ext_old=ext_old,
                                                varkey_list = varlist_list,
                                                dir_data = dir_output_data_era5,
                                                dir_output = dir_output,
                                                time_slice = slice(date_min, date_max))

ext_old.save(filepath=ext_file_old,path_style=path_style)


#%% initialize mdu file and update settings
mdu_file = os.path.join(dir_output, f'{model_name}.mdu')
mdu = hcdfm.FMModel()

mdu.geometry.netfile = netfile #TODO: path is windows/unix dependent #TODO: providing os.path.basename(netfile) raises "ValidationError: 1 validation error for Geometry - netfile:   File: `C:\SnapVolumesTemp\MountPoints\{45c63495-0000-0000-0000-100000000000}\{79DE0690-9470-4166-B9EE-4548DC416BBD}\SVROOT\DATA\dfm_tools\tests\examples\Bonaire_net.nc` not found, skipped parsing." (wrong current directory)
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
if 'nudge_salinity_temperature' in [x.quantity for x in ext_old.forcing]:
    # initialsalinity/initialtemperature gives 33.8ppt uniform (so is not read)
    # the only possible way is with iniwithnudge=2 and nudge_salinity_temperature in old extfile.
    # for other variables we can use initialtracerbndNO3 for instance
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

#TODO: workaround to make paths relative until https://github.com/Deltares/HYDROLIB-core/issues/532 is implemented
#TODO: currently only works with path_style='windows' (same OS as IDE)
dfmt.make_paths_relative(mdu_file)

nproc = 1 # number of processes
dimrset_folder = None # r"c:\Program Files\Deltares\Delft3D FM Suite 2023.03 HMWQ\plugins\DeltaShell.Dimr\kernels" #alternatively r"p:\d-hydro\dimrset\weekly\2.25.17.78708"
dfmt.create_model_exec_files(file_mdu=mdu_file, nproc=nproc, dimrset_folder=dimrset_folder)

