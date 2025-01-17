# -*- coding: utf-8 -*-
"""
Created on Thu Dec  8 20:54:39 2022

@author: veenstra

"""

import meshkernel
import xarray as xr
import matplotlib.pyplot as plt
plt.close('all')
import contextily as ctx
import dfm_tools as dfmt
import datetime as dt

#general settings
lon_min,lon_max = -5,2
lat_min,lat_max = 48.5,51.2
lon_res,lat_res = 0.2,0.2

figsize = (10,4)
crs = 'EPSG:4326'


"""
Make a regular (potentially rotated) rectilinear grid. First generate a curvilinear grid than convert the curvilinear grid into unstructured grid. The steps are the following:
- curvilinear_make_uniform, see the following notebook: https://github.com/Deltares/MeshKernelPy/blob/AddCurvilinearGridSupport/docs/examples/04_curvilineargrid_basics.ipynb
- curvilinear_convert_to_mesh2d: https://github.com/Deltares/MeshKernelPy/blob/118cb4953c4e95d5b18ed283bb37f391134b2bb2/meshkernel/meshkernel.py#L1399 

"""

#TODO: maybe use make_basegrid and refine_basegrid functions from dfmt.meshkernel_helpers
# Create an instance of MakeGridParameters and set the values
make_grid_parameters = meshkernel.MakeGridParameters(angle=0.0, #TODO: does non-zero result in an orthogonal spherical grid?
                                                     origin_x=lon_min,
                                                     origin_y=lat_min,
                                                     upper_right_x=lon_max, #TODO: angle and upper_right_x cannot be combined: https://github.com/Deltares/MeshKernelPy/issues/74
                                                     upper_right_y=lat_max,
                                                     block_size_x=lon_res,
                                                     block_size_y=lat_res)

mk = meshkernel.MeshKernel(projection=meshkernel.ProjectionType.SPHERICAL)
mk.curvilinear_compute_rectangular_grid_on_extension(make_grid_parameters)
mk.curvilinear_convert_to_mesh2d() #convert to ugrid/mesh2d

mesh2d_basegrid = mk.mesh2d_get() #in case of curvi grid: mk.curvilinear_convert_to_mesh2d()
fig, ax = plt.subplots(figsize=figsize)
mesh2d_basegrid.plot_edges(ax,linewidth=0.8)
ctx.add_basemap(ax=ax, crs=crs, attribution=False)


"""
Mesh refinement in MeshKernelPy with bathymetry samples and plot result
"""
# connect to bathymetry dataset
file_nc_bathy = r'p:\metocean-data\open\GEBCO\2021\GEBCO_2021.nc'
data_bathy = xr.open_dataset(file_nc_bathy).elevation
# alternatively you can connect to ETOPO, for which there is also a 15s (15 arcseconds) resolution dataset available
# file_nc_bathy = "https://www.ngdc.noaa.gov/thredds/dodsC/global/ETOPO2022/30s/30s_surface_elev_netcdf/ETOPO_2022_v1_30s_N90W180_surface.nc"
# data_bathy = xr.open_dataset(file_nc_bathy).z

# subset bathy to area of interest 
data_bathy_sel = data_bathy.sel(lon=slice(lon_min-1, lon_max+1), lat=slice(lat_min-1, lat_max+1))

# plot bathymetry
fig, ax = plt.subplots(figsize=figsize)
data_bathy_sel.plot(ax=ax, center=False)
ctx.add_basemap(ax=ax, crs=crs, attribution=False)

#convert bathy data to GriddedSamples
lon_np = data_bathy_sel.lon.to_numpy()
lat_np = data_bathy_sel.lat.to_numpy()
values_np = data_bathy_sel.to_numpy().flatten()
gridded_samples = meshkernel.GriddedSamples(x_coordinates=lon_np,y_coordinates=lat_np,values=values_np)

#refinement
mesh_refinement_parameters = meshkernel.MeshRefinementParameters(#refine_intersected=False, #TODO: what does this do?
                                                                 #use_mass_center_when_refining=False, #TODO: what does this do?
                                                                 #account_for_samples_outside_face=False, #outsidecell argument for --refine?
                                                                 #directional_refinement=False, #TODO: what does this do?
                                                                 #max_refinement_iterations=5, #TODO: default is 10
                                                                 min_edge_size=1000, #always in meters
                                                                 refinement_type=meshkernel.RefinementType(1), #Wavecourant/1,
                                                                 connect_hanging_nodes=True, #set to False to do multiple refinement steps (e.g. for multiple regions)
                                                                 smoothing_iterations=2, #TODO: consider 3 (gtsm value)
                                                                 max_courant_time=150) #TODO: consider 200 (gtsm value), is this value read?

mk.mesh2d_refine_based_on_gridded_samples(gridded_samples=gridded_samples,
                                          mesh_refinement_params=mesh_refinement_parameters,
                                          use_nodal_refinement=True) #TODO: what does this do?

mesh2d_refinedgrid = mk.mesh2d_get()
fig, ax = plt.subplots(figsize=figsize)
mesh2d_refinedgrid.plot_edges(ax,linewidth=0.8)
ctx.add_basemap(ax=ax, crs=crs, attribution=False)


"""
delete (landward) part of grid with polygon and plot result
"""

print('deleting landward grid with coastlines (takes a while)...')
dtstart = dt.datetime.now()
dfmt.meshkernel_delete_withcoastlines(mk=mk, res='h', min_area=100)
print(f'done in {(dt.datetime.now()-dtstart).total_seconds():.2f} sec')

mesh2d_noland = mk.mesh2d_get()
fig, ax = plt.subplots(figsize=figsize)
mesh2d_noland.plot_edges(ax,linewidth=0.8)
dfmt.plot_coastlines(ax=ax, res='h', min_area=100)
ctx.add_basemap(ax=ax, crs=crs, attribution=False)


"""
remove disconnected parts of grid (from non-contiguous to contiguous)
"""

mk.mesh2d_remove_disconnected_regions()

mesh2d_contiguous = mk.mesh2d_get()
fig, ax = plt.subplots(figsize=figsize)
mesh2d_contiguous.plot_edges(ax,linewidth=0.8)
dfmt.plot_coastlines(ax=ax, res='h', min_area=100)
ctx.add_basemap(ax=ax, crs=crs, attribution=False)


"""
convert meshkernel grid to xugrid, interp bathymetry, plot and save to *_net.nc
"""

xu_grid_uds = dfmt.meshkernel_to_UgridDataset(mk, crs=crs)

fig, ax = plt.subplots(figsize=figsize)
xu_grid_uds.grid.plot(ax=ax,linewidth=0.8) #TODO: maybe make uds instead of ds (but then bathy interpolation goes wrong)
dfmt.plot_coastlines(ax=ax, res='h', min_area=100)
ctx.add_basemap(ax=ax, crs=crs, attribution=False)

#interp bathy
data_bathy_interp = data_bathy_sel.interp(lon=xu_grid_uds.obj.mesh2d_node_x, lat=xu_grid_uds.obj.mesh2d_node_y) #interpolates lon/lat gebcodata to mesh2d_nNodes dimension #TODO: if these come from xu_grid_uds, the mesh2d_node_z var has no ugrid accessor since the dims are lat/lon instead of mesh2d_nNodes
xu_grid_uds['mesh2d_node_z'] = data_bathy_interp.clip(max=10)
#TODO: alternatively do this with TODO: mk.mesh2d_averaging_interpolation() or mk.mesh2d_triangulation_interpolation() >> bilinear would be faster

fig, ax = plt.subplots(figsize=figsize)
xu_grid_uds.mesh2d_node_z.ugrid.plot(ax=ax,center=False)
ctx.add_basemap(ax=ax, crs=crs, attribution=False)

#write xugrid grid to netcdf
netfile = 'englishchannel_net.nc'
xu_grid_uds.ugrid.to_netcdf(netfile)



"""
get illegal cells
"""

illegalcells_gdf = dfmt.meshkernel_get_illegalcells(mk=mk)


"""
others
"""


#TODO: small flow links in resulting grid

#TODO: cleanup grid necessary?
print('mk.mesh2d_get_obtuse_triangles_mass_centers()')
print(mk.mesh2d_get_obtuse_triangles_mass_centers().values)
print('mk.mesh2d_get_orthogonality()')
print(mk.mesh2d_get_orthogonality().values.max())
print('mk.mesh2d_get_hanging_edges()')
print(mk.mesh2d_get_hanging_edges())
mk.mesh2d_delete_hanging_edges()
mk_ortho = mk.mesh2d_get_orthogonality()
mk_ortho_vals = mk_ortho.values
mk_ortho_vals[mk_ortho_vals==mk_ortho.geometry_separator] = 0 # or np.nan, but results in invisible edges (grey would be better)
xu_grid_uds['orthogonality'] = xr.DataArray(mk_ortho_vals, dims=xu_grid_uds.grid.edge_dimension)
fig, ax = plt.subplots(figsize=figsize)
xu_grid_uds['orthogonality'].ugrid.plot(robust=True)
