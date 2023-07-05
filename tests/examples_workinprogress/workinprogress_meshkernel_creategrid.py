# -*- coding: utf-8 -*-
"""
Created on Thu Dec  8 20:54:39 2022

@author: veenstra

"""

import meshkernel
import xarray as xr
import matplotlib.pyplot as plt
plt.close('all')
import numpy as np
import contextily as ctx
import dfm_tools as dfmt

#TODO: maybe use make_basegrid and refine_basegrid functions from dfmt.meshkernel_helpers

#general settings
lon_min,lon_max = -6,2
lat_min,lat_max = 48.5,51.2
lon_res,lat_res = 0.2,0.2
is_geographic = True #True for spherical grids

figsize = (10,4)
crs = 'EPSG:4326'


"""
Make a regular (potentially rotated) rectilinear grid. First generate a curvilinear grid than convert the curvilinear grid into unstructured grid. The steps are the following:
- curvilinear_make_uniform, see the following notebook: https://github.com/Deltares/MeshKernelPy/blob/AddCurvilinearGridSupport/docs/examples/04_curvilineargrid_basics.ipynb
- curvilinear_convert_to_mesh2d: https://github.com/Deltares/MeshKernelPy/blob/118cb4953c4e95d5b18ed283bb37f391134b2bb2/meshkernel/meshkernel.py#L1399 

"""

# Create an instance of MakeGridParameters and set the values
make_grid_parameters = meshkernel.MakeGridParameters(angle=0.0, #TODO: does non-zero result in an orthogonal spherical grid?
                                                     origin_x=lon_min,
                                                     origin_y=lat_min,
                                                     upper_right_x=lon_max, #TODO: angle and upper_right_x cannot be combined: https://github.com/Deltares/MeshKernelPy/issues/74
                                                     upper_right_y=lat_max,
                                                     block_size_x=lon_res,
                                                     block_size_y=lat_res)

grid_in_pol = False
# If a polygon is provided it will be used in the generation of the curvilinear grid. The polygon must be closed
if grid_in_pol: #can be used instead of origin_x/origin_y and num_x/num_y
    pol_x = np.array([-6,-4,0,-6], dtype=np.double)
    pol_y = np.array([48,51,49.5,48], dtype=np.double)
    geometry_list = meshkernel.GeometryList(pol_x, pol_y)
else:
    geometry_list = None


mk = meshkernel.MeshKernel(is_geographic=is_geographic)
mk.curvilinear_make_uniform_on_extension(make_grid_parameters) #TODO: geometry_list is possible with curvilinear_make_uniform, but not for this function: https://github.com/Deltares/MeshKernelPy/issues/76
mk.curvilinear_convert_to_mesh2d() #convert to ugrid/mesh2d
mesh2d_basegrid = mk.mesh2d_get() #in case of curvi grid: mk.curvilinear_convert_to_mesh2d()
fig, ax = plt.subplots(figsize=figsize)
mesh2d_basegrid.plot_edges(ax,linewidth=0.8)
ctx.add_basemap(ax=ax, crs=crs, attribution=False)


"""
Mesh refinement in MeshKernelPy with bathymetry samples and plot result
"""
#select and plot bathy
file_nc_bathy = r'p:\metocean-data\open\GEBCO\2021\GEBCO_2021.nc'
data_bathy = xr.open_dataset(file_nc_bathy)
data_bathy_sel = data_bathy.sel(lon=slice(lon_min-1,lon_max+1),lat=slice(lat_min-1,lat_max+1))

fig, ax = plt.subplots(figsize=figsize)
data_bathy_sel.elevation.plot(ax=ax, center=False)
ctx.add_basemap(ax=ax, crs=crs, attribution=False)

#convert bathy data to GriddedSamples
lon_np = data_bathy_sel.lon.to_numpy()
lat_np = data_bathy_sel.lat.to_numpy()
values_np = data_bathy_sel.elevation.to_numpy().flatten().astype('float') #TODO: astype to avoid "TypeError: incompatible types, c_short_Array_74880 instance instead of LP_c_double instance"
gridded_samples = meshkernel.GriddedSamples(x_coordinates=lon_np,y_coordinates=lat_np,values=values_np) #TODO: does not result in refinement


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

#line, = ax.plot([], [],'o-') # empty line
# #linebuilder = dfmt.LineBuilder(line) #TODO: this makes it possible to interactively click a line in the bedlevel figure. Use linebuilder.line_array as alternative line_array
# delete_pol = np.array([[ 1.91741935, 49.76580645],
#                         [ 0.20387097, 49.9       ],
#                         [-0.25032258, 48.71290323],
#                         [ 1.92774194, 48.59935484]])
dfmt.meshkernel_delete_withcoastlines(mk=mk, res='h', min_area=1000)


mesh2d_noland = mk.mesh2d_get()
fig, ax = plt.subplots(figsize=figsize)
mesh2d_noland.plot_edges(ax,linewidth=0.8)
# xlim,ylim = ax.get_xlim(),ax.get_ylim() #get x/ylims before ldb plotting changes it
# for iP, pol_del in enumerate(pol_ldb_list):
#     ax.plot(pol_del['x'],pol_del['y'],'-r')
# ax.set_xlim(xlim)
# ax.set_ylim(ylim)
dfmt.plot_coastlines(ax=ax, res='h', min_area=1000)
ctx.add_basemap(ax=ax, crs=crs, attribution=False)



"""
convert meshkernel grid to xugrid, interp bathymetry, plot and save to *_net.nc
"""

xu_grid_uds = dfmt.meshkernel_to_UgridDataset(mk, remove_noncontiguous=True, is_geographic=is_geographic, crs=crs) #TODO: put remove_noncontiguous in meshkernel?: https://github.com/Deltares/MeshKernelPy/issues/44

fig, ax = plt.subplots(figsize=figsize)
xu_grid_uds.grid.plot(ax=ax) #TODO: maybe make uds instead of ds (but then bathy interpolation goes wrong)
ctx.add_basemap(ax=ax, crs=crs, attribution=False)

#interp bathy
data_bathy_interp = data_bathy_sel.interp(lon=xu_grid_uds.obj.mesh2d_node_x, lat=xu_grid_uds.obj.mesh2d_node_y).reset_coords(['lat','lon']) #interpolates lon/lat gebcodata to mesh2d_nNodes dimension #TODO: if these come from xu_grid_uds, the mesh2d_node_z var has no ugrid accessor since the dims are lat/lon instead of mesh2d_nNodes
xu_grid_uds['mesh2d_node_z'] = data_bathy_interp.elevation.clip(max=10)
#TODO: alternatively do this with TODO: mk.mesh2d_averaging_interpolation() or mk.mesh2d_triangulation_interpolation() >> bilinear would be faster

fig, ax = plt.subplots(figsize=figsize)
xu_grid_uds.mesh2d_node_z.ugrid.plot(ax=ax,center=False)
ctx.add_basemap(ax=ax, crs=crs, attribution=False)

#write xugrid grid to netcdf
xu_grid_uds.ugrid.to_netcdf('test_net.nc')

#TODO: update https://github.com/Deltares/dfm_tools/issues/217

#TODO: there is (was) a link missing, maybe due to ldb?

