# -*- coding: utf-8 -*-
"""
Created on Thu Dec  8 20:54:39 2022

@author: veenstra

"""

import meshkernel
import xarray as xr
import xugrid as xu
import matplotlib.pyplot as plt
plt.close('all')
import numpy as np
import contextily as ctx
import dfm_tools as dfmt
import pandas as pd
from pathlib import Path
try: #0.3.1 release
    from hydrolib.core.io.polyfile.models import PolyFile
except: #main branch and next release #TODO: move to easy imports after https://github.com/Deltares/HYDROLIB-core/issues/410
    from hydrolib.core.dflowfm.polyfile.models import PolyFile

#general settings
lon_min,lon_max = -6,2
lat_min,lat_max = 48.5,51.2
lon_res,lat_res = 0.2,0.2
num_x = int(np.ceil((lon_max-lon_min)/lon_res))
num_y = int(np.ceil((lat_max-lat_min)/lat_res))

figsize = (10,4)
crs = 'EPSG:4326'


"""
Make a regular (potentially rotated) rectilinear grid. First generate a curvilinear grid than convert the curvilinear grid into unstructured grid. The steps are the following:
- curvilinear_make_uniform, see the following notebook: https://github.com/Deltares/MeshKernelPy/blob/AddCurvilinearGridSupport/docs/examples/04_curvilineargrid_basics.ipynb
- curvilinear_convert_to_mesh2d: https://github.com/Deltares/MeshKernelPy/blob/118cb4953c4e95d5b18ed283bb37f391134b2bb2/meshkernel/meshkernel.py#L1399 

"""

# Create an instance of MakeGridParameters and set the values
make_grid_parameters = meshkernel.MakeGridParameters() #TODO: contains default values that are maybe not intuitive
make_grid_parameters.num_columns = num_x
make_grid_parameters.num_rows = num_y
make_grid_parameters.angle = 0.0 #TODO: does non-zero result in an orthogonal spherical grid?
#make_grid_parameters.block_size = 0.0
make_grid_parameters.origin_x = lon_min
make_grid_parameters.origin_y = lat_min
make_grid_parameters.block_size_x = lon_res
make_grid_parameters.block_size_y = lat_res

grid_in_pol = False
# A polygon must to be provided. If empty it will not be used. If a polygon is provided it will be used in the generation of the curvilinear grid. The polygon must be closed
if grid_in_pol: #can be used instead of origin_x/origin_y and num_x/num_y
    pol_x = np.array([-6,-4,0,-6], dtype=np.double)
    pol_y = np.array([48,51,49.5,48], dtype=np.double)
else:
    pol_x = np.empty(0, dtype=np.double)
    pol_y = np.empty(0, dtype=np.double)
geometry_list = meshkernel.GeometryList(pol_x, pol_y)

mk1 = meshkernel.MeshKernel()
mk1.curvilinear_make_uniform(make_grid_parameters, geometry_list) #TODO: make geometry_list argument optional
mk1.curvilinear_convert_to_mesh2d() #convert to ugrid/mesh2d
mesh2d_grid1 = mk1.mesh2d_get() #in case of curvi grid: mk.curvilinear_convert_to_mesh2d()
fig, ax = plt.subplots(figsize=figsize)
mesh2d_grid1.plot_edges(ax,linewidth=1.2)
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

#convert bathy data to geomlist
samp_x,samp_y = np.meshgrid(data_bathy_sel.lon.to_numpy(),data_bathy_sel.lat.to_numpy())
samp_z = data_bathy_sel.elevation.to_numpy().astype(float) #TODO: without .astype(float), meshkernelpy generates "TypeError: incompatible types, c_short_Array_27120 instance instead of LP_c_double instance"
samp_x = samp_x.ravel()
samp_y = samp_y.ravel()
samp_z = samp_z.ravel()
geomlist = meshkernel.GeometryList(x_coordinates=samp_x, y_coordinates=samp_y, values=samp_z) #TODO: does not check if lenghts of input array is equal (samp_z[1:])

#refinement
mesh_refinement_parameters = meshkernel.MeshRefinementParameters(refine_intersected=False, #TODO: provide defaults for several arguments, so less arguments are required
                                                                 use_mass_center_when_refining=False, #TODO: what does this do?
                                                                 min_face_size=0.01, #TODO: size in meters would be more convenient
                                                                 refinement_type=meshkernel.RefinementType(1), #Wavecourant/1,
                                                                 connect_hanging_nodes=False, #set to False to do multiple refinement steps (e.g. for multiple regions)
                                                                 account_for_samples_outside_face=False, #outsidecell argument for --refine?
                                                                 max_refinement_iterations=5,
                                                                 ) #TODO: missing the arguments dtmax (necessary?), hmin (min_face_size but then in meters instead of degrees), smoothiters (currently refinement is patchy along coastlines), spherical 1/0 (necessary?)

mk2 = meshkernel.MeshKernel()
mk2.curvilinear_make_uniform(make_grid_parameters, geometry_list)
mk2.curvilinear_convert_to_mesh2d() #convert to ugrid/mesh2d
mk2.mesh2d_refine_based_on_samples(samples=geomlist,
                                   relative_search_radius=0.5, ##TODO: bilin interp is preferred, but this is currently not supported (samples have to be ravelled)
                                   minimum_num_samples=3,
                                   mesh_refinement_params=mesh_refinement_parameters,
                                   )

mesh2d_grid2 = mk2.mesh2d_get()
fig, ax = plt.subplots(figsize=figsize)
mesh2d_grid2.plot_edges(ax,linewidth=1.2)
ctx.add_basemap(ax=ax, crs=crs, attribution=False)


"""
delete (landward) part of grid with polygon and plot result
"""
delete_with_ldb = True
if delete_with_ldb: #landboundary file has to contain closed polygons
    print('reading ldb')
    file_ldb = r'p:\1230882-emodnet_hrsm\global_tide_surge_model\trunk\gtsm5.0\GSHHS_intermediate.ldb'
    pol_ldb = PolyFile(Path(file_ldb))
    print('converting ldb')
    pol_ldb_list = [dfmt.pointlike_to_DataFrame(x) for x in pol_ldb.objects] #TODO, this is really slow, speed up possible?
    pol_ldb_list = [x for x in pol_ldb_list if len(x)>1000] #filter only large polygons for performance
    print('done with ldb')
else:
    #line, = ax.plot([], [],'o-') # empty line
    #linebuilder = dfmt.LineBuilder(line) #this makes it possible to interactively click a line in the bedlevel figure. Use linebuilder.line_array as alternative line_array
    delete_pol = np.array([[ 1.91741935, 49.76580645],
                           [ 0.20387097, 49.9       ],
                           [-0.25032258, 48.71290323],
                           [ 1.92774194, 48.59935484]])
    delete_pol = np.concatenate([delete_pol,delete_pol[[0],:]],axis=0) #close polygon
    pol_ldb_list = [pd.DataFrame(delete_pol,columns=['x','y'])]

for iP, pol_del in enumerate(pol_ldb_list): #TODO: also possible without loop? >> geometry_separator=-999.9 so that value can be used to concat polygons
    delete_pol_geom = meshkernel.GeometryList(x_coordinates=pol_del['x'].to_numpy(), y_coordinates=pol_del['y'].to_numpy()) #TODO: .copy()/to_numpy() makes the array contiguous in memory, which is necessary for meshkernel.mesh2d_delete()
    mk2.mesh2d_delete(geometry_list=delete_pol_geom, 
                      delete_option=meshkernel.DeleteMeshOption(2), #ALL_COMPLETE_FACES/2: Delete all faces of which the complete face is inside the polygon
                      invert_deletion=False)

mesh2d_grid3 = mk2.mesh2d_get()
fig, ax = plt.subplots(figsize=figsize)
mesh2d_grid3.plot_edges(ax,linewidth=1.2)
xlim,ylim = ax.get_xlim(),ax.get_ylim() #get x/ylims before ldb plotting changes it
for iP, pol_del in enumerate(pol_ldb_list):
    ax.plot(pol_del['x'],pol_del['y'],'-r')
ax.set_xlim(xlim)
ax.set_ylim(ylim)
ctx.add_basemap(ax=ax, crs=crs, attribution=False)


"""
convert meshkernel grid to xugrid, plot and save to *_net.nc
"""
#generate face_node_connectivity array to make conversion to xugrid possible #TODO: move this to meshkernelpy internal? >> is/willbe available in xugrid
num_faces = len(mesh2d_grid2.nodes_per_face)
num_face_nodes_max = np.max(mesh2d_grid2.nodes_per_face)
face_node_connectivity = np.full(shape=(num_faces,num_face_nodes_max), dtype=np.int32, fill_value=-1)
index_in_mesh2d = 0
for face_index, num_face_nodes in enumerate(mesh2d_grid2.nodes_per_face):
    range_face_node_index = index_in_mesh2d+num_face_nodes
    face_node_connectivity[face_index,range(num_face_nodes)] = mesh2d_grid2.face_nodes[index_in_mesh2d:range_face_node_index]
    index_in_mesh2d = index_in_mesh2d + num_face_nodes

#convert to xugrid grid and plot
xu_grid = xu.Ugrid2d(node_x=mesh2d_grid2.node_x,
                     node_y=mesh2d_grid2.node_y,
                     fill_value=-1,
                     face_node_connectivity=face_node_connectivity)

fig, ax = plt.subplots(figsize=figsize)
if 0: #TODO: disabled since it generates "AttributeError: 'Ugrid2d' object has no attribute 'plot'" (was not the case before) https://github.com/Deltares/xugrid/issues/39
    xu_grid.plot(ax=ax)
    ctx.add_basemap(ax=ax, crs=crs, attribution=False)

#write xugrid grid to netcdf
xu_grid.to_dataset().to_netcdf('test_net.nc')

#TODO: update https://github.com/Deltares/dfm_tools/issues/217 and https://github.com/Deltares/MeshKernelPy/issues/28

