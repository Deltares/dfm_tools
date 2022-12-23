# -*- coding: utf-8 -*-
"""
Created on Thu Dec  8 20:54:39 2022

@author: veenstra

"""

import meshkernel

import matplotlib.pyplot as plt
plt.close('all')
import numpy as np




## Make a uniform rectilinear grid
"""
Make a regular and rotated grid. This is possible in MeshKernelPy but I need to merge the related pull request in master. The same workflow of interactor is followed here: first generate a curvilinear grid than convert the curvilinear grid into unstructured grid. The steps are the following:
- curvilinear_make_uniform, see the following notebook: https://github.com/Deltares/MeshKernelPy/blob/AddCurvilinearGridSupport/docs/examples/04_curvilineargrid_basics.ipynb
- curvilinear_convert_to_mesh2d: https://github.com/Deltares/MeshKernelPy/blob/118cb4953c4e95d5b18ed283bb37f391134b2bb2/meshkernel/meshkernel.py#L1399 

"""

# Create an instance of MakeGridParameters and set the values
make_grid_parameters = meshkernel.MakeGridParameters()
make_grid_parameters.num_columns = 20
make_grid_parameters.num_rows = 20
make_grid_parameters.angle = 0.0
make_grid_parameters.block_size = 0.0
make_grid_parameters.origin_x = 0.0
make_grid_parameters.origin_y = 0.0
make_grid_parameters.block_size_x = 0.2
make_grid_parameters.block_size_y = 0.2



# A polygon must to be provided. If empty it will not be used.
node_x = np.empty(0, dtype=np.double)
node_y = np.empty(0, dtype=np.double)
geometry_list = meshkernel.GeometryList(node_x, node_y)

mk = meshkernel.MeshKernel()
mk.curvilinear_make_uniform(make_grid_parameters, geometry_list)
curvilinear_grid = mk.curvilineargrid_get()
fig, ax = plt.subplots()
curvilinear_grid.plot_edges(ax)
mk.curvilinear_convert_to_mesh2d() #convert to ugrid/mesh2d



#If a polygon is provided it will be used in the generation of the curvilinear grid. The polygon must be closed
node_x = np.array([2.5,5.5,3.5,0.5,2.5], dtype=np.double)
node_y = np.array([0.5,3.0,5.0,2.5,0.5], dtype=np.double)
geometry_list = meshkernel.GeometryList(node_x, node_y)

mk = meshkernel.MeshKernel()
mk.curvilinear_make_uniform(make_grid_parameters, geometry_list)
curvilinear_grid = mk.curvilineargrid_get()
fig, ax = plt.subplots()
curvilinear_grid.plot_edges(ax)
mk.curvilinear_convert_to_mesh2d() #convert to ugrid/mesh2d



# refine with samples. 
"""
Mesh refinement in MeshKernelPy (where the samples parameter contains the sample coordinates and the sample depths used for mesh refinement):
- https://github.com/Deltares/MeshKernelPy/blob/6a5110d75ea2a86a27f63da7c7196e008bbce3dc/meshkernel/meshkernel.py#L451


"""

mesh_refinement_parameters = meshkernel.MeshRefinementParameters(refine_intersected=True, #TODO: possible to make it not deliberate
                                                                 use_mass_center_when_refining=False,
                                                                 min_face_size=0.05,
                                                                 refinement_type=meshkernel.RefinementType(1), #Wavecourant/1,
                                                                 connect_hanging_nodes=False,
                                                                 account_for_samples_outside_face=False,
                                                                 max_refinement_iterations=3,
                                                                 )


samp_x = np.arange(0,5,0.1)
samp_y = np.arange(0,5,0.1)
samp_x,samp_y = np.meshgrid(samp_x,samp_y)
samp_z = np.random.uniform(low=0.5, high=13.3, size=samp_y.shape)

samp_x = samp_x.ravel()
samp_y = samp_y.ravel()
samp_z = samp_z.ravel()

geomlist = meshkernel.GeometryList(x_coordinates=samp_x, y_coordinates=samp_y, values=samp_z)

mk.mesh2d_refine_based_on_samples(samples=geomlist, #TODO: provide sample set (GeometryList type is expected)
                                  relative_search_radius=1,
                                  minimum_num_samples=3,
                                  mesh_refinement_params=mesh_refinement_parameters,
                                  )

mk.curvilinear_make_uniform(make_grid_parameters, geometry_list)
curvilinear_grid = mk.curvilineargrid_get()
fig, ax = plt.subplots()
curvilinear_grid.plot_edges(ax)

#TODO: write to netcdf with: xugrid or ugridpy (https://github.com/Deltares/UGridPy/blob/main/docs/examples/05_working_with_meshkernel.ipynb)




